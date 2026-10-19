//! MetaboliteDetector desktop app (egui). Analysis lives in `md-core`; this crate only draws.

use std::path::PathBuf;
use std::sync::mpsc;

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, Points, VLine};
use md_core::deconv::{Component, DeconvSettings, deconvolute};
use md_core::peaks::{IonTraces, PeakSettings};
use md_core::{Chromatogram, decimate};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("MetaboliteDetector")
            .with_inner_size([1200.0, 800.0]),
        ..Default::default()
    };
    eframe::run_native(
        "MetaboliteDetector",
        options,
        Box::new(|cc| {
            let mut app = App::default();
            if let Some(path) = std::env::args_os().nth(1) {
                app.open(&cc.egui_ctx, path.into());
            }
            Ok(Box::new(app))
        }),
    )
}

#[derive(Default)]
struct App {
    run: Option<Run>,
    loading: Option<(PathBuf, mpsc::Receiver<Result<Run, String>>)>,
    error: Option<String>,
}

struct Run {
    path: PathBuf,
    data: Chromatogram,
    /// TIC x (minutes) and y, cached once per file.
    tic_min: Vec<f64>,
    tic: Vec<f64>,
    /// Deconvoluted substances with the old default settings, ordered by retention time.
    components: Vec<Component>,
    scan: usize,
    /// Selected component; its deconvoluted spectrum replaces the raw one.
    component: Option<usize>,
}

impl App {
    /// Parses on a worker thread so huge files never freeze the window.
    fn open(&mut self, ctx: &egui::Context, path: PathBuf) {
        let (tx, rx) = mpsc::channel();
        let (ctx, worker_path) = (ctx.clone(), path.clone());
        std::thread::spawn(move || {
            let result = md_core::andi::read(&worker_path)
                .map(|data| {
                    let traces = IonTraces::new(&data);
                    let peaks = traces.detect_peaks(&PeakSettings::default());
                    let mut components = deconvolute(&traces, &peaks, &DeconvSettings::default());
                    components.sort_by_key(|k| k.apex_scan);
                    Run {
                        tic_min: data.scan_times_s.iter().map(|t| t / 60.0).collect(),
                        tic: data.tic(),
                        components,
                        data,
                        path: worker_path.clone(),
                        scan: 0,
                        component: None,
                    }
                })
                .map_err(|e| format!("{}: {e}", worker_path.display()));
            tx.send(result).ok();
            ctx.request_repaint();
        });
        self.loading = Some((path, rx));
    }

    fn poll_loading(&mut self) {
        let Some((_, rx)) = &self.loading else { return };
        match rx.try_recv() {
            Ok(Ok(run)) => (self.run, self.error) = (Some(run), None),
            Ok(Err(e)) => self.error = Some(e),
            Err(mpsc::TryRecvError::Empty) => return,
            Err(mpsc::TryRecvError::Disconnected) => self.error = Some("loader crashed".into()),
        }
        self.loading = None;
    }
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.poll_loading();
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open NetCDF (ANDI-MS)…").clicked()
                        && let Some(path) = rfd::FileDialog::new()
                            .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                            .pick_file()
                    {
                        self.open(ui.ctx(), path);
                    }
                });
            });
        });

        egui::Panel::bottom("status").show(ui, |ui| {
            match (&self.loading, &self.error, &self.run) {
                (Some((path, _)), _, _) => {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label(format!("Loading {}…", path.display()))
                    })
                    .inner
                }
                (None, Some(e), _) => ui.colored_label(ui.visuals().error_fg_color, e),
                (None, None, Some(r)) => {
                    let t = r.data.scan_times_s[r.scan];
                    let shown = match r.component {
                        Some(k) => format!(
                            "component {} of {} ({} ions)",
                            k + 1,
                            r.components.len(),
                            r.components[k].spectrum.len()
                        ),
                        None => format!("{} components (click a dot)", r.components.len()),
                    };
                    ui.label(format!(
                        "{} — {} scans — scan {} at {:.3} min (←/→ to step) — {shown}",
                        r.path.display(),
                        r.data.scan_count(),
                        r.scan + 1,
                        t / 60.0
                    ))
                }
                (None, None, None) => {
                    ui.label("Open an ANDI-MS NetCDF chromatogram via File → Open.")
                }
            }
        });

        egui::CentralPanel::default().show(ui, |ui| {
            let Some(run) = self.run.as_mut().filter(|r| r.data.scan_count() > 0) else {
                return;
            };
            let last = run.data.scan_count() - 1;
            ui.input(|i| {
                if i.key_pressed(egui::Key::ArrowRight) {
                    (run.scan, run.component) = ((run.scan + 1).min(last), None);
                }
                if i.key_pressed(egui::Key::ArrowLeft) {
                    (run.scan, run.component) = (run.scan.saturating_sub(1), None);
                }
            });

            let half = ui.available_height() / 2.0;
            let selected_min = run.data.scan_times_s[run.scan] / 60.0;
            let tic = Plot::new("tic")
                .height(half)
                .x_axis_label("Retention time (min)")
                .y_axis_label("TIC")
                .show(ui, |p| {
                    // Draw at most ~2 points per pixel of the visible range, whatever the file size.
                    let (lo, hi) = if p.auto_bounds().x {
                        (f64::NEG_INFINITY, f64::INFINITY)
                    } else {
                        (p.plot_bounds().min()[0], p.plot_bounds().max()[0])
                    };
                    let width = p.response().rect.width().max(100.0) as usize;
                    p.line(Line::new(
                        "TIC",
                        decimate(&run.tic_min, &run.tic, lo, hi, width),
                    ));
                    let apex = |k: &Component| {
                        let s = k.apex_scan as usize;
                        [run.tic_min[s], run.tic[s]]
                    };
                    p.points(
                        Points::new(
                            "components",
                            run.components.iter().map(apex).collect::<Vec<_>>(),
                        )
                        .radius(3.0),
                    );
                    p.vline(VLine::new("selected scan", selected_min));
                    p.pointer_coordinate()
                });
            if tic.response.clicked()
                && let Some(s) = tic.inner.and_then(|pt| run.data.scan_at(pt.x * 60.0))
            {
                // A click within 3 scans of a component apex selects that component.
                run.component = (0..run.components.len())
                    .min_by_key(|k| run.components[*k].apex_scan.abs_diff(s as u32))
                    .filter(|k| run.components[*k].apex_scan.abs_diff(s as u32) <= 3);
                run.scan = run
                    .component
                    .map_or(s, |k| run.components[k].apex_scan as usize);
            }

            let (name, bars) = match run.component {
                Some(k) => (
                    "deconvoluted spectrum",
                    run.components[k]
                        .spectrum
                        .iter()
                        .map(|(m, h)| Bar::new(*m as f64, *h as f64).width(0.3))
                        .collect(),
                ),
                None => {
                    let (mz, intensity) = run.data.spectrum(run.scan);
                    let bars = mz
                        .iter()
                        .zip(intensity)
                        .map(|(m, i)| Bar::new(*m as f64, *i as f64).width(0.3))
                        .collect();
                    ("spectrum", bars)
                }
            };
            Plot::new("spectrum")
                .x_axis_label("m/z")
                .y_axis_label("Intensity")
                .show(ui, |p| p.bar_chart(BarChart::new(name, bars)));
        });
    }
}
