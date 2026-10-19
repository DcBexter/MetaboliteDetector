//! MetaboliteDetector desktop app (egui). Analysis lives in `md-core`; this crate only draws.

use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};
use std::time::Instant;

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, Points, Text, VLine};
use md_core::deconv::{Component, DeconvSettings, deconvolute};
use md_core::peaks::{IonTraces, PeakSettings};
use md_core::{Chromatogram, decimate};

mod calibration;
mod identify;

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
            for path in std::env::args_os().skip(1).map(PathBuf::from) {
                if identify::is_library(&path) {
                    app.open_library(&cc.egui_ctx, path);
                } else {
                    app.open(&cc.egui_ctx, path);
                }
            }
            Ok(Box::new(app))
        }),
    )
}

#[derive(Clone, Copy, PartialEq)]
enum Trace {
    Tic,
    Sic,
    Both,
}

struct App {
    runs: Vec<Run>,
    /// Index into `runs` whose scan cursor and spectrum are shown.
    active: usize,
    loading: Vec<(PathBuf, mpsc::Receiver<Result<Run, String>>)>,
    load_started: Option<Instant>,
    load_note: String,
    error: Option<String>,
    trace: Trace,
    mz: f32,
    tolerance: f32,
    ri: calibration::Calibration,
    library: Option<Arc<identify::Library>>,
    library_job: Option<identify::LibraryJob>,
}

impl Default for App {
    fn default() -> Self {
        Self {
            runs: Vec::new(),
            active: 0,
            loading: Vec::new(),
            load_started: None,
            load_note: String::new(),
            error: None,
            trace: Trace::Tic,
            // TMS fragment; many runs start scanning above m/z 73.
            mz: 147.0,
            tolerance: 0.5,
            ri: calibration::Calibration::default(),
            library: None,
            library_job: None,
        }
    }
}

/// (m/z, tolerance) of a SIC.
type SicKey = (f32, f32);

struct Run {
    path: PathBuf,
    name: String,
    data: Arc<Chromatogram>,
    /// Scan times in minutes and TIC, cached once per file.
    minutes: Vec<f64>,
    tic: Vec<f64>,
    /// SIC for the (m/z, tolerance) in `sic_for`; recomputed on a worker thread when that changes.
    sic: Vec<f64>,
    sic_for: Option<SicKey>,
    sic_job: Option<(SicKey, mpsc::Receiver<Vec<f64>>)>,
    /// Deconvoluted substances with the old default settings, ordered by retention time.
    components: Vec<Component>,
    scan: usize,
    /// Selected component; its deconvoluted spectrum replaces the raw one.
    component: Option<usize>,
    visible: bool,
    /// Library hits of `components` once a library is open.
    ident: Option<identify::Identified>,
    ident_job: Option<identify::IdentJob>,
}

/// Runs `job` on its own thread and wakes the UI when it is done, so the window never freezes.
fn spawn<T: Send + 'static>(
    ctx: &egui::Context,
    job: impl FnOnce() -> T + Send + 'static,
) -> mpsc::Receiver<T> {
    let (tx, rx) = mpsc::channel();
    let ctx = ctx.clone();
    std::thread::spawn(move || {
        tx.send(job()).ok();
        ctx.request_repaint();
    });
    rx
}

/// Old MetaboliteDetector caches (`.bin` + `.idx`) by extension, everything else as ANDI-MS NetCDF.
fn read_chromatogram(p: &Path) -> Result<md_core::Chromatogram, String> {
    match p.extension().is_some_and(|e| e.eq_ignore_ascii_case("bin")) {
        true => md_core::diskscan::read(p).map_err(|e| e.to_string()),
        false => md_core::andi::read(p).map_err(|e| e.to_string()),
    }
}

/// Same colour per run in every plot and in the run list (egui_plot's golden-ratio hues).
fn color(i: usize) -> egui::Color32 {
    egui::ecolor::Hsva::new(i as f32 * 0.618_034, 0.85, 0.5, 1.0).into()
}

impl App {
    /// Every file loads on its own thread, so many runs load in parallel.
    fn open(&mut self, ctx: &egui::Context, path: PathBuf) {
        let p = path.clone();
        let rx = spawn(ctx, move || {
            read_chromatogram(&p)
                .and_then(|data| match data.scan_count() {
                    0 => Err("no scans".into()),
                    _ => {
                        let traces = IonTraces::new(&data);
                        let peaks = traces.detect_peaks(&PeakSettings::default());
                        let mut components =
                            deconvolute(&traces, &peaks, &DeconvSettings::default());
                        components.sort_by_key(|k| k.apex_scan);
                        Ok(Run {
                            name: p.file_stem().unwrap_or_default().to_string_lossy().into(),
                            minutes: data.scan_times_s.iter().map(|t| t / 60.0).collect(),
                            tic: data.tic(),
                            data: Arc::new(data),
                            path: p.clone(),
                            sic: Vec::new(),
                            sic_for: None,
                            sic_job: None,
                            components,
                            scan: 0,
                            component: None,
                            visible: true,
                            ident: None,
                            ident_job: None,
                        })
                    }
                })
                .map_err(|e| format!("{}: {e}", p.display()))
        });
        self.load_started.get_or_insert_with(Instant::now);
        self.loading.push((path, rx));
    }

    fn poll_loading(&mut self) {
        let (count, active) = (
            self.runs.len(),
            self.runs.get(self.active).map(|r| r.path.clone()),
        );
        self.loading.retain(|(_, rx)| {
            match rx.try_recv() {
                Ok(Ok(run)) => self.runs.push(run),
                Ok(Err(e)) => self.error = Some(e),
                Err(mpsc::TryRecvError::Empty) => return true,
                Err(mpsc::TryRecvError::Disconnected) => self.error = Some("loader crashed".into()),
            }
            false
        });
        // Runs finish in any order; sorting by path keeps list order and colours stable.
        if self.runs.len() != count {
            self.runs.sort_by(|a, b| a.path.cmp(&b.path));
            self.active = active
                .and_then(|p| self.runs.iter().position(|r| r.path == p))
                .unwrap_or(0);
        }
        if self.loading.is_empty()
            && let Some(t) = self.load_started.take()
        {
            self.load_note = format!(
                "{} runs open, last batch loaded in {:.2?}",
                self.runs.len(),
                t.elapsed()
            );
        }
    }

    /// Collects finished SICs and starts one job per run whose SIC is stale.
    fn update_sics(&mut self, ctx: &egui::Context) {
        let key = (self.mz, self.tolerance);
        for run in &mut self.runs {
            if let Some((k, rx)) = &run.sic_job
                && let Ok(sic) = rx.try_recv()
            {
                (run.sic, run.sic_for) = (sic, Some(*k));
                run.sic_job = None;
            }
            if run.sic_for != Some(key) && run.sic_job.as_ref().map(|j| j.0) != Some(key) {
                let data = run.data.clone();
                run.sic_job = Some((key, spawn(ctx, move || data.sic(key.0, key.1))));
            }
        }
    }
}

/// One chromatogram plot with a line per run and `marks` as dots; returns the clicked retention
/// time (min).
#[allow(clippy::too_many_arguments)]
fn trace_plot(
    ui: &mut egui::Ui,
    id: &str,
    label: &str,
    height: f32,
    runs: &[Run],
    cursor_min: f64,
    trace: fn(&Run) -> &[f64],
    marks: &[[f64; 2]],
    labels: &[identify::Label],
) -> Option<f64> {
    let link = egui::Id::unique("chromatograms");
    let plot = Plot::new(id)
        .height(height)
        // Same width in both plots, so their linked x axes line up.
        .y_axis_min_width(90.0)
        .link_axis(link, [true, false])
        .link_cursor(link, [true, false])
        .x_axis_label("Retention time (min)")
        .y_axis_label(label)
        .show(ui, |p| {
            // Draw at most ~2 points per pixel of the visible range, whatever the file size.
            let (lo, hi) = if p.auto_bounds().x {
                (f64::NEG_INFINITY, f64::INFINITY)
            } else {
                (p.plot_bounds().min()[0], p.plot_bounds().max()[0])
            };
            let width = p.response().rect.width().max(100.0) as usize;
            for (i, run) in runs.iter().enumerate() {
                let ys = trace(run);
                // A SIC is empty until its first job finishes.
                if run.visible && ys.len() == run.minutes.len() {
                    let points = decimate(&run.minutes, ys, lo, hi, width);
                    p.line(Line::new(run.name.as_str(), points).color(color(i)));
                }
            }
            if !marks.is_empty() {
                p.points(Points::new("components", marks.to_vec()).radius(3.0));
            }
            for (at, name, color) in labels {
                let text = egui::RichText::new(name).color(*color);
                p.text(Text::new("hits", (*at).into(), text).anchor(egui::Align2::LEFT_BOTTOM));
            }
            p.vline(VLine::new("", cursor_min));
            p.pointer_coordinate()
        });
    plot.response.clicked().then_some(plot.inner?.x)
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.poll_loading();
        self.update_identification(ui.ctx());
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| {
                    if ui.button("Open chromatograms…").clicked()
                        && let Some(paths) = rfd::FileDialog::new()
                            .add_filter("NetCDF (ANDI-MS)", &["cdf", "CDF", "nc"])
                            .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
                            .pick_files()
                    {
                        for path in paths {
                            self.open(ui.ctx(), path);
                        }
                    }
                    if let Some(r) = self.runs.get(self.active)
                        && ui.button("Export as NetCDF…").clicked()
                        && let Some(path) = rfd::FileDialog::new()
                            .add_filter("NetCDF (ANDI-MS)", &["cdf"])
                            .set_file_name(format!("{}.cdf", r.name))
                            .save_file()
                    {
                        // shortcut: writes on the UI thread (~50 ms per example run); spawn if runs grow.
                        self.error = md_core::andi::write(&path, &r.data)
                            .err()
                            .map(|e| format!("{}: {e}", path.display()));
                    }
                    self.library_menu(ui);
                    if ui.button("Close all").clicked() {
                        (self.runs, self.active) = (Vec::new(), 0);
                    }
                });
                ui.menu_button("Analysis", |ui| {
                    if ui.button("RI calibration…").clicked() {
                        self.ri.open = true;
                    }
                });
            });
        });

        self.ri.window(ui.ctx(), self.runs.get(self.active));

        egui::Panel::bottom("status").show(ui, |ui| {
            if !self.loading.is_empty() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Loading {} file(s)…", self.loading.len()));
                });
            } else if let Some(e) = &self.error {
                ui.colored_label(ui.visuals().error_fg_color, e);
            } else if let Some(r) = self.runs.get(self.active) {
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
                    "{} — {} scans — scan {} at {:.3} min{} (←/→ to step) — {shown} — {}",
                    r.path.display(),
                    r.data.scan_count(),
                    r.scan + 1,
                    r.minutes[r.scan],
                    self.ri.ri_note(r, r.scan),
                    self.load_note
                ));
            } else {
                ui.label("Open one or more chromatograms (NetCDF or old .bin) via File → Open.");
            }
        });

        if !self.runs.is_empty() {
            egui::Panel::left("runs").show(ui, |ui| {
                let mut close = None;
                for (i, run) in self.runs.iter_mut().enumerate() {
                    ui.horizontal(|ui| {
                        if ui.small_button("×").on_hover_text("Close").clicked() {
                            close = Some(i);
                        }
                        ui.checkbox(&mut run.visible, "")
                            .on_hover_text("Show in plots");
                        ui.colored_label(color(i), "■");
                        if ui
                            .selectable_label(i == self.active, &run.name)
                            .on_hover_text(run.path.display().to_string())
                            .clicked()
                        {
                            self.active = i;
                        }
                    });
                }
                if let Some(i) = close {
                    self.runs.remove(i);
                    if self.active > i || self.active >= self.runs.len() {
                        self.active = self.active.saturating_sub(1);
                    }
                }
            });
        }

        if let Some(run) = self.runs.get(self.active) {
            identify::hits_panel(ui, run);
        }

        egui::CentralPanel::default().show(ui, |ui| {
            ui.horizontal(|ui| {
                ui.selectable_value(&mut self.trace, Trace::Tic, "TIC");
                ui.selectable_value(&mut self.trace, Trace::Sic, "SIC");
                ui.selectable_value(&mut self.trace, Trace::Both, "TIC + SIC");
                ui.add_enabled_ui(self.trace != Trace::Tic, |ui| {
                    ui.label("m/z");
                    ui.add(
                        egui::DragValue::new(&mut self.mz)
                            .range(1.0..=2000.0)
                            .max_decimals(2),
                    );
                    ui.label("±");
                    ui.add(
                        egui::DragValue::new(&mut self.tolerance)
                            .range(0.0..=10.0)
                            .speed(0.05)
                            .max_decimals(2),
                    );
                });
            });
            if self.trace != Trace::Tic {
                self.update_sics(ui.ctx());
            }
            let Some(run) = self.runs.get_mut(self.active) else {
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
            let cursor = run.minutes[run.scan];

            // Both traces: a third each, so their axes keep enough room for tick labels.
            let plots = if self.trace == Trace::Both { 3.0 } else { 2.0 };
            let height = ui.available_height() / plots;
            let mut clicked = None;
            if self.trace != Trace::Sic {
                // Components of the active run as dots on its TIC.
                let marks: Vec<[f64; 2]> = run
                    .components
                    .iter()
                    .map(|k| {
                        [
                            run.minutes[k.apex_scan as usize],
                            run.tic[k.apex_scan as usize],
                        ]
                    })
                    .collect();
                let labels = identify::labels(run);
                clicked = trace_plot(
                    ui,
                    "tic",
                    "TIC",
                    height,
                    &self.runs,
                    cursor,
                    |r| &r.tic,
                    &marks,
                    &labels,
                );
            }
            if self.trace != Trace::Tic {
                let label = format!("SIC m/z {} ± {}", self.mz, self.tolerance);
                let sic = trace_plot(
                    ui,
                    "sic",
                    &label,
                    height,
                    &self.runs,
                    cursor,
                    |r| &r.sic,
                    &[],
                    &[],
                );
                clicked = clicked.or(sic);
            }

            let run = &mut self.runs[self.active];
            if let Some(s) = clicked.and_then(|x| run.data.scan_at(x * 60.0)) {
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
                .y_axis_label(format!("Intensity — {}", run.name))
                .show(ui, |p| p.bar_chart(BarChart::new(name, bars)));
        });
    }
}
