//! TOF calibration window: the m/z calibration of FastFlight raw data. A FastFlight run is summed
//! per TOF bin; a click on a peak adds a reference ion (its TOF bin from the peak centre, its m/z
//! typed in), the polynomial is fitted through the ions and saved as a `.mzcalib` file, which
//! opening such a run needs. A saved file can be loaded to look at or correct its ions.

use std::path::PathBuf;
use std::sync::mpsc;

use eframe::egui;
use egui_plot::{Line, Plot, VLine};
use md_core::decimate;
use md_core::tof::{TofCalibration, TofError, peak_bin};

use crate::layout::ToolWindow;
use crate::{poll_job, spawn};

/// One reference ion: a peak of the run (TOF bin) and the m/z it stands for. Rows that are not
/// `used` stay in the table but out of the fit.
#[derive(Debug, Clone, PartialEq)]
struct Row {
    bin: f64,
    mz: f64,
    used: bool,
}

type Summed = (PathBuf, Vec<(f32, f64)>);

pub struct TofWindow {
    pub open: bool,
    /// The run being calibrated and its summed bin spectrum (bins and their intensity).
    run: Option<Summed>,
    xs: Vec<f64>,
    ys: Vec<f64>,
    job: Option<mpsc::Receiver<Result<Summed, String>>>,
    rows: Vec<Row>,
    degree: usize,
    /// Bins around a click searched for the peak.
    reach: f32,
    note: Option<String>,
    error: Option<String>,
}

impl Default for TofWindow {
    fn default() -> Self {
        Self {
            open: false,
            run: None,
            xs: Vec::new(),
            ys: Vec::new(),
            job: None,
            rows: Vec::new(),
            degree: 2,
            reach: 20.0,
            note: None,
            error: None,
        }
    }
}

impl TofWindow {
    pub fn window(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.job
            && let Some(result) = poll_job(rx, "FastFlight import").map(|r| r.and_then(|x| x))
        {
            match result {
                Ok(run) => self.set_run(run),
                Err(e) => self.error = Some(e),
            }
            self.job = None;
        }
        let mut open = self.open;
        egui::Window::new("TOF calibration")
            .tool_window(ctx, 11, [520.0, 560.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui));
        self.open = open;
    }

    fn set_run(&mut self, run: Summed) {
        (self.xs, self.ys) = run.1.iter().map(|(b, v)| (f64::from(*b), *v)).unzip();
        self.run = Some(run);
    }

    fn contents(&mut self, ui: &mut egui::Ui) {
        self.buttons(ui);
        if let Some(e) = &self.error {
            let mut clear = false;
            ui.horizontal_wrapped(|ui| {
                ui.colored_label(ui.visuals().error_fg_color, e);
                clear = ui.small_button("Clear").clicked();
            });
            if clear {
                self.error = None;
            }
        }
        if let Some(n) = &self.note {
            ui.label(n);
        }
        if self.job.is_some() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Adding up the scans of the run…");
            });
        }
        self.plot(ui);
        let fit = self.fit();
        self.table(ui, fit.as_ref().ok());
        self.result(ui, &fit);
    }

    fn buttons(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Open run…")
                .on_hover_text(
                    "Choose the MsData.FFT of a FastFlight run; its peaks are shown per TOF bin",
                )
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("FastFlight raw data (MsData.FFT)", &["FFT", "fft"])
                    .pick_file()
            {
                self.error = None;
                // The run is summed on its own thread: raw runs are large.
                let ctx = ui.ctx().clone();
                self.job = Some(spawn(&ctx, move || {
                    md_core::fastflight::bin_spectrum(&path)
                        .map(|s| (path, s))
                        .map_err(|e| e.to_string())
                }));
            }
            if ui
                .button("Load .mzcalib…")
                .on_hover_text(
                    "Read the reference ions of a saved calibration to check or correct them",
                )
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("TOF calibration", &["mzcalib"])
                    .pick_file()
            {
                self.load(&path);
            }
        });
        if let Some((path, _)) = &self.run {
            ui.label(format!("Run: {}", path.parent().unwrap_or(path).display()));
        }
    }

    fn load(&mut self, path: &std::path::Path) {
        let text = std::fs::read_to_string(path).map_err(|e| e.to_string());
        let loaded = text.and_then(|t| {
            let calibration = TofCalibration::parse(&t).map_err(|e| e.to_string())?;
            Ok((calibration, TofCalibration::parse_pairs(&t)))
        });
        match loaded {
            Ok((calibration, pairs)) => {
                self.degree = (calibration.coefficients.len() - 1).clamp(1, 6);
                self.rows = pairs
                    .into_iter()
                    .map(|(bin, mz)| Row {
                        bin: f64::from(bin),
                        mz,
                        used: true,
                    })
                    .collect();
                self.note = Some(match self.rows.is_empty() {
                    true => "The file lists no reference ions; add them from the run.".into(),
                    false => format!("{} reference ions loaded", self.rows.len()),
                });
                self.error = None;
            }
            Err(e) => self.error = Some(format!("{}: {e}", path.display())),
        }
    }

    /// The summed spectrum with a line at each reference ion; a click adds an ion at the peak there.
    fn plot(&mut self, ui: &mut egui::Ui) {
        if self.run.is_none() {
            ui.label(
                "Open a FastFlight run, then click its peaks to add reference ions. \
                 Or type TOF bins and m/z values into the table.",
            );
            return;
        }
        let (xs, ys, rows) = (&self.xs, &self.ys, &self.rows);
        let shown = Plot::new("tof bins")
            .height(180.0)
            .x_axis_label("TOF bin")
            .y_axis_label("Intensity (sum of all scans)")
            .y_axis_formatter(|m, _| crate::view::compact(m.value, m.step_size))
            .show(ui, |p| {
                let width = p.response().rect.width().max(100.0) as usize;
                let (lo, hi) = match p.auto_bounds().x {
                    true => (f64::NEG_INFINITY, f64::INFINITY),
                    false => {
                        let b = p.plot_bounds();
                        (b.min()[0], b.max()[0])
                    }
                };
                p.line(Line::new("run", decimate(xs, ys, lo, hi, width)));
                for r in rows.iter().filter(|r| r.used) {
                    p.vline(VLine::new("", r.bin));
                }
                p.pointer_coordinate()
            });
        if shown.response.clicked()
            && !shown.response.double_clicked()
            && let Some(at) = shown.inner
            && let Some((_, spectrum)) = &self.run
        {
            match peak_bin(spectrum, at.x as f32, self.reach) {
                Some(bin) => self.rows.push(Row {
                    bin: f64::from(bin),
                    mz: 0.0,
                    used: true,
                }),
                None => self.note = Some("No peak within reach of the click".into()),
            }
        }
        ui.horizontal(|ui| {
            let label = ui.label("Reach (bins)");
            ui.add(egui::DragValue::new(&mut self.reach).range(2.0..=500.0))
                .labelled_by(label.id)
                .on_hover_text("How far from a click the peak is searched, in TOF bins");
        });
    }

    fn fit(&self) -> Result<TofCalibration, TofError> {
        let pairs: Vec<(f32, f64)> = (self.rows.iter().filter(|r| r.used))
            .map(|r| (r.bin as f32, r.mz))
            .collect();
        TofCalibration::fit(&pairs, self.degree)
    }

    /// Editable reference ions with the distance of each from the fitted line.
    fn table(&mut self, ui: &mut egui::Ui, fit: Option<&TofCalibration>) {
        ui.separator();
        ui.horizontal(|ui| {
            ui.strong("Reference ions");
            if ui
                .small_button("Add ion")
                .on_hover_text("Add a row to fill in by hand")
                .clicked()
            {
                self.rows.push(Row {
                    bin: 0.0,
                    mz: 0.0,
                    used: true,
                });
            }
        });
        let mut remove = None;
        egui::Grid::new("tof ions").striped(true).show(ui, |ui| {
            for h in ["Use", "TOF bin", "m/z", "Error (mDa)", ""] {
                ui.label(h);
            }
            ui.end_row();
            for (i, r) in self.rows.iter_mut().enumerate() {
                ui.checkbox(&mut r.used, "")
                    .on_hover_text("Leave this ion out of the fit when unchecked");
                ui.add(egui::DragValue::new(&mut r.bin).speed(0.1).max_decimals(2))
                    .on_hover_text("TOF bin of the peak");
                ui.add(egui::DragValue::new(&mut r.mz).speed(0.01).max_decimals(4))
                    .on_hover_text("Mass of the reference ion");
                match fit {
                    Some(f) if r.used => {
                        let err = (f64::from(f.mz(r.bin as f32)) - r.mz) * 1000.0;
                        ui.label(format!("{err:+.1}"));
                    }
                    _ => {
                        ui.label("");
                    }
                }
                if ui
                    .small_button("×")
                    .on_hover_text("Remove this ion")
                    .clicked()
                {
                    remove = Some(i);
                }
                ui.end_row();
            }
        });
        if let Some(i) = remove {
            self.rows.remove(i);
        }
    }

    fn result(&mut self, ui: &mut egui::Ui, fit: &Result<TofCalibration, TofError>) {
        ui.separator();
        ui.horizontal(|ui| {
            let label = ui.label("Polynomial order");
            ui.add(egui::DragValue::new(&mut self.degree).range(1..=6))
                .labelled_by(label.id)
                .on_hover_text("Degree of the polynomial from TOF bin to m/z; needs one ion more");
        });
        match fit {
            Err(e) => {
                ui.label(e.to_string());
            }
            Ok(c) => {
                let terms: Vec<String> = c
                    .coefficients
                    .iter()
                    .enumerate()
                    .map(|(i, v)| format!("a{i} = {v:.6e}"))
                    .collect();
                ui.label(terms.join(",  "));
            }
        }
        let ready = fit.as_ref().ok();
        if ui
            .add_enabled(ready.is_some(), egui::Button::new("Save .mzcalib…"))
            .on_hover_text(
                "Write the calibration; open FastFlight runs from the same folder with it",
            )
            .on_disabled_hover_text("Add ions until the polynomial can be fitted")
            .clicked()
            && let Some(c) = ready
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("TOF calibration", &["mzcalib"])
                .set_file_name("calibration.mzcalib")
                .save_file()
        {
            let pairs: Vec<(f32, f64)> = (self.rows.iter().filter(|r| r.used))
                .map(|r| (r.bin as f32, r.mz))
                .collect();
            match c.write(&path, &pairs) {
                Ok(()) => self.note = Some(format!("Saved {}", path.display())),
                Err(e) => self.error = Some(e.to_string()),
            }
        }
    }
}

#[cfg(test)]
impl TofWindow {
    /// A summed run with two peaks and three reference ions, for the screenshot and overlap tests.
    pub fn demo() -> Self {
        let spectrum: Vec<(f32, f64)> = (0..3000)
            .map(|b| {
                let peak = |c: f64, h: f64| h * (-((f64::from(b) - c) / 6.0).powi(2) / 2.0).exp();
                (b as f32, 20.0 + peak(900.0, 5e5) + peak(2100.0, 2e5))
            })
            .collect();
        let mut w = Self {
            open: true,
            ..Self::default()
        };
        w.set_run(("MsData.FFT".into(), spectrum));
        w.rows = [(900.0, 73.0), (1500.0, 147.0), (2100.0, 281.0)]
            .map(|(bin, mz)| Row {
                bin,
                mz,
                used: true,
            })
            .to_vec();
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rows_fit_into_a_calibration() {
        let mut w = TofWindow::default();
        let c = TofCalibration {
            coefficients: vec![5.0, 0.01, 4e-8],
        };
        w.rows = [800.0, 1500.0, 2600.0, 4000.0]
            .map(|b| Row {
                bin: b,
                mz: f64::from(c.mz(b as f32)),
                used: true,
            })
            .to_vec();
        assert!(w.fit().is_ok());
        w.rows[3].used = false;
        w.rows[2].used = false;
        assert!(matches!(w.fit(), Err(TofError::NotEnoughPairs)));
    }
}
