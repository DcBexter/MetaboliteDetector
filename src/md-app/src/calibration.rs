//! RI calibration window: alkanes found in one run, or a `.calib` file; applies to every open run.

use std::sync::mpsc;

use eframe::egui;
use md_core::ri::{RiCalibration, marker_peaks, series_points};

use crate::{Run, spawn};

pub struct Calibration {
    pub cal: Option<RiCalibration>,
    /// Where `cal` came from, shown in the window.
    source: String,
    pub open: bool,
    mz: f32,
    /// Alkane carbon numbers: C`first` … C`last` get RI 100 · n.
    first: u32,
    last: u32,
    /// Matched detector settings of the old calibration dialog (no defaults there; these find
    /// the markers of 20 Hz example runs).
    density: usize,
    threshold: i32,
    /// Running alkane search and the source it will set.
    job: Option<(String, mpsc::Receiver<Result<RiCalibration, String>>)>,
    error: Option<String>,
}

impl Default for Calibration {
    fn default() -> Self {
        Self {
            cal: None,
            source: String::new(),
            open: false,
            // Alkane fragment C4H9+.
            mz: 57.0,
            first: 10,
            last: 36,
            density: 40,
            threshold: 100,
            job: None,
            error: None,
        }
    }
}

impl Calibration {
    /// ", RI …" for the scan, or nothing without a calibration.
    pub fn ri_note(&self, run: &Run, scan: usize) -> String {
        self.cal
            .as_ref()
            .map(|c| format!(", RI {:.1}", c.ri(run.data.scan_times_s[scan])))
            .unwrap_or_default()
    }

    fn set(&mut self, cal: RiCalibration, source: String) {
        (self.cal, self.source, self.error) = (Some(cal), source, None);
    }

    pub fn window(&mut self, ctx: &egui::Context, run: Option<&Run>) {
        if let Some((source, rx)) = &self.job
            && let Ok(result) = rx.try_recv()
        {
            match result {
                Ok(cal) => self.set(cal, source.clone()),
                Err(e) => self.error = Some(e),
            }
            self.job = None;
        }
        let mut open = self.open;
        egui::Window::new("RI calibration")
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui, run));
        self.open = open;
    }

    fn contents(&mut self, ui: &mut egui::Ui, run: Option<&Run>) {
        egui::Grid::new("alkane settings").show(ui, |ui| {
            ui.label("Marker ion m/z");
            ui.add(egui::DragValue::new(&mut self.mz).range(1.0..=2000.0));
            ui.end_row();
            ui.label("Alkanes C");
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.first).range(1..=self.last));
                ui.label("to C");
                ui.add(egui::DragValue::new(&mut self.last).range(self.first..=100));
            });
            ui.end_row();
            ui.label("Filter width (scans)");
            ui.add(egui::DragValue::new(&mut self.density).range(1..=1000));
            ui.end_row();
            ui.label("Threshold");
            ui.add(egui::DragValue::new(&mut self.threshold).range(0..=i32::MAX));
            ui.end_row();
        });
        ui.horizontal(|ui| {
            if let Some(run) = run
                && ui
                    .add_enabled(self.job.is_none(), egui::Button::new(format!("Find alkanes in {}", run.name)))
                    .on_hover_text("The highest peaks of the marker ion, in time order, become C first … C last")
                    .clicked()
            {
                let (data, mz, density, threshold) = (run.data.clone(), self.mz, self.density, self.threshold);
                let ris: Vec<f64> = (self.first..=self.last).map(|n| 100.0 * n as f64).collect();
                let source = format!("{}, m/z {mz}, C{}–C{}", run.name, self.first, self.last);
                let job = spawn(ui.ctx(), move || {
                    let peaks = marker_peaks(&data, mz, density, threshold);
                    RiCalibration::new(series_points(&peaks, &ris)).map_err(|e| e.to_string())
                });
                self.job = Some((source, job));
            }
            if ui.button("Load .calib…").clicked()
                && let Some(path) = calib_dialog().pick_file()
            {
                match RiCalibration::read_calib(&path) {
                    Ok(cal) => self.set(cal, path.display().to_string()),
                    Err(e) => self.error = Some(format!("{}: {e}", path.display())),
                }
            }
            if let Some(cal) = &self.cal {
                if ui.button("Save .calib…").clicked()
                    && let Some(path) = calib_dialog().save_file()
                    && let Err(e) = std::fs::write(&path, cal.to_calib())
                {
                    self.error = Some(format!("{}: {e}", path.display()));
                }
                if ui.button("Clear").clicked() {
                    self.cal = None;
                }
            }
        });
        if let Some(e) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, e);
        }
        let Some(cal) = &self.cal else {
            ui.label("Not calibrated: components have no RI.");
            return;
        };
        ui.label(format!(
            "{} points from {}",
            cal.points().len(),
            self.source
        ));
        egui::ScrollArea::vertical()
            .max_height(300.0)
            .show(ui, |ui| {
                egui::Grid::new("calibration points")
                    .striped(true)
                    .show(ui, |ui| {
                        ui.strong("RT (min)");
                        ui.strong("RI");
                        ui.end_row();
                        for (t, ri) in cal.points() {
                            ui.label(format!("{:.3}", t / 60.0));
                            ui.label(format!("{ri:.1}"));
                            ui.end_row();
                        }
                    });
            });
    }
}

fn calib_dialog() -> rfd::FileDialog {
    rfd::FileDialog::new().add_filter("RI calibration", &["calib"])
}
