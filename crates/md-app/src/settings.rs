//! Settings of the old version (`AN_*`, `LS_*`, `QUANT_*`, …): a dialog with the old defaults, load and
//! save of the old settings file, and the background re-analysis of open runs when they change.

use std::io::Read;
use std::ops::RangeInclusive;
use std::sync::mpsc;

use eframe::egui;
use md_core::Chromatogram;
use md_core::deconv::{Component, DeconvSettings, deconvolute};
use md_core::peaks::{IonTraces, Peak, PeakSettings};
use md_core::settings::{GcmsSettings, MAX_LEN};

use crate::{Run, spawn};

/// What the components of a run were computed with.
pub type Key = (PeakSettings, DeconvSettings);

#[derive(Default)]
pub struct Settings {
    /// In use: new runs open with it and open runs are re-analysed to match it.
    pub applied: GcmsSettings,
    /// Being edited in the dialog until "Apply".
    draft: GcmsSettings,
    pub open: bool,
    note: Option<Result<String, String>>,
}

type Analyzed = Result<(Vec<Component>, Vec<u32>, Vec<Peak>), String>;

/// Re-analysis state of one run.
#[derive(Default)]
pub struct Analysis {
    /// Settings behind `Run::components`.
    pub key: Option<Key>,
    job: Option<(Key, mpsc::Receiver<Analyzed>)>,
}

impl Analysis {
    pub fn new(key: Key) -> Self {
        Self {
            key: Some(key),
            job: None,
        }
    }

    pub fn busy(&self) -> bool {
        self.job.is_some()
    }
}

/// Peak detection and deconvolution of one run, components ordered by retention time.
pub fn analyze(data: &Chromatogram, key: &Key) -> Analyzed {
    let traces = IonTraces::new(data).map_err(|e| e.to_string())?;
    let peaks = traces.detect_peaks(&key.0);
    let mut components = deconvolute(&traces, &peaks, &key.1);
    components.sort_by_key(|k| k.apex_scan);
    let masses = peaks
        .iter()
        .map(|p| traces.masses[p.row as usize])
        .collect();
    Ok((components, masses, peaks))
}

impl Settings {
    pub fn key(&self) -> Key {
        (self.applied.peak(), self.applied.deconv())
    }

    /// Collects finished analyses and starts one for every run whose components are stale.
    /// Returns the error of an analysis that failed (the run keeps its old components).
    pub fn update_runs(&self, ctx: &egui::Context, runs: &mut [Run]) -> Option<String> {
        let wanted = self.key();
        let mut error = None;
        for run in runs {
            let a = &mut run.analysis;
            if let Some((key, rx)) = &a.job {
                match rx.try_recv() {
                    Ok(result) => {
                        match result {
                            Ok((components, masses, peaks)) => {
                                (run.components, run.peak_masses, run.peaks) =
                                    (components, masses, peaks);
                                // Selection and library hits belong to the old components.
                                (run.component, run.ident, run.ident_job) = (None, None, None);
                                run.rev += 1;
                            }
                            Err(e) => error = Some(format!("{}: {e}", run.name)),
                        }
                        a.key = Some(key.clone());
                        a.job = None;
                    }
                    Err(mpsc::TryRecvError::Empty) => {}
                    Err(mpsc::TryRecvError::Disconnected) => {
                        // Panicked: keep the old components and do not retry these settings.
                        error = Some(format!("{}: analysis crashed", run.name));
                        a.key = Some(key.clone());
                        a.job = None;
                    }
                }
            }
            if a.key.as_ref() == Some(&wanted) {
                // Settings changed back while a job for other ones ran: its result is not wanted.
                a.job = None;
            } else if a.job.as_ref().map(|j| &j.0) != Some(&wanted) {
                let (data, key) = (run.data.clone(), wanted.clone());
                a.job = Some((wanted.clone(), spawn(ctx, move || analyze(&data, &key))));
            }
        }
        error
    }

    pub fn window(&mut self, ctx: &egui::Context) {
        let mut open = self.open;
        egui::Window::new("Settings")
            .open(&mut open)
            .default_width(480.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical()
                    .max_height(ui.available_height() - 60.0)
                    .show(ui, |ui| self.fields(ui));
                ui.separator();
                self.buttons(ui);
            });
        self.open = open;
    }

    fn buttons(&mut self, ui: &mut egui::Ui) {
        let valid = self.draft.validate();
        ui.horizontal(|ui| {
            let changed = self.draft != self.applied;
            if ui
                .add_enabled(valid.is_ok() && changed, egui::Button::new("Apply"))
                .on_hover_text("Re-runs peak detection, deconvolution and the library search")
                .clicked()
            {
                self.applied = self.draft.clone();
                self.note = None;
            }
            if ui.button("Old defaults").clicked() {
                self.draft = GcmsSettings::default();
            }
            if ui.button("Load…").clicked()
                && let Some(path) = dialog().pick_file()
            {
                self.note = Some(self.load(&path));
            }
            if ui
                .add_enabled(valid.is_ok(), egui::Button::new("Save…"))
                .clicked()
                && let Some(path) = dialog().set_file_name("GCMSSettings").save_file()
            {
                self.note = Some(
                    std::fs::write(&path, self.draft.to_bytes())
                        .map(|()| format!("Saved {}", path.display()))
                        .map_err(|e| format!("{}: {e}", path.display())),
                );
            }
        });
        match (&valid, &self.note) {
            (Err(e), _) | (_, Some(Err(e))) => {
                ui.colored_label(ui.visuals().error_fg_color, e);
            }
            (_, Some(Ok(m))) => {
                ui.label(m);
            }
            _ => {}
        }
    }

    /// Reads an old settings file into the dialog and applies it. The file is untrusted:
    /// `from_bytes` checks size, checksum and values, and reading stops at `MAX_LEN`.
    fn load(&mut self, path: &std::path::Path) -> Result<String, String> {
        let err = |e: &dyn std::fmt::Display| format!("{}: {e}", path.display());
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .and_then(|f| f.take(MAX_LEN + 1).read_to_end(&mut bytes))
            .map_err(|e| err(&e))?;
        if bytes.len() as u64 > MAX_LEN {
            return Err(err(&"file too large for a settings file"));
        }
        self.draft = self.applied.from_bytes(&bytes).map_err(|e| err(&e))?;
        self.applied = self.draft.clone();
        Ok(format!("Loaded {}", path.display()))
    }

    fn fields(&mut self, ui: &mut egui::Ui) {
        let d = &mut self.draft;
        ui.strong("Peak detection and deconvolution (AN_*)");
        egui::Grid::new("an").num_columns(2).show(ui, |ui| {
            f64_row(
                ui,
                "AN_PEAK_THRESHOLD_BEGIN",
                "Smoothed slope that starts a peak.",
                &mut d.an_peak_threshold_begin,
            );
            f64_row(
                ui,
                "AN_PEAK_THRESHOLD_END",
                "Falling flank ends once the slope stays above this.",
                &mut d.an_peak_threshold_end,
            );
            range_row(
                ui,
                "AN_MIN_PEAK_HEIGHT",
                "Minimum height above the baseline in noise units (S/N).",
                &mut d.an_min_peak_height,
                0.0..=1e9,
            );
            i32_row(
                ui,
                "AN_NUMBER_BINS",
                "Sub-scan steps of the precise apex.",
                &mut d.an_number_bins,
                1..=100,
            );
            range_row(
                ui,
                "AN_DECONVOLUTION_WIDTH",
                "Width in scans of the matched filter over the apex histogram.",
                &mut d.an_deconvolution_width,
                0.1..=100.0,
            );
            i32_row(
                ui,
                "AN_MIN_PEAK_NUMBER",
                "Minimum ion peaks per component.",
                &mut d.an_min_peak_number,
                1..=10_000,
            );
            range_row(
                ui,
                "AN_MIN_HEIGHT_BP",
                "Ions below this fraction of the base peak leave the spectrum (0 = keep all).",
                &mut d.an_min_height_bp,
                0.0..=1.0,
            );
            bool_row(
                ui,
                "AN_MODELED_VALUES_FOR_INTEGRATION",
                "Integrate the fitted model instead of the raw trace.",
                &mut d.an_modeled_values_for_integration,
            );
        });
        ranges(
            ui,
            "QUANT_EXCL_IONS",
            "Ion ranges (m/z from–to) not used for the model peak.",
            &mut d.quant_excl_ions,
        );

        ui.add_space(8.0);
        ui.strong("Library search (LS_*)");
        egui::Grid::new("ls").num_columns(2).show(ui, |ui| {
            range_row(
                ui,
                "LS_RI_DIFF",
                "Sigma of the RI score; differences above 3 sigma score 0.",
                &mut d.ls_ri_diff,
                1e-6..=1e9,
            );
            range_row(
                ui,
                "LS_PURE_FACTOR",
                "Weight of the pure spectrum score.",
                &mut d.ls_pure_factor,
                0.0..=1e9,
            );
            range_row(
                ui,
                "LS_IMPURE_FACTOR",
                "Weight of the impure spectrum score.",
                &mut d.ls_impure_factor,
                0.0..=1e9,
            );
            f64_row(
                ui,
                "LS_SCALE_FACTOR_ION",
                "Peaks are weighted as mass^a · intensity^b: a.",
                &mut d.ls_scale_factor_ion,
            );
            f64_row(
                ui,
                "LS_SCALE_FACTOR_VALUE",
                "… and b.",
                &mut d.ls_scale_factor_value,
            );
            f64_row(
                ui,
                "LS_SCALE_FACTOR_NO_MATCH",
                "Factor for peaks without a partner in the other spectrum.",
                &mut d.ls_scale_factor_no_match,
            );
            bool_row(
                ui,
                "LS_USE_RI",
                "RI score counts in the overall score.",
                &mut d.ls_use_ri,
            );
            bool_row(
                ui,
                "LS_USE_SPEC",
                "Spectrum score counts in the overall score.",
                &mut d.ls_use_spec,
            );
        });
        ranges(
            ui,
            "LS_MASS_FILTER",
            "Mass ranges (m/z from–to) ignored in both spectra.",
            &mut d.ls_mass_filter,
        );

        ui.add_space(8.0);
        egui::CollapsingHeader::new("Stored in the file, not used yet").show(ui, |ui| {
            egui::Grid::new("stored").num_columns(2).show(ui, |ui| {
                f64_row(
                    ui,
                    "SF_PEAK_THRESHOLD_BEGIN",
                    "",
                    &mut d.sf_peak_threshold_begin,
                );
                f64_row(
                    ui,
                    "SF_PEAK_THRESHOLD_END",
                    "",
                    &mut d.sf_peak_threshold_end,
                );
                f64_row(ui, "SF_MAX_BASE_DEV", "", &mut d.sf_max_base_dev);
                f64_row(ui, "SF_FWHM", "", &mut d.sf_fwhm);
                f64_row(
                    ui,
                    "SF_CENTROID_THRESHOLD_BEGIN",
                    "",
                    &mut d.sf_centroid_threshold_begin,
                );
                f64_row(
                    ui,
                    "SF_CENTROID_THRESHOLD_END",
                    "",
                    &mut d.sf_centroid_threshold_end,
                );
                f64_row(
                    ui,
                    "SF_CENTROID_MAX_BASE_DEV",
                    "",
                    &mut d.sf_centroid_max_base_dev,
                );
                f64_row(
                    ui,
                    "AN_MAX_BASE_DEV",
                    "The old peak detector does not read it.",
                    &mut d.an_max_base_dev,
                );
                let mut res = d.ls_resolution as f64;
                f64_row(ui, "LS_RESOLUTION", "", &mut res);
                d.ls_resolution = res as f32;
                bool_row(ui, "LS_USE_SCALED_LIB", "", &mut d.ls_use_scaled_lib);
                i32_row(
                    ui,
                    "QUANT_ION_NUMBER",
                    "",
                    &mut d.quant_ion_number,
                    0..=10_000,
                );
                f64_row(ui, "QUANT_MIN_DIST_IONS", "", &mut d.quant_min_dist_ions);
                f64_row(ui, "QUANT_QUALITY_CUTOFF", "", &mut d.quant_quality_cutoff);
            });
        });
    }
}

/// The old file has no extension convention.
fn dialog() -> rfd::FileDialog {
    rfd::FileDialog::new()
}

fn f64_row(ui: &mut egui::Ui, name: &str, hint: &str, v: &mut f64) {
    range_row(ui, name, hint, v, -1e9..=1e9);
}

fn range_row(ui: &mut egui::Ui, name: &str, hint: &str, v: &mut f64, range: RangeInclusive<f64>) {
    ui.label(name).on_hover_text(hint);
    ui.add(
        egui::DragValue::new(v)
            .range(range)
            .speed(0.1)
            .max_decimals(6),
    );
    ui.end_row();
}

fn i32_row(ui: &mut egui::Ui, name: &str, hint: &str, v: &mut i32, range: RangeInclusive<i32>) {
    ui.label(name).on_hover_text(hint);
    ui.add(egui::DragValue::new(v).range(range));
    ui.end_row();
}

fn bool_row(ui: &mut egui::Ui, name: &str, hint: &str, v: &mut bool) {
    ui.label(name).on_hover_text(hint);
    ui.checkbox(v, "");
    ui.end_row();
}

/// Editable list of inclusive m/z ranges.
fn ranges(ui: &mut egui::Ui, name: &str, hint: &str, list: &mut Vec<(f32, f32)>) {
    ui.add_space(4.0);
    ui.horizontal(|ui| {
        ui.label(name).on_hover_text(hint);
        if ui.small_button("+").on_hover_text("Add a range").clicked() {
            list.push((0.0, 0.0));
        }
    });
    let mut remove = None;
    for (i, (lo, hi)) in list.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            ui.add(
                egui::DragValue::new(lo)
                    .range(0.0..=100_000.0)
                    .max_decimals(2),
            );
            ui.label("–");
            ui.add(
                egui::DragValue::new(hi)
                    .range(0.0..=100_000.0)
                    .max_decimals(2),
            );
            if ui.small_button("×").clicked() {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
}
