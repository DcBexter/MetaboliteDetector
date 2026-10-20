//! Settings of the old version: a dialog with plain names (the old `AN_*`, `LS_*`, … keys on hover),
//! the old defaults, load and save of the old settings file, and the background re-analysis of open
//! runs when they change.

use std::io::Read;
use std::ops::RangeInclusive;
use std::path::Path;
use std::sync::mpsc;

use crate::layout::ToolWindow;
use eframe::egui;
use md_core::Chromatogram;
use md_core::deconv::{AUTO_WIDTH_S, Component, DeconvSettings, auto_width, deconvolute};
use md_core::peaks::{Detector, IonTraces, Peak, PeakSettings};
use md_core::settings::{GcmsSettings, MAX_LEN};

use crate::{Run, spawn};

/// What the components of a run were computed with.
pub type Key = (PeakSettings, DeconvSettings);

/// The old Appearance options (QSettings `GUI/axes_label` and `GUI/peak_label`, both on by default).
/// They only change the drawing, so they apply at once.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Appearance {
    /// Axis titles of the plots.
    pub axes_labels: bool,
    /// The m/z of the highest peaks in view of the spectrum.
    pub peak_labels: bool,
}

impl Default for Appearance {
    fn default() -> Self {
        Self {
            axes_labels: true,
            peak_labels: true,
        }
    }
}

/// How many of the highest peaks in view get a label (old `SpectrumPlot`: 5).
pub const LABELLED_PEAKS: usize = 5;

/// A library the user wants to choose from the Settings window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Pick {
    Library,
    Nist,
}

/// (LS_USE_SPEC, LS_USE_RI) and the name of each score choice of the old dialog.
const SCORES: [(bool, bool, &str); 3] = [
    (true, true, "Combined"),
    (true, false, "Spectrum"),
    (false, true, "RI"),
];

#[derive(Default)]
pub struct Settings {
    /// In use: new runs open with it and open runs are re-analysed to match it.
    pub applied: GcmsSettings,
    /// Being edited in the dialog until "Apply".
    draft: GcmsSettings,
    pub appearance: Appearance,
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
    let mut traces = IonTraces::new(data).map_err(|e| e.to_string())?;
    let peaks = traces.detect_all_peaks(&key.0);
    let mut components = deconvolute(&traces, &peaks, &data.scan_times_s, &key.1);
    components.sort_by_key(|k| k.apex_scan);
    let masses = peaks
        .iter()
        .map(|p| traces.masses[p.row as usize])
        .collect();
    Ok((components, masses, peaks))
}

impl Settings {
    /// The saved settings at start, checked by `persist`.
    pub fn restore(&mut self, s: GcmsSettings) {
        self.draft = s.clone();
        self.applied = s;
    }

    /// Opens the dialog with a message about the saved or imported settings.
    pub fn show_note(&mut self, note: Result<String, String>) {
        (self.note, self.open) = (Some(note), true);
    }

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

    /// The window; `libraries` are the paths of the open compound and NIST library, `times` the scan
    /// times of the shown run. Returns the library the user wants to choose.
    pub fn window(
        &mut self,
        ctx: &egui::Context,
        libraries: [Option<&Path>; 2],
        mode: &mut crate::guided::Mode,
        times: Option<&[f64]>,
    ) -> Option<Pick> {
        let mut open = self.open;
        let mut pick = None;
        egui::Window::new("Settings")
            .tool_window(ctx, 0, [560.0, 680.0])
            .open(&mut open)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical()
                    .max_height(ui.available_height() - 60.0)
                    .show(ui, |ui| {
                        pick = self.libraries(ui, libraries);
                        ui.strong(crate::guided::text::MODE);
                        crate::guided::switch(ui, mode);
                        ui.add_space(8.0);
                        // An invalid value must not hide under a closed "Advanced" (#303).
                        let invalid = self.draft.validate().is_err();
                        crate::guided::advanced_open(ui, "settings", invalid, |ui| {
                            self.fields(ui, times)
                        });
                    });
                ui.separator();
                self.buttons(ui, times);
            });
        self.open = open;
        pick
    }

    /// Libraries and appearance: they take effect at once, without "Apply".
    fn libraries(&mut self, ui: &mut egui::Ui, paths: [Option<&Path>; 2]) -> Option<Pick> {
        let mut pick = None;
        ui.strong("Libraries");
        egui::Grid::new("libraries").num_columns(3).show(ui, |ui| {
            for (name, path, which, hint) in [
                (
                    "Compound library",
                    paths[0],
                    Pick::Library,
                    "Components are identified against it. Opened again at the next start.",
                ),
                (
                    "NIST library",
                    paths[1],
                    Pick::Nist,
                    "Used by \"NIST search\" on the selected component. Opened again at the next start.",
                ),
            ] {
                ui.label(name).on_hover_text(hint);
                match path {
                    Some(p) => ui.label(p.display().to_string()),
                    None => ui.weak("none"),
                };
                if ui.button("Choose…").on_hover_text(hint).clicked() {
                    pick = Some(which);
                }
                ui.end_row();
            }
        });
        ui.add_space(8.0);
        ui.strong("Appearance");
        let a = &mut self.appearance;
        ui.checkbox(&mut a.axes_labels, "Axes labels").on_hover_text(
            "Axis titles on the spectrum, the Detected Peaks plot and the y axis of the chromatograms",
        );
        ui.checkbox(&mut a.peak_labels, "Label highest peaks")
            .on_hover_text(format!(
                "The {LABELLED_PEAKS} highest peaks in view of the spectrum show their m/z"
            ));
        ui.add_space(8.0);
        pick
    }

    fn buttons(&mut self, ui: &mut egui::Ui, times: Option<&[f64]>) {
        let valid = self.draft.validate();
        let fix = "Fix the value shown in red first (it names the old key)";
        ui.horizontal(|ui| {
            let changed = self.draft != self.applied;
            if ui
                .add_enabled(valid.is_ok() && changed, egui::Button::new("Apply"))
                .on_hover_text("Re-runs peak detection, deconvolution and the library search")
                .on_disabled_hover_text(if valid.is_ok() {
                    "Nothing changed"
                } else {
                    fix
                })
                .clicked()
            {
                self.applied = self.draft.clone();
                self.note = None;
            }
            if ui
                .button("Old defaults")
                .on_hover_text("Fill in the old version's defaults, with the deconvolution width on Auto; Apply uses them")
                .clicked()
            {
                self.draft = GcmsSettings::default();
            }
            if ui
                .button("Load…")
                .on_hover_text("Read an old GCMSSettings file and apply it")
                .clicked()
                && let Some(path) = dialog().pick_file()
            {
                self.note = Some(self.load(&path));
            }
            if ui
                .add_enabled(valid.is_ok(), egui::Button::new("Save…"))
                .on_hover_text("Write these values as an old GCMSSettings file")
                .on_disabled_hover_text(fix)
                .clicked()
                && let Some(path) = dialog().set_file_name("GCMSSettings").save_file()
            {
                // The old file has no Auto: it gets the width of the shown run.
                let mut file = self.draft.clone();
                (file.an_deconvolution_width).get_or_insert_with(|| auto_width(times.unwrap_or_default()));
                self.note = Some(
                    std::fs::write(&path, file.to_bytes())
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

    fn fields(&mut self, ui: &mut egui::Ui, times: Option<&[f64]>) {
        let d = &mut self.draft;
        ui.strong("Peak detection and deconvolution");
        egui::Grid::new("an").num_columns(2).show(ui, |ui| {
            let hint = "Default: the old version's slope detector. \
                Matched Filtering: the old Calibration Setter's matched peak detector.";
            let label = ui.label("Peak detection method").on_hover_text(hint);
            egui::ComboBox::from_id_salt("detector")
                .selected_text(match d.peak_detector {
                    Detector::Default => "Default",
                    Detector::Matched { .. } => "Matched Filtering",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut d.peak_detector, Detector::Default, "Default");
                    // Choosing it again keeps the width and threshold.
                    let matched = matches!(d.peak_detector, Detector::Matched { .. });
                    if ui.selectable_label(matched, "Matched Filtering").clicked() && !matched {
                        d.peak_detector = Detector::MATCHED;
                    }
                })
                .response
                .labelled_by(label.id)
                .on_hover_text(hint);
            ui.end_row();
            if let Detector::Matched { width, threshold } = &mut d.peak_detector {
                row(
                    ui,
                    "Peak width",
                    "",
                    "Scans; the filter density is 4 × width. The slopes and Minimum S/N below are not used.",
                    egui::DragValue::new(width).range(1..=1000),
                );
                row(
                    ui,
                    "Peak threshold",
                    "",
                    "A peak is a region where the filtered trace lies above this.",
                    egui::DragValue::new(threshold).range(-1_000_000_000..=1_000_000_000),
                );
            }
            row(
                ui,
                "Peak start slope",
                "AN_PEAK_THRESHOLD_BEGIN",
                "Smoothed slope that starts a peak.",
                num(&mut d.an_peak_threshold_begin, ANY),
            );
            row(
                ui,
                "Peak end slope",
                "AN_PEAK_THRESHOLD_END",
                "The falling flank ends once the slope stays above this.",
                num(&mut d.an_peak_threshold_end, ANY),
            );
            row(
                ui,
                "Minimum S/N",
                "AN_MIN_PEAK_HEIGHT",
                "Minimum height above the baseline, in noise units.",
                num(&mut d.an_min_peak_height, 0.0..=1e9),
            );
            row(
                ui,
                "Apex steps per scan",
                "AN_NUMBER_BINS",
                "Sub-scan steps for locating the precise apex.",
                egui::DragValue::new(&mut d.an_number_bins).range(1..=100),
            );
            width_row(ui, &mut d.an_deconvolution_width, times);
            row(
                ui,
                "Minimum ions per component",
                "AN_MIN_PEAK_NUMBER",
                "Components with fewer ion peaks are dropped.",
                egui::DragValue::new(&mut d.an_min_peak_number).range(1..=10_000),
            );
            row(
                ui,
                "Minimum share of base peak",
                "AN_MIN_HEIGHT_BP",
                "Ions below this fraction of the base peak leave the spectrum (0 = keep all).",
                num(&mut d.an_min_height_bp, 0.0..=1.0),
            );
            row(
                ui,
                "Integrate fitted model",
                "AN_MODELED_VALUES_FOR_INTEGRATION",
                "Integrate the fitted model instead of the raw trace.",
                egui::Checkbox::without_text(&mut d.an_modeled_values_for_integration),
            );
        });
        ranges(
            ui,
            "Ions left out of the model peak",
            &explain(
                "Ion ranges (m/z from–to) not used for the model peak.",
                "QUANT_EXCL_IONS",
            ),
            &mut d.quant_excl_ions,
        );

        ui.add_space(8.0);
        ui.strong("Library search");
        egui::Grid::new("ls").num_columns(2).show(ui, |ui| {
            row(
                ui,
                "RI tolerance (sigma)",
                "LS_RI_DIFF",
                "Sigma of the RI score; differences above 3 sigma score 0.",
                num(&mut d.ls_ri_diff, 1e-6..=1e9),
            );
            row(
                ui,
                "Pure spectrum weight",
                "LS_PURE_FACTOR",
                "Weight of the pure spectrum score.",
                num(&mut d.ls_pure_factor, 0.0..=1e9),
            );
            row(
                ui,
                "Impure spectrum weight",
                "LS_IMPURE_FACTOR",
                "Weight of the impure spectrum score.",
                num(&mut d.ls_impure_factor, 0.0..=1e9),
            );
            row(
                ui,
                "Mass exponent",
                "LS_SCALE_FACTOR_ION",
                "Peaks are weighted as mass^a · intensity^b: this is a.",
                num(&mut d.ls_scale_factor_ion, ANY),
            );
            row(
                ui,
                "Intensity exponent",
                "LS_SCALE_FACTOR_VALUE",
                "Peaks are weighted as mass^a · intensity^b: this is b.",
                num(&mut d.ls_scale_factor_value, ANY),
            );
            row(
                ui,
                "Unmatched peak factor",
                "LS_SCALE_FACTOR_NO_MATCH",
                "Factor for peaks without a partner in the other spectrum.",
                num(&mut d.ls_scale_factor_no_match, ANY),
            );
            let hint = explain(
                "What the overall score is made of. Combined: spectrum score^(2/3) × RI score^(1/3).",
                "LS_USE_SPEC and LS_USE_RI",
            );
            let label = ui.label("Score").on_hover_text(&hint);
            let now = (d.ls_use_spec, d.ls_use_ri);
            // An old file can switch both off: then every hit scores 1.
            let name = (SCORES.iter().find(|m| (m.0, m.1) == now))
                .map_or("None (every hit scores 1)", |m| m.2);
            egui::ComboBox::from_id_salt("score")
                .selected_text(name)
                .show_ui(ui, |ui| {
                    for (spec, ri, name) in SCORES {
                        if ui.selectable_label(now == (spec, ri), name).clicked() {
                            (d.ls_use_spec, d.ls_use_ri) = (spec, ri);
                        }
                    }
                })
                .response
                .labelled_by(label.id)
                .on_hover_text(hint);
            ui.end_row();
            row(
                ui,
                "Cutoff score",
                "LS_CUTOFF",
                "A hit below this overall score counts as no match: in the component list, Find substance and the batch.",
                num(&mut d.ls_cutoff, 0.0..=1.0),
            );
            row(
                ui,
                "Number of fragments",
                "SQL_MAX_NUMBER_FRAGMENTS",
                "NIST search: how many of the component's strongest fragments are compared.",
                num(&mut d.sql_max_number_fragments, 1..=100),
            );
            row(
                ui,
                "Minimum number of identical frag.",
                "SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS",
                "NIST search: a library entry is scored only if it shares more than this many of them.",
                num(&mut d.sql_minimal_number_equal_fragments, 0..=100),
            );
        });
        ranges(
            ui,
            "Ignored masses",
            &explain(
                "Mass ranges (m/z from–to) ignored in both spectra.",
                "LS_MASS_FILTER",
            ),
            &mut d.ls_mass_filter,
        );

        ui.add_space(8.0);
        ui.strong("Centroiding (export)");
        egui::Grid::new("centroid").num_columns(2).show(ui, |ui| {
            let used = "Finds the peaks of a profile spectrum when it is exported as centroids.";
            for (name, key, v) in [
                (
                    "Centroid: peak start",
                    "SF_CENTROID_THRESHOLD_BEGIN",
                    &mut d.sf_centroid_threshold_begin,
                ),
                (
                    "Centroid: peak end",
                    "SF_CENTROID_THRESHOLD_END",
                    &mut d.sf_centroid_threshold_end,
                ),
                (
                    "Centroid: baseline deviation",
                    "SF_CENTROID_MAX_BASE_DEV",
                    &mut d.sf_centroid_max_base_dev,
                ),
            ] {
                row(ui, name, key, used, num(v, ANY));
            }
        });

        ui.add_space(8.0);
        egui::CollapsingHeader::new("Stored in the file, not used").show(ui, |ui| {
            let unused = "Kept in the settings file; this version does not use it.";
            // RE of 3.3 (docs/reverse-engineering/algorithms.md §4b.4): only saved and shown there too.
            let dead = "Kept in the settings file; 3.3 does not use it either.";
            egui::Grid::new("stored").num_columns(2).show(ui, |ui| {
                for (name, key, v) in [
                    (
                        "Spectrum filter: peak start",
                        "SF_PEAK_THRESHOLD_BEGIN",
                        &mut d.sf_peak_threshold_begin,
                    ),
                    (
                        "Spectrum filter: peak end",
                        "SF_PEAK_THRESHOLD_END",
                        &mut d.sf_peak_threshold_end,
                    ),
                    (
                        "Spectrum filter: baseline deviation",
                        "SF_MAX_BASE_DEV",
                        &mut d.sf_max_base_dev,
                    ),
                ] {
                    row(ui, name, key, dead, num(v, ANY));
                }
                row(
                    ui,
                    "Spectrum filter: FWHM",
                    "SF_FWHM",
                    "3.3 centroids with the points at or above this share of a peak's maximum; this version always uses 0.5.",
                    num(&mut d.sf_fwhm, ANY),
                );
                let stored: [(&str, &str, &mut f64); 3] = [
                    (
                        "Peak baseline deviation",
                        "AN_MAX_BASE_DEV",
                        &mut d.an_max_base_dev,
                    ),
                    (
                        "Quant ion minimum distance",
                        "QUANT_MIN_DIST_IONS",
                        &mut d.quant_min_dist_ions,
                    ),
                    (
                        "Quant quality cutoff",
                        "QUANT_QUALITY_CUTOFF",
                        &mut d.quant_quality_cutoff,
                    ),
                ];
                for (name, key, v) in stored {
                    row(ui, name, key, unused, num(v, ANY));
                }
                row(
                    ui,
                    "Search resolution",
                    "LS_RESOLUTION",
                    unused,
                    num(&mut d.ls_resolution, -1e9..=1e9),
                );
                row(
                    ui,
                    "Use scaled library",
                    "LS_USE_SCALED_LIB",
                    dead,
                    egui::Checkbox::without_text(&mut d.ls_use_scaled_lib),
                );
                row(
                    ui,
                    "Quant ions",
                    "QUANT_ION_NUMBER",
                    unused,
                    egui::DragValue::new(&mut d.quant_ion_number).range(0..=10_000),
                );
            });
        });
    }
}

/// The old file has no extension convention.
fn dialog() -> rfd::FileDialog {
    rfd::FileDialog::new()
}

const ANY: RangeInclusive<f64> = -1e9..=1e9;

/// What a setting does, then its key in the old 3.3 settings file.
fn explain(hint: &str, key: &str) -> String {
    match key {
        "" => hint.to_owned(),
        _ => format!("{hint}\nOld key: {key}"),
    }
}

/// One grid row: the plain name, then the field it names; both explain on hover.
pub(crate) fn row(ui: &mut egui::Ui, name: &str, key: &str, hint: &str, field: impl egui::Widget) {
    let hint = explain(hint, key);
    let label = ui.label(name).on_hover_text(&hint);
    ui.add(field).labelled_by(label.id).on_hover_text(hint);
    ui.end_row();
}

/// The deconvolution width row: Auto with the scans it gives in the shown run, or a set width.
fn width_row(ui: &mut egui::Ui, width: &mut Option<f64>, times: Option<&[f64]>) {
    let hint = explain(
        &format!(
            "Width in scans of the matched filter over the apex histogram. Auto: {AUTO_WIDTH_S} s at \
             each run's scan rate, at least 1 scan, at most 100 and a tenth of the run (3.3 always used 1). A width from a loaded or \
             imported 3.3 file is a set one; Save writes Auto as the width of the shown run (1 without one)."
        ),
        "AN_DECONVOLUTION_WIDTH",
    );
    let label = ui.label("Deconvolution width").on_hover_text(&hint);
    ui.horizontal(|ui| {
        let mut auto = width.is_none();
        let check = ui.checkbox(&mut auto, "Auto");
        if check.labelled_by(label.id).on_hover_text(&hint).changed() {
            // A set width starts from the automatic one.
            *width = (!auto).then(|| auto_width(times.unwrap_or_default()));
        }
        match width {
            Some(w) => ui.add(num(w, 0.1..=100.0)).labelled_by(label.id),
            None => match times {
                Some(t) => ui.label(format!("{} scans in the shown run", auto_width(t))),
                None => ui.weak(format!("{AUTO_WIDTH_S} s of each run")),
            },
        }
        .on_hover_text(&hint);
    });
    ui.end_row();
}

fn num<N: egui::emath::Numeric>(v: &mut N, range: RangeInclusive<N>) -> egui::DragValue<'_> {
    egui::DragValue::new(v)
        .range(range)
        .speed(0.1)
        .max_decimals(6)
}

/// Editable list of inclusive m/z ranges.
pub(crate) fn ranges(ui: &mut egui::Ui, name: &str, hint: &str, list: &mut Vec<(f32, f32)>) {
    ui.add_space(4.0);
    let list_name = ui.horizontal(|ui| {
        let list_name = ui.label(name).on_hover_text(hint);
        if ui.small_button("+").on_hover_text("Add a range").clicked() {
            list.push((0.0, 0.0));
        }
        list_name.id
    });
    let list_name = list_name.inner;
    let mut remove = None;
    for (i, (lo, hi)) in list.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            let from = ui.label("from");
            ui.add(
                egui::DragValue::new(lo)
                    .range(0.0..=100_000.0)
                    .max_decimals(2),
            )
            .labelled_by(list_name)
            .labelled_by(from.id)
            .on_hover_text(hint);
            let to = ui.label("to");
            ui.add(
                egui::DragValue::new(hi)
                    .range(0.0..=100_000.0)
                    .max_decimals(2),
            )
            .labelled_by(list_name)
            .labelled_by(to.id)
            .on_hover_text(hint);
            if ui
                .small_button("×")
                .on_hover_text("Remove this range")
                .clicked()
            {
                remove = Some(i);
            }
        });
    }
    if let Some(i) = remove {
        list.remove(i);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use egui_kittest::Harness;
    use egui_kittest::kittest::Queryable;

    /// #303: in Guided mode "Advanced" is closed, but an invalid draft opens it.
    #[test]
    fn invalid_value_opens_closed_advanced() {
        let s = Settings {
            open: true,
            ..Default::default()
        };
        let mut h = Harness::new_ui_state(
            |ui, (s, mode): &mut (Settings, crate::guided::Mode)| {
                let ctx = ui.ctx().clone();
                let ctx = &ctx;
                crate::guided::publish(ctx, *mode);
                s.window(ctx, [None, None], mode, None);
            },
            (s, crate::guided::Mode::Guided),
        );
        h.run();
        assert!(
            h.query_by_label("Deconvolution width").is_none(),
            "Advanced starts closed"
        );
        h.state_mut().0.draft.sf_fwhm = f64::NAN;
        h.run();
        h.get_by_label("Deconvolution width");
    }

    /// #279: Auto shows the width of the shown run; unticking it sets that width.
    #[test]
    fn auto_width_shows_the_scans_and_unticks_to_them() {
        let times: Vec<f64> = (0..400).map(|i| i as f64 / 20.0).collect();
        let mut h = Harness::new_ui_state(
            |ui, s: &mut Settings| s.fields(ui, Some(&times)),
            Settings::default(),
        );
        h.run();
        h.get_by_label("10 scans in the shown run");
        // The row label names the Auto checkbox.
        h.get_by_role_and_label(egui::accesskit::Role::CheckBox, "Deconvolution width")
            .click();
        h.run();
        assert_eq!(h.state().draft.an_deconvolution_width, Some(10.0));
    }
}
