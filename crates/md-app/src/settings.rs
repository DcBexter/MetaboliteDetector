//! Settings of the old version: a dialog with plain names (the old `AN_*`, `LS_*`, … keys on hover),
//! the old defaults, load and save of the old settings file, and the background re-analysis of open
//! runs when they change.

use std::io::Read;
use std::ops::RangeInclusive;
use std::path::Path;

use crate::errors;
use crate::lang::{Msg, msg, t, tf};
use crate::layout::ToolWindow;
use eframe::egui;
use md_core::Chromatogram;
use md_core::analysis::{PeakPass, analyze_run};
use md_core::deconv::{AUTO_WIDTH_S, Component, DeconvSettings, auto_width};
use md_core::job::Cancel;
use md_core::peaks::{Detector, Peak, PeakSettings};
use md_core::settings::{ANY, GcmsSettings, MAX_LEN, keys, limits};
use md_core::view::LABELLED_PEAKS;

use crate::{Job, Run};

/// What the analysis of a run depends on: peak detection, deconvolution and the centroid `fwhm`.
#[derive(Debug, Clone, PartialEq)]
pub struct Key {
    pub peaks: PeakSettings,
    pub deconv: DeconvSettings,
    pub fwhm: f64,
}

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

/// A library the user wants to choose from the Settings window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Pick {
    Library,
    Nist,
}

/// (LS_USE_SPEC, LS_USE_RI) and the name of each score choice of the old dialog.
fn scores() -> [(bool, bool, &'static str); 3] {
    [
        (true, true, t("Combined", "Kombiniert")),
        (true, false, t("Spectrum", "Spektrum")),
        (false, true, "RI"),
    ]
}

#[derive(Default)]
pub struct Settings {
    /// In use: new runs open with it and open runs are re-analysed to match it.
    pub applied: GcmsSettings,
    /// Being edited in the dialog until "Apply".
    draft: GcmsSettings,
    pub appearance: Appearance,
    pub open: bool,
    note: Option<Result<Msg, Msg>>,
}

type Analyzed = Result<(Vec<Component>, Vec<u32>, Vec<Peak>), String>;

/// Re-analysis state of one run.
#[derive(Default)]
pub struct Analysis {
    /// Settings behind `Run::components`.
    pub key: Option<Key>,
    job: Option<(Key, Job<Analyzed>)>,
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
pub fn analyze(data: &Chromatogram, key: &Key, cancel: &Cancel) -> Analyzed {
    let a = analyze_run(
        data,
        &key.peaks,
        key.fwhm,
        &key.deconv,
        PeakPass::Reordered,
        cancel,
    )
    .map_err(|e| e.to_string())?;
    Ok((a.components, a.masses, a.peaks))
}

impl Settings {
    /// The saved settings at start, checked by `persist`.
    pub fn restore(&mut self, s: GcmsSettings) {
        self.draft = s.clone();
        self.applied = s;
    }

    /// Opens the dialog with a message about the saved or imported settings.
    pub fn show_note(&mut self, note: Result<Msg, Msg>) {
        (self.note, self.open) = (Some(note), true);
    }

    pub fn key(&self) -> Key {
        let a = self.applied.analysis();
        Key {
            peaks: a.peaks,
            deconv: a.deconv,
            fwhm: a.fwhm,
        }
    }

    /// Collects finished analyses and starts one for every run whose components are stale.
    /// Returns the error of an analysis that failed (the run keeps its old components).
    pub fn update_runs(&self, ctx: &egui::Context, runs: &mut [Run]) -> Option<Msg> {
        let wanted = self.key();
        let correct = self.applied.baseline_correction;
        let mut error = None;
        for run in runs {
            // Runs without stored means (CDF, other caches) stay as they are, components included.
            if !run.data.baseline.means.is_empty() && run.data.baseline.corrected != correct {
                // Copies the run only while a worker still reads the old data.
                std::sync::Arc::make_mut(&mut run.data).set_baseline_correction(correct);
                // Components, SIC and the shown spectrum come from the corrected values; the TIC does not.
                (run.analysis.key, run.analysis.job) = (None, None);
                (run.sic_for, run.sic_job) = (None, None);
            }
            let a = &mut run.analysis;
            if let Some((key, job)) = &a.job {
                let what = || format!("{}: {}", run.name, t("analysis", "Analyse"));
                match job.poll(what) {
                    Some(Ok(result)) => {
                        match result {
                            Ok((components, masses, peaks)) => {
                                (run.components, run.peak_masses, run.peaks) =
                                    (components, masses, peaks);
                                // Selection and library hits belong to the old components.
                                (run.component, run.ident, run.ident_job) = (None, None, None);
                                run.rev += 1;
                            }
                            Err(e) => error = Some(format!("{}: {e}", run.name).into()),
                        }
                        a.key = Some(key.clone());
                        a.job = None;
                    }
                    None => {}
                    Some(Err(e)) => {
                        // Panicked: keep the old components and do not retry these settings.
                        error = Some(e);
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
                a.job = Some((
                    wanted.clone(),
                    Job::spawn(ctx, move |_, cancel| analyze(&data, &key, cancel)),
                ));
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
        egui::Window::new(t("Settings", "Einstellungen"))
            .id(egui::Id::unique("Settings"))
            .tool_window(ctx, 0, [560.0, 680.0])
            .open(&mut open)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical()
                    .max_height(ui.available_height() - 60.0)
                    .show(ui, |ui| {
                        pick = self.libraries(ui, libraries);
                        ui.strong(crate::guided::text::mode());
                        crate::guided::switch(ui, mode);
                        crate::lang::switch(ui);
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
        ui.strong(t("Libraries", "Bibliotheken"));
        egui::Grid::new("libraries").num_columns(3).show(ui, |ui| {
            for (name, path, which, hint) in [
                (
                    t("Compound library", "Substanzbibliothek"),
                    paths[0],
                    Pick::Library,
                    t(
                        "Components are identified against it. Opened again at the next start.",
                        "Die Komponenten werden mit ihr identifiziert. Beim nächsten Start wieder geöffnet.",
                    ),
                ),
                (
                    t("NIST library", "NIST-Bibliothek"),
                    paths[1],
                    Pick::Nist,
                    t(
                        "Used by \"NIST search\" on the selected component. Opened again at the next start.",
                        "Für die „NIST-Suche“ zur gewählten Komponente. Beim nächsten Start wieder geöffnet.",
                    ),
                ),
            ] {
                ui.label(name).on_hover_text(hint);
                match path {
                    Some(p) => ui.label(p.display().to_string()),
                    None => ui.weak(t("none", "keine")),
                };
                if ui.button(t("Choose…", "Wählen…")).on_hover_text(hint).clicked() {
                    pick = Some(which);
                }
                ui.end_row();
            }
        });
        ui.add_space(8.0);
        ui.strong(t("Appearance", "Darstellung"));
        let a = &mut self.appearance;
        ui.checkbox(&mut a.axes_labels, t("Axes labels", "Achsentitel"))
            .on_hover_text(t(
                "Axis titles on the spectrum, the Detected Peaks plot and the y axis of the chromatograms",
                "Achsentitel am Spektrum, am Plot der erkannten Peaks und an der y-Achse der Chromatogramme",
            ));
        ui.checkbox(
            &mut a.peak_labels,
            t("Label highest peaks", "Höchste Peaks beschriften"),
        )
        .on_hover_text(tf!(
            "The {LABELLED_PEAKS} highest peaks in view of the spectrum show their m/z",
            "Die {LABELLED_PEAKS} höchsten sichtbaren Peaks des Spektrums zeigen ihr m/z"
        ));
        ui.add_space(8.0);
        pick
    }

    fn buttons(&mut self, ui: &mut egui::Ui, times: Option<&[f64]>) {
        let valid = self.draft.validate();
        let fix = t(
            "Fix the value shown in red first (it names the old key)",
            "Erst den rot angezeigten Wert korrigieren (er nennt den alten Schlüssel)",
        );
        ui.horizontal(|ui| {
            let changed = self.draft != self.applied;
            if ui
                .add_enabled(valid.is_ok() && changed, egui::Button::new(t("Apply", "Übernehmen")))
                .on_hover_text(t(
                    "Re-runs peak detection, deconvolution and the library search",
                    "Führt Peakerkennung, Dekonvolution und Bibliothekssuche neu aus",
                ))
                .on_disabled_hover_text(if valid.is_ok() {
                    t("Nothing changed", "Nichts geändert")
                } else {
                    fix
                })
                .clicked()
            {
                self.applied = self.draft.clone();
                self.note = None;
            }
            if ui
                .button(t("Old defaults", "Alte Standardwerte"))
                .on_hover_text(t(
                    "Fill in the old version's defaults, with the deconvolution width on Auto; Apply uses them",
                    "Die Standardwerte der alten Version eintragen, Dekonvolutionsbreite auf Auto; Übernehmen verwendet sie",
                ))
                .clicked()
            {
                self.draft = GcmsSettings::default();
            }
            if ui
                .button(t("Load…", "Laden…"))
                .on_hover_text(t(
                    "Read an old GCMSSettings file and apply it",
                    "Eine alte GCMSSettings-Datei lesen und übernehmen",
                ))
                .clicked()
                && let Some(path) = dialog().pick_file()
            {
                self.note = Some(self.load(&path));
            }
            if ui
                .add_enabled(valid.is_ok(), egui::Button::new(t("Save…", "Speichern…")))
                .on_hover_text(t(
                    "Write these values as an old GCMSSettings file",
                    "Diese Werte als alte GCMSSettings-Datei schreiben",
                ))
                .on_disabled_hover_text(fix)
                .clicked()
                && let Some(path) = dialog().set_file_name("GCMSSettings").save_file()
            {
                // The old file has no Auto: it gets the width of the shown run.
                let mut file = self.draft.clone();
                (file.an_deconvolution_width).get_or_insert_with(|| auto_width(times.unwrap_or_default()));
                self.note = Some(
                    std::fs::write(&path, file.to_bytes())
                        .map(|()| msg!("Saved {}", "Gespeichert: {}", path.display()))
                        .map_err(|e| format!("{}: {e}", path.display()).into()),
                );
            }
        });
        let error = ui.visuals().error_fg_color;
        match (&valid, &self.note) {
            (Err(e), _) => {
                ui.colored_label(error, errors::text(e));
            }
            (_, Some(Err(e))) => {
                ui.colored_label(error, &**e);
            }
            (_, Some(Ok(m))) => {
                ui.label(&**m);
            }
            _ => {}
        }
    }

    /// Reads an old settings file into the dialog and applies it. The file is untrusted:
    /// `from_bytes` checks size, checksum and values, and reading stops at `MAX_LEN`.
    fn load(&mut self, path: &std::path::Path) -> Result<Msg, Msg> {
        let err = |e: &dyn std::fmt::Display| Msg::from(format!("{}: {e}", path.display()));
        let mut bytes = Vec::new();
        std::fs::File::open(path)
            .and_then(|f| f.take(MAX_LEN + 1).read_to_end(&mut bytes))
            .map_err(|e| err(&e))?;
        if bytes.len() as u64 > MAX_LEN {
            return Err(msg!(
                "{}: file too large for a settings file",
                "{}: Datei zu groß für eine Einstellungsdatei",
                path.display()
            ));
        }
        self.draft = self.applied.from_bytes(&bytes).map_err(|e| err(&e))?;
        self.applied = self.draft.clone();
        Ok(msg!("Loaded {}", "Geladen: {}", path.display()))
    }

    fn fields(&mut self, ui: &mut egui::Ui, times: Option<&[f64]>) {
        let d = &mut self.draft;
        ui.strong(t(
            "Peak detection and deconvolution",
            "Peakerkennung und Dekonvolution",
        ));
        egui::Grid::new("an").num_columns(2).show(ui, |ui| {
            let hint = t(
                "Default: the old version's slope detector. \
                 Matched Filtering: the old Calibration Setter's matched peak detector.",
                "Standard: der Steigungsdetektor der alten Version. \
                 Matched Filtering: der Matched-Peak-Detektor des alten Calibration Setter.",
            );
            let label = ui.label(name_of(keys::PEAK_DETECTOR)).on_hover_text(hint);
            egui::ComboBox::from_id_salt("detector")
                .selected_text(match d.peak_detector {
                    Detector::Default => t("Default", "Standard"),
                    Detector::Matched { .. } => "Matched Filtering",
                })
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut d.peak_detector,
                        Detector::Default,
                        t("Default", "Standard"),
                    );
                    // Choosing it again keeps the width and threshold.
                    let matched = matches!(d.peak_detector, Detector::Matched { .. });
                    if ui.selectable_label(matched, t("Matched Filtering", "Matched Filtering")).clicked() && !matched {
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
                    t("Peak width", "Peakbreite"),
                    "",
                    t(
                        "Scans; the filter density is 4 × width. The slopes and Minimum S/N below are not used.",
                        "Scans; die Filterdichte ist 4 × Breite. Die Steigungen und das Minimum-S/N unten werden nicht verwendet.",
                    ),
                    egui::DragValue::new(width).range(limit::<u32>(keys::PEAK_DETECTOR)),
                );
                row(
                    ui,
                    t("Peak threshold", "Peakschwelle"),
                    "",
                    t(
                        "A peak is a region where the filtered trace lies above this.",
                        "Ein Peak ist ein Bereich, in dem die gefilterte Spur darüber liegt.",
                    ),
                    egui::DragValue::new(threshold).range(bounds::<i32>(ANY)),
                );
            }
            row(
                ui,
                name_of(keys::AN_PEAK_THRESHOLD_BEGIN),
                keys::AN_PEAK_THRESHOLD_BEGIN,
                t(
                    "Smoothed slope that starts a peak.",
                    "Geglättete Steigung, mit der ein Peak beginnt.",
                ),
                num(&mut d.an_peak_threshold_begin, limit(keys::AN_PEAK_THRESHOLD_BEGIN)),
            );
            row(
                ui,
                name_of(keys::AN_PEAK_THRESHOLD_END),
                keys::AN_PEAK_THRESHOLD_END,
                t(
                    "The falling flank ends once the slope stays above this.",
                    "Die fallende Flanke endet, sobald die Steigung darüber bleibt.",
                ),
                num(&mut d.an_peak_threshold_end, limit(keys::AN_PEAK_THRESHOLD_END)),
            );
            row(
                ui,
                name_of(keys::AN_MIN_PEAK_HEIGHT),
                keys::AN_MIN_PEAK_HEIGHT,
                t(
                    "Minimum height above the baseline, in noise units.",
                    "Mindesthöhe über der Basislinie, in Rauscheinheiten.",
                ),
                num(&mut d.an_min_peak_height, limit(keys::AN_MIN_PEAK_HEIGHT)),
            );
            row(
                ui,
                name_of(keys::AN_NUMBER_BINS),
                keys::AN_NUMBER_BINS,
                t(
                    "Sub-scan steps for locating the precise apex.",
                    "Teilschritte innerhalb eines Scans, um den Apex genau zu finden.",
                ),
                egui::DragValue::new(&mut d.an_number_bins).range(limit::<i32>(keys::AN_NUMBER_BINS)),
            );
            width_row(ui, &mut d.an_deconvolution_width, times);
            row(
                ui,
                name_of(keys::AN_MIN_PEAK_NUMBER),
                keys::AN_MIN_PEAK_NUMBER,
                t(
                    "Components with fewer ion peaks are dropped.",
                    "Komponenten mit weniger Ionen-Peaks werden verworfen.",
                ),
                egui::DragValue::new(&mut d.an_min_peak_number).range(limit::<i32>(keys::AN_MIN_PEAK_NUMBER)),
            );
            row(
                ui,
                name_of(keys::AN_MIN_HEIGHT_BP),
                keys::AN_MIN_HEIGHT_BP,
                t(
                    "Ions below this fraction of the base peak leave the spectrum (0 = keep all).",
                    "Ionen unter diesem Anteil am Basispeak fallen aus dem Spektrum (0 = alle behalten).",
                ),
                num(&mut d.an_min_height_bp, limit(keys::AN_MIN_HEIGHT_BP)),
            );
            row(
                ui,
                name_of(keys::AN_MODELED_VALUES_FOR_INTEGRATION),
                keys::AN_MODELED_VALUES_FOR_INTEGRATION,
                t(
                    "Integrate the fitted model instead of the raw trace.",
                    "Das angepasste Modell statt der Rohspur integrieren.",
                ),
                egui::Checkbox::without_text(&mut d.an_modeled_values_for_integration),
            );
            row(
                ui,
                name_of(keys::BASELINE_CORRECTION),
                keys::BASELINE_CORRECTION,
                t(
                    "Subtract each scan's stored baseline from its intensities. Only old type 1 caches (FastFlight .bin/.idx) store one; other runs stay unchanged.",
                    "Die gespeicherte Basislinie jedes Scans von seinen Intensitäten abziehen. Nur alte Typ-1-Caches (FastFlight .bin/.idx) speichern eine; andere Läufe bleiben unverändert.",
                ),
                egui::Checkbox::without_text(&mut d.baseline_correction),
            );
        });
        ranges(
            ui,
            keys::QUANT_EXCL_IONS,
            name_of(keys::QUANT_EXCL_IONS),
            &explain(
                t(
                    "Ion ranges (m/z from–to) not used for the model peak.",
                    "Ionenbereiche (m/z von–bis), die nicht für den Modellpeak verwendet werden.",
                ),
                keys::QUANT_EXCL_IONS,
            ),
            &mut d.quant_excl_ions,
        );

        ui.add_space(8.0);
        ui.strong(t("Library search", "Bibliothekssuche"));
        egui::Grid::new("ls").num_columns(2).show(ui, |ui| {
            row(
                ui,
                name_of(keys::LS_RI_DIFF),
                keys::LS_RI_DIFF,
                t(
                    "Sigma of the RI score; differences above 3 sigma score 0.",
                    "Sigma des RI-Scores; Abweichungen über 3 Sigma ergeben 0.",
                ),
                num(&mut d.ls_ri_diff, limit(keys::LS_RI_DIFF)),
            );
            row(
                ui,
                name_of(keys::LS_PURE_FACTOR),
                keys::LS_PURE_FACTOR,
                t(
                    "Weight of the pure spectrum score.",
                    "Gewicht des Scores des reinen Spektrums.",
                ),
                num(&mut d.ls_pure_factor, limit(keys::LS_PURE_FACTOR)),
            );
            row(
                ui,
                name_of(keys::LS_IMPURE_FACTOR),
                keys::LS_IMPURE_FACTOR,
                t(
                    "Weight of the impure spectrum score.",
                    "Gewicht des Scores des unreinen Spektrums.",
                ),
                num(&mut d.ls_impure_factor, limit(keys::LS_IMPURE_FACTOR)),
            );
            row(
                ui,
                name_of(keys::LS_SCALE_FACTOR_ION),
                keys::LS_SCALE_FACTOR_ION,
                t(
                    "Peaks are weighted as mass^a · intensity^b: this is a.",
                    "Peaks werden als Masse^a · Intensität^b gewichtet: dies ist a.",
                ),
                num(&mut d.ls_scale_factor_ion, limit(keys::LS_SCALE_FACTOR_ION)),
            );
            row(
                ui,
                name_of(keys::LS_SCALE_FACTOR_VALUE),
                keys::LS_SCALE_FACTOR_VALUE,
                t(
                    "Peaks are weighted as mass^a · intensity^b: this is b.",
                    "Peaks werden als Masse^a · Intensität^b gewichtet: dies ist b.",
                ),
                num(&mut d.ls_scale_factor_value, limit(keys::LS_SCALE_FACTOR_VALUE)),
            );
            row(
                ui,
                name_of(keys::LS_SCALE_FACTOR_NO_MATCH),
                keys::LS_SCALE_FACTOR_NO_MATCH,
                t(
                    "Factor for peaks without a partner in the other spectrum.",
                    "Faktor für Peaks ohne Partner im anderen Spektrum.",
                ),
                num(&mut d.ls_scale_factor_no_match, limit(keys::LS_SCALE_FACTOR_NO_MATCH)),
            );
            let hint = explain(
                t(
                    "What the overall score is made of. Combined: spectrum score^(2/3) × RI score^(1/3).",
                    "Woraus sich der Gesamt-Score zusammensetzt. Kombiniert: Spektrum-Score^(2/3) × RI-Score^(1/3).",
                ),
                &format!(
                    "{} {} {}",
                    keys::LS_USE_SPEC,
                    t("and", "und"),
                    keys::LS_USE_RI
                ),
            );
            let label = ui.label(name_of(keys::LS_USE_SPEC)).on_hover_text(&hint);
            let now = (d.ls_use_spec, d.ls_use_ri);
            // An old file can switch both off: then every hit scores 1.
            let name = (scores().iter().find(|m| (m.0, m.1) == now)).map_or(
                t("None (every hit scores 1)", "Keiner (jeder Treffer ergibt 1)"),
                |m| m.2,
            );
            egui::ComboBox::from_id_salt("score")
                .selected_text(name)
                .show_ui(ui, |ui| {
                    for (spec, ri, name) in scores() {
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
                name_of(keys::LS_CUTOFF),
                keys::LS_CUTOFF,
                t(
                    "A hit below this overall score counts as no match: in the component list, Find substance and the batch.",
                    "Ein Treffer unter diesem Gesamt-Score gilt als kein Treffer: in der Komponentenliste, bei Substanz finden und im Batch.",
                ),
                num(&mut d.ls_cutoff, limit(keys::LS_CUTOFF)),
            );
            row(
                ui,
                name_of(keys::SQL_MAX_NUMBER_FRAGMENTS),
                keys::SQL_MAX_NUMBER_FRAGMENTS,
                t(
                    "NIST search: how many of the component's strongest fragments are compared.",
                    "NIST-Suche: wie viele der stärksten Fragmente der Komponente verglichen werden.",
                ),
                num(&mut d.sql_max_number_fragments, limit(keys::SQL_MAX_NUMBER_FRAGMENTS)),
            );
            row(
                ui,
                name_of(keys::SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS),
                keys::SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS,
                t(
                    "NIST search: a library entry is scored only if it shares more than this many of them.",
                    "NIST-Suche: ein Bibliothekseintrag wird nur bewertet, wenn er mehr als so viele davon teilt.",
                ),
                num(&mut d.sql_minimal_number_equal_fragments, limit(keys::SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS)),
            );
        });
        ranges(
            ui,
            keys::LS_MASS_FILTER,
            name_of(keys::LS_MASS_FILTER),
            &explain(
                t(
                    "Mass ranges (m/z from–to) ignored in both spectra.",
                    "Massenbereiche (m/z von–bis), die in beiden Spektren ignoriert werden.",
                ),
                keys::LS_MASS_FILTER,
            ),
            &mut d.ls_mass_filter,
        );

        ui.add_space(8.0);
        ui.strong(t("Centroiding", "Zentroidierung"));
        egui::Grid::new("centroid").num_columns(2).show(ui, |ui| {
            let used = t(
                "Finds the peaks of a profile spectrum when it is exported as centroids.",
                "Findet die Peaks eines Profilspektrums, wenn es als Zentroide exportiert wird.",
            );
            for (name, key, v) in [
                (
                    name_of(keys::SF_CENTROID_THRESHOLD_BEGIN),
                    keys::SF_CENTROID_THRESHOLD_BEGIN,
                    &mut d.sf_centroid_threshold_begin,
                ),
                (
                    name_of(keys::SF_CENTROID_THRESHOLD_END),
                    keys::SF_CENTROID_THRESHOLD_END,
                    &mut d.sf_centroid_threshold_end,
                ),
                (
                    name_of(keys::SF_CENTROID_MAX_BASE_DEV),
                    keys::SF_CENTROID_MAX_BASE_DEV,
                    &mut d.sf_centroid_max_base_dev,
                ),
            ] {
                row(ui, name, key, used, num(v, limit(key)));
            }
            row(
                ui,
                name_of(keys::SF_FWHM),
                keys::SF_FWHM,
                t(
                    "Points at or above this share of a peak's maximum make up its centroid mass, \
                     in every centroiding (high-resolution CDF import, centroid export, MSL export, TOF calibration).",
                    "Punkte ab diesem Anteil am Maximum eines Peaks bilden seine Zentroidmasse, \
                     bei jeder Zentroidierung (hochauflösender CDF-Import, Zentroid-Export, MSL-Export, TOF-Kalibrierung).",
                ),
                num(&mut d.sf_fwhm, limit(keys::SF_FWHM)).speed(0.01),
            );
        });

        ui.add_space(8.0);
        egui::CollapsingHeader::new(t(
            "Stored in the file, not used",
            "In der Datei gespeichert, nicht verwendet",
        ))
        .id_salt("Stored in the file, not used")
        .show(ui, |ui| {
            let unused = t(
                "Kept in the settings file; this version does not use it.",
                "In der Einstellungsdatei behalten; diese Version verwendet es nicht.",
            );
            // RE of 3.3 (docs/reverse-engineering/algorithms.md §4b.4): only saved and shown there too.
            let dead = t(
                "Kept in the settings file; 3.3 does not use it either.",
                "In der Einstellungsdatei behalten; auch 3.3 verwendet es nicht.",
            );
            egui::Grid::new("stored").num_columns(2).show(ui, |ui| {
                for (name, key, v) in [
                    (
                        name_of(keys::SF_PEAK_THRESHOLD_BEGIN),
                        keys::SF_PEAK_THRESHOLD_BEGIN,
                        &mut d.sf_peak_threshold_begin,
                    ),
                    (
                        name_of(keys::SF_PEAK_THRESHOLD_END),
                        keys::SF_PEAK_THRESHOLD_END,
                        &mut d.sf_peak_threshold_end,
                    ),
                    (
                        name_of(keys::SF_MAX_BASE_DEV),
                        keys::SF_MAX_BASE_DEV,
                        &mut d.sf_max_base_dev,
                    ),
                ] {
                    row(ui, name, key, dead, num(v, limit(key)));
                }
                let stored: [(&str, &str, &mut f64); 3] = [
                    (
                        name_of(keys::AN_MAX_BASE_DEV),
                        keys::AN_MAX_BASE_DEV,
                        &mut d.an_max_base_dev,
                    ),
                    (
                        name_of(keys::QUANT_MIN_DIST_IONS),
                        keys::QUANT_MIN_DIST_IONS,
                        &mut d.quant_min_dist_ions,
                    ),
                    (
                        name_of(keys::QUANT_QUALITY_CUTOFF),
                        keys::QUANT_QUALITY_CUTOFF,
                        &mut d.quant_quality_cutoff,
                    ),
                ];
                for (name, key, v) in stored {
                    row(ui, name, key, unused, num(v, limit(key)));
                }
                row(
                    ui,
                    name_of(keys::LS_RESOLUTION),
                    keys::LS_RESOLUTION,
                    unused,
                    num(&mut d.ls_resolution, limit(keys::LS_RESOLUTION)),
                );
                row(
                    ui,
                    name_of(keys::LS_USE_SCALED_LIB),
                    keys::LS_USE_SCALED_LIB,
                    dead,
                    egui::Checkbox::without_text(&mut d.ls_use_scaled_lib),
                );
                row(
                    ui,
                    name_of(keys::QUANT_ION_NUMBER),
                    keys::QUANT_ION_NUMBER,
                    unused,
                    egui::DragValue::new(&mut d.quant_ion_number)
                        .range(limit::<i32>(keys::QUANT_ION_NUMBER)),
                );
            });
        });
    }
}

/// The old file has no extension convention.
fn dialog() -> rfd::FileDialog {
    rfd::FileDialog::new()
}

/// Display text of a setting from the settings table, in the current language.
pub(crate) fn name_of(key: &str) -> &'static str {
    let (en, de) = md_core::settings::label(key);
    t(en, de)
}

/// The range of an input as the settings table has it (`md_core::settings::KEYS`), so window and
/// `validate` cannot disagree.
fn limit<N: egui::emath::Numeric>(key: &str) -> RangeInclusive<N> {
    bounds(limits(key))
}

fn bounds<N: egui::emath::Numeric>((lo, hi): (f64, f64)) -> RangeInclusive<N> {
    N::from_f64(lo)..=N::from_f64(hi)
}

/// What a setting does, then its key in the old 3.3 settings file.
fn explain(hint: &str, key: &str) -> String {
    match key {
        "" => hint.to_owned(),
        _ => tf!("{hint}\nOld key: {key}", "{hint}\nAlter Schlüssel: {key}"),
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
        &tf!(
            "Width in scans of the matched filter over the apex histogram. Auto: {AUTO_WIDTH_S} s at \
             each run's scan rate, at least 1 scan, at most 100 and a tenth of the run (3.3 always used 1). A width from a loaded or \
             imported 3.3 file is a set one; Save writes Auto as the width of the shown run (1 without one).",
            "Breite in Scans des Matched-Filters über dem Apex-Histogramm. Auto: {AUTO_WIDTH_S} s bei \
             der Scanrate jedes Laufs, mindestens 1 Scan, höchstens 100 und ein Zehntel des Laufs (3.3 nahm immer 1). Eine Breite \
             aus einer geladenen oder importierten 3.3-Datei ist fest gesetzt; Speichern schreibt Auto als Breite des gezeigten Laufs (1 ohne Lauf)."
        ),
        keys::AN_DECONVOLUTION_WIDTH,
    );
    let label = ui
        .label(name_of(keys::AN_DECONVOLUTION_WIDTH))
        .on_hover_text(&hint);
    ui.horizontal(|ui| {
        let mut auto = width.is_none();
        let check = ui.checkbox(&mut auto, t("Auto", "Auto"));
        if check.labelled_by(label.id).on_hover_text(&hint).changed() {
            // A set width starts from the automatic one.
            *width = (!auto).then(|| auto_width(times.unwrap_or_default()));
        }
        match width {
            Some(w) => ui
                .add(num(w, limit(keys::AN_DECONVOLUTION_WIDTH)))
                .labelled_by(label.id),
            None => match times {
                Some(s) => ui.label(tf!(
                    "{} scans in the shown run",
                    "{} Scans im gezeigten Lauf",
                    auto_width(s)
                )),
                None => ui.weak(tf!(
                    "{AUTO_WIDTH_S} s of each run",
                    "{AUTO_WIDTH_S} s jedes Laufs"
                )),
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
pub(crate) fn ranges(
    ui: &mut egui::Ui,
    key: &str,
    name: &str,
    hint: &str,
    list: &mut Vec<(f32, f32)>,
) {
    let bounds = limit::<f32>(key);
    ui.add_space(4.0);
    let list_name = ui.horizontal(|ui| {
        let list_name = ui.label(name).on_hover_text(hint);
        if ui
            .small_button("+")
            .on_hover_text(t("Add a range", "Einen Bereich hinzufügen"))
            .clicked()
        {
            list.push((0.0, 0.0));
        }
        list_name.id
    });
    let list_name = list_name.inner;
    let mut remove = None;
    for (i, (lo, hi)) in list.iter_mut().enumerate() {
        ui.horizontal(|ui| {
            let from = ui.label(t("from", "von"));
            ui.add(
                egui::DragValue::new(lo)
                    .range(bounds.clone())
                    .max_decimals(2),
            )
            .labelled_by(list_name)
            .labelled_by(from.id)
            .on_hover_text(hint);
            let to = ui.label(t("to", "bis"));
            ui.add(
                egui::DragValue::new(hi)
                    .range(bounds.clone())
                    .max_decimals(2),
            )
            .labelled_by(list_name)
            .labelled_by(to.id)
            .on_hover_text(hint);
            if ui
                .small_button("×")
                .on_hover_text(t("Remove this range", "Diesen Bereich entfernen"))
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

    #[test]
    fn a_cancelled_analysis_is_an_error_not_a_result() {
        let run = crate::results::tests::run();
        let cancel = Cancel::default();
        let key = Settings::default().key();
        assert!(analyze(&run.data, &key, &cancel).is_ok());
        cancel.cancel();
        assert_eq!(
            analyze(&run.data, &key, &cancel).err().as_deref(),
            Some("cancelled")
        );
    }

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
