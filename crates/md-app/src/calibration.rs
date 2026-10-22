//! RI calibration window: alkanes or a calibration mix from the library found in one run, or a `.calib`
//! file, in an editable RT → RI table. The calibration applies to every open run unless a run has its own;
//! it is optionally shifted per run onto a reference substance from the library (RT shift correction).
//! "Retention time only" is the old emergency exit without an alkane mix.

use std::collections::{BTreeSet, HashMap};
use std::path::PathBuf;
use std::sync::Arc;

use crate::lang::{Msg, msg, t, tf};
use crate::layout::ToolWindow;
use eframe::egui;
use md_core::deconv::Component;
use md_core::library::LibraryCompound;
use md_core::ri::{RiCalibration, marker_peaks, match_mix, series_points};
use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};

use crate::identify::Library;
use crate::{Job, Run};

/// One line of the calibration table; lines that are not `used` (a mix substance not found, or switched
/// off by hand) stay out of the calibration.
#[derive(Debug, Clone, PartialEq)]
struct Row {
    name: String,
    ri: f64,
    rt_min: f64,
    used: bool,
}

/// The calibration that RT-only alignment replaced.
struct Saved {
    cal: Option<Arc<RiCalibration>>,
    source: Msg,
    rows: Vec<Row>,
    applied: Vec<Row>,
    error: Option<Msg>,
}

type RowsJob = (Msg, Job<Result<Vec<Row>, Msg>>);

pub struct Calibration {
    /// Shared with the library search jobs; a new `Arc` on every change re-identifies all runs.
    pub cal: Option<Arc<RiCalibration>>,
    /// Where `cal` came from, shown in the window.
    source: Msg,
    pub open: bool,
    /// The table being edited and the state `cal` was last built from; an edit that is finished
    /// (no field dragged or focused) builds a new `cal`.
    rows: Vec<Row>,
    applied: Vec<Row>,
    /// Calibrations of single runs (by file), which win over `cal`.
    own: HashMap<PathBuf, Arc<RiCalibration>>,
    /// RT-only alignment: `cal` is the identity and the library is searched by its retention times.
    pub rt_only: bool,
    saved: Option<Saved>,
    /// Calibration mix: names of library compounds with RI, and the text filtering the picker.
    mix: BTreeSet<String>,
    mix_filter: String,
    mix_job: Option<RowsJob>,
    mz: f32,
    /// Alkane carbon numbers: C`first` … C`last` get RI 100 · n.
    first: u32,
    last: u32,
    /// Matched detector settings of the old calibration dialog (no defaults there; these find
    /// the markers of 20 Hz example runs).
    density: usize,
    threshold: i32,
    /// Library compound (by name) whose position shifts the calibration per run, and the minimum
    /// overall score it must reach there (old spin box; default `LS_CUTOFF`).
    pub reference: Option<String>,
    pub ref_cutoff: f64,
    /// Cutoff being edited; applied when the drag ends, so dragging does not search on every frame.
    cutoff_edit: f64,
    /// Running alkane search and the source it will set.
    job: Option<(Msg, Job<Result<RiCalibration, Msg>>)>,
    error: Option<Msg>,
}

impl Default for Calibration {
    fn default() -> Self {
        Self {
            cal: None,
            source: Msg::default(),
            open: false,
            rows: Vec::new(),
            applied: Vec::new(),
            own: HashMap::new(),
            rt_only: false,
            saved: None,
            mix: BTreeSet::new(),
            mix_filter: String::new(),
            mix_job: None,
            // Alkane fragment C4H9+.
            mz: 57.0,
            first: 10,
            last: 36,
            density: 40,
            threshold: 100,
            reference: None,
            ref_cutoff: 0.8,
            cutoff_edit: 0.8,
            job: None,
            error: None,
        }
    }
}

impl Calibration {
    /// The calibration of the run before any shift: its own, else the global one (RT only: the identity).
    pub fn base(&self, run: &Run) -> Option<&Arc<RiCalibration>> {
        let own = self.own.get(&run.path).filter(|_| !self.rt_only);
        own.or(self.cal.as_ref())
    }

    /// The run's calibration: shifted onto the reference when it was found there, else `base`.
    /// While a new search runs (library, calibration or reference changed) its shift is not known yet.
    /// None in RT-only mode, where "RI" is a time.
    pub fn for_run<'a>(&'a self, run: &'a Run) -> Option<&'a RiCalibration> {
        let ptr = |c: Option<&Arc<RiCalibration>>| c.map(Arc::as_ptr);
        let base = self.base(run).filter(|_| !self.rt_only);
        (run.ident.as_ref())
            .filter(|i| run.ident_job.is_none() && ptr(i.0.1.as_ref()) == ptr(base))
            .and_then(|i| i.2.as_ref())
            .or(base.map(|c| &**c))
    }

    /// What the batch and MID scans search with: the calibration, not while the library is rebuilt
    /// after a switch of RT-only (its anchors would not match), and the reference with its cutoff.
    pub fn for_search(
        &self,
        lib: Option<&Library>,
    ) -> (Option<&RiCalibration>, Option<(&str, f64)>) {
        let cal = (self.cal.as_deref())
            .filter(|_| lib.is_some_and(|l| l.prepared.is_rt_only() == self.rt_only));
        let reference = (self.reference.as_deref()).map(|n| (n, self.ref_cutoff));
        (cal, reference)
    }

    /// ", RI …" for the scan, or nothing without a calibration.
    pub fn ri_note(&self, run: &Run, scan: usize) -> String {
        self.for_run(run)
            .map(|c| format!(", RI {:.1}", c.ri(run.data.scan_times_s[scan])))
            .unwrap_or_default()
    }

    /// What a library compound is anchored to: its RI, or in RT-only mode its time in s. Negative: none.
    pub fn anchor(&self, c: &LibraryCompound) -> f64 {
        anchor(self.rt_only, c)
    }

    /// Gives one run its own calibration (one the old program saved with it).
    pub fn set_own(&mut self, run: PathBuf, cal: RiCalibration) {
        self.own.insert(run, Arc::new(cal));
    }

    /// Also drops a running search, so its late result cannot replace this calibration.
    fn set(&mut self, cal: RiCalibration, source: Msg) {
        self.rows = (cal.points().iter())
            .map(|(t, ri)| Row {
                name: String::new(),
                ri: *ri,
                rt_min: t / 60.0,
                used: true,
            })
            .collect();
        self.applied = self.rows.clone();
        (self.cal, self.source, self.error) = (Some(Arc::new(cal)), source, None);
        (self.job, self.mix_job) = (None, None);
    }

    /// Takes the table of a mix search; `cal` follows once the table is valid (see `apply_edits`).
    fn set_rows(&mut self, rows: Vec<Row>, source: Msg) {
        (self.rows, self.applied, self.source, self.error) = (rows, Vec::new(), source, None);
        (self.job, self.mix_job) = (None, None);
    }

    /// Builds `cal` from an edited table once no field is being dragged or typed in, so that a change
    /// takes effect at once but not on every frame of a drag. An invalid table keeps the old `cal`.
    fn apply_edits(&mut self) {
        if self.rows == self.applied {
            return;
        }
        // A time of 0 is a line just added by hand (no marker elutes then), not yet a point.
        let used = (self.rows.iter().filter(|r| r.used && r.rt_min > 0.0))
            .map(|r| (r.rt_min * 60.0, r.ri));
        match RiCalibration::new(used) {
            Ok(cal) => (self.cal, self.error) = (Some(Arc::new(cal)), None),
            Err(e) => self.error = Some(e.to_string().into()),
        }
        self.applied = self.rows.clone();
    }

    /// Adds an empty line by hand (markers without an alkane mix, e.g. FAMEs). It joins the calibration
    /// once it is valid, see `apply_edits`.
    fn add_point(&mut self) {
        if self.rows.is_empty() {
            self.source = msg!("entered by hand", "von Hand eingegeben");
        }
        self.rows.push(Row {
            name: String::new(),
            ri: 0.0,
            rt_min: 0.0,
            used: true,
        });
    }

    fn clear(&mut self) {
        (self.cal, self.job, self.mix_job) = (None, None, None);
        (self.rows, self.applied, self.rt_only, self.error) = (Vec::new(), Vec::new(), false, None);
        (self.own, self.saved) = (HashMap::new(), None);
    }

    /// Switches RT-only alignment. The identity maps time to "RI" (the library is searched by its
    /// retention times); the calibration in use before is kept and comes back when it is switched off.
    fn set_rt_only(&mut self, on: bool) {
        self.rt_only = on;
        if on {
            let saved = Saved {
                cal: self.cal.take(),
                source: std::mem::take(&mut self.source),
                rows: std::mem::take(&mut self.rows),
                applied: std::mem::take(&mut self.applied),
                error: self.error.take(),
            };
            self.set(
                RiCalibration::identity(),
                msg!("retention time only", "nur Retentionszeit"),
            );
            self.saved = Some(saved);
        } else if let Some(s) = self.saved.take() {
            (self.cal, self.source, self.rows) = (s.cal, s.source, s.rows);
            (self.applied, self.error) = (s.applied, s.error);
        }
    }

    pub fn window(
        &mut self,
        ctx: &egui::Context,
        runs: &[Run],
        active: usize,
        lib: Option<&Library>,
    ) {
        if let Some((source, job)) = &self.job
            && let Some(result) = job
                .poll(|| t("RI calibration", "RI-Kalibrierung").into())
                .map(|r| r.and_then(|x| x))
        {
            match result {
                Ok(cal) => self.set(cal, source.clone()),
                Err(e) => self.error = Some(e),
            }
            self.job = None;
        }
        if let Some((source, job)) = &self.mix_job
            && let Some(result) = job
                .poll(|| t("calibration mix", "Kalibriermischung").into())
                .map(|r| r.and_then(|x| x))
        {
            match result {
                Ok(rows) => self.set_rows(rows, source.clone()),
                Err(e) => self.error = Some(e),
            }
            self.mix_job = None;
        }
        let mut open = self.open;
        egui::Window::new(t("RI calibration", "RI-Kalibrierung"))
            .id(egui::Id::unique("RI calibration"))
            .tool_window(ctx, 1, [380.0, 420.0])
            .open(&mut open)
            .show(ctx, |ui| {
                self.contents(ui, runs.get(active));
                self.mix_contents(ui, runs.get(active), lib);
                self.table(ui);
                self.run_contents(ui, runs);
                self.shift_contents(ui, runs, lib);
            });
        self.open = open;
    }

    /// Why the last search, load or table edit failed.
    pub fn error_label(&self, ui: &mut egui::Ui) {
        if let Some(e) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, &**e);
        }
    }

    pub fn contents(&mut self, ui: &mut egui::Ui, run: Option<&Run>) {
        // First, so the guided step that finds the alkanes shows why it failed too.
        self.error_label(ui);
        let was = self.rt_only;
        ui.checkbox(
            &mut self.rt_only,
            t(
                "Use retention time for alignment (not recommended)",
                "Retentionszeit zum Abgleich verwenden (nicht empfohlen)",
            ),
        )
        .on_hover_text(t(
            "Only if you have no n-alkane mix. Alignment and identification by retention time \
             are not reproducible and the results cannot be compared between instruments and labs.",
            "Nur ohne n-Alkan-Mischung. Abgleich und Identifizierung über die Retentionszeit sind \
             nicht reproduzierbar, und die Ergebnisse lassen sich nicht zwischen Geräten und Laboren vergleichen.",
        ));
        if self.rt_only != was {
            self.set_rt_only(self.rt_only);
        }
        if self.rt_only {
            ui.colored_label(
                ui.visuals().warn_fg_color,
                t(
                    "Retention times are compared directly. Library entries need a retention time; \
                     components get no retention index.",
                    "Retentionszeiten werden direkt verglichen. Bibliothekseinträge brauchen eine \
                     Retentionszeit; Komponenten erhalten keinen Retentionsindex.",
                ),
            );
            return;
        }
        egui::Grid::new("alkane settings").show(ui, |ui| {
            let label = ui.label(t("Marker ion m/z", "Marker-Ion m/z"));
            ui.add(egui::DragValue::new(&mut self.mz).range(1.0..=2000.0))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "The ion whose peaks mark the alkanes, e.g. 57 for n-alkanes",
                    "Das Ion, dessen Peaks die Alkane markieren, z. B. 57 für n-Alkane",
                ));
            ui.end_row();
            let label = ui.label(t("Alkanes C", "Alkane C"));
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut self.first).range(1..=self.last))
                    .labelled_by(label.id)
                    .on_hover_text(t(
                        "Carbon number of the first alkane in the mix",
                        "Kohlenstoffzahl des ersten Alkans der Mischung",
                    ));
                let to = ui.label(t("to C", "bis C"));
                ui.add(egui::DragValue::new(&mut self.last).range(self.first..=100))
                    .labelled_by(to.id)
                    .on_hover_text(t(
                        "Carbon number of the last alkane in the mix",
                        "Kohlenstoffzahl des letzten Alkans der Mischung",
                    ));
            });
            ui.end_row();
        });
        // The matched detector's settings for the marker peaks.
        crate::guided::advanced(ui, "alkane filter", |ui| {
            egui::Grid::new("alkane filter").show(ui, |ui| {
            let label = ui.label(t("Filter width (scans)", "Filterbreite (Scans)"));
            ui.add(egui::DragValue::new(&mut self.density).range(1..=1000))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Peak width the filter is tuned to (sigma about width / 4 scans); match it to the alkane peaks",
                    "Peakbreite, auf die der Filter abgestimmt ist (Sigma etwa Breite / 4 Scans); an die Alkan-Peaks anpassen",
                ));
            ui.end_row();
            let label = ui.label(t("Threshold", "Schwelle"));
            ui.add(egui::DragValue::new(&mut self.threshold).range(0..=i32::MAX))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Peaks are kept where the filtered marker-ion trace rises above this value (not the raw intensity)",
                    "Peaks bleiben, wo die gefilterte Marker-Ion-Spur diesen Wert übersteigt (nicht die Rohintensität)",
                ));
            ui.end_row();
            });
        });
        ui.horizontal(|ui| {
            if run.is_none() {
                ui.weak(t(
                    "Open a run to find the alkanes in it.",
                    "Einen Lauf öffnen, um die Alkane darin zu finden.",
                ));
            }
            if let Some(run) = run
                && ui
                    .add_enabled(self.job.is_none(), egui::Button::new(tf!("Find alkanes in {}", "Alkane in {} finden", run.name)))
                    .on_hover_text(t(
                        "The highest peaks of the marker ion, in time order, become C first … C last",
                        "Die höchsten Peaks des Marker-Ions werden in zeitlicher Reihenfolge zu C erstes … C letztes",
                    ))
                    .clicked()
            {
                let (data, mz, density, threshold) = (run.data.clone(), self.mz, self.density, self.threshold);
                let ris: Vec<f64> = (self.first..=self.last).map(|n| 100.0 * n as f64).collect();
                let source = Msg::from(format!("{}, m/z {mz}, C{}–C{}", run.name, self.first, self.last));
                let job = Job::spawn(ui.ctx(), move |_, cancel| {
                    let peaks = marker_peaks(&data, mz, density, threshold);
                    if cancel.is_cancelled() {
                        return Err(msg!("Cancelled.", "Abgebrochen."));
                    }
                    RiCalibration::new(series_points(&peaks, &ris)).map_err(|e| e.to_string().into())
                });
                self.job = Some((source, job));
            }
            if self.job.is_some() {
                ui.spinner();
                if ui.small_button(t("Cancel", "Abbrechen")).clicked() {
                    self.job = None;
                }
            }
            if ui.button(t("Load .calib…", ".calib laden…")).clicked()
                && let Some(path) = calib_dialog().pick_file()
            {
                match RiCalibration::read_calib(&path) {
                    Ok(cal) => self.set(cal, path.display().to_string().into()),
                    Err(e) => self.error = Some(format!("{}: {e}", path.display()).into()),
                }
            }
            if let Some(cal) = &self.cal {
                if ui.button(t("Save .calib…", ".calib speichern…")).clicked()
                    && let Some(path) = calib_dialog().save_file()
                    && let Err(e) = std::fs::write(&path, cal.to_calib())
                {
                    self.error = Some(format!("{}: {e}", path.display()).into());
                }
                if ui.button(t("Clear", "Leeren")).clicked() {
                    self.clear();
                }
            }
        });
    }

    /// Calibration mix (old wizard page "RI Calibration Mix"): library compounds with RI that are found
    /// in a run by their spectra; their times in that run and their library RIs make the table.
    fn mix_contents(&mut self, ui: &mut egui::Ui, run: Option<&Run>, lib: Option<&Library>) {
        let Some(lib) = lib.filter(|_| !self.rt_only) else {
            return;
        };
        egui::CollapsingHeader::new(t(
            "Calibration mix from the library",
            "Kalibriermischung aus der Bibliothek",
        ))
        .id_salt("Calibration mix from the library")
        .show(ui, |ui| {
            ui.horizontal(|ui| {
                let add = ui.label(t("Add", "Hinzufügen"));
                ui.text_edit_singleline(&mut self.mix_filter)
                    .labelled_by(add.id)
                    .on_hover_text(t(
                        "Part of a library compound name; only compounds with an RI are offered",
                        "Teil eines Substanznamens der Bibliothek; nur Substanzen mit RI werden angeboten",
                    ));
            });
            if self.mix.is_empty() && self.mix_filter.is_empty() {
                ui.weak(t(
                    "Type a name to pick the substances of your calibration mix.",
                    "Einen Namen eingeben, um die Substanzen Ihrer Kalibriermischung zu wählen.",
                ));
            }
            let filter = self.mix_filter.to_lowercase();
            if !filter.is_empty() {
                let found: Vec<_> = (lib.compounds.iter())
                    .filter(|c| c.spectrum.ri >= 0.0 && !self.mix.contains(&c.name))
                    .filter(|c| c.name.to_lowercase().contains(&filter))
                    .take(8)
                    .collect();
                for c in found {
                    if ui
                        .small_button(format!("+ {} (RI {:.1})", c.name, c.spectrum.ri))
                        .clicked()
                    {
                        self.mix.insert(c.name.clone());
                    }
                }
            }
            let mut remove = None;
            for name in &self.mix {
                if ui.small_button(format!("− {name}")).clicked() {
                    remove = Some(name.clone());
                }
            }
            if let Some(name) = remove {
                self.mix.remove(&name);
            }
            let Some(run) = run else {
                return;
            };
            let button = egui::Button::new(tf!("Find the mix in {}", "Mischung in {} finden", run.name));
            let searching = self.mix_job.is_some();
            if ui
                .add_enabled(!self.mix.is_empty() && !searching, button)
                .on_hover_text(t(
                    "Going through the mix by RI, each substance takes the next component in time \
                     order that the library search finds for it",
                    "Die Mischung wird nach RI durchgegangen; jede Substanz erhält die nächste Komponente \
                     in zeitlicher Reihenfolge, die die Bibliothekssuche für sie findet",
                ))
                .clicked()
            {
                let mix = (lib.compounds.iter())
                    .filter(|c| c.spectrum.ri >= 0.0 && self.mix.contains(&c.name))
                    .cloned()
                    .collect();
                let (comps, data, s) = (
                    run.components.clone(),
                    run.data.clone(),
                    lib.settings.clone(),
                );
                let source = msg!(
                    "{}, {} mix substances",
                    "{}, {} Substanzen der Mischung",
                    run.name,
                    self.mix.len()
                );
                let job = Job::spawn(ui.ctx(), move |_, cancel| {
                    match cancel.is_cancelled() {
                        true => Err(msg!("Cancelled.", "Abgebrochen.")),
                        false => Ok(find_mix(mix, &comps, &data.scan_times_s, &s)),
                    }
                });
                self.mix_job = Some((source, job));
            }
            if searching {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(t("Searching the mix…", "Mischung wird gesucht…"));
                    if ui.small_button(t("Cancel", "Abbrechen")).clicked() {
                        self.mix_job = None;
                    }
                });
            }
        });
    }

    /// The editable table: a change takes effect when the field is left or the drag ends.
    pub fn table(&mut self, ui: &mut egui::Ui) {
        if self.rt_only {
            return;
        }
        if self.rows.is_empty() {
            ui.label(t(
                "Not calibrated: components have no RI, the library search uses the spectrum only.",
                "Nicht kalibriert: Komponenten haben keinen RI, die Bibliothekssuche nutzt nur das Spektrum.",
            ));
            self.add_point_button(ui);
            return;
        }
        ui.label(tf!(
            "{} of {} points from {}",
            "{} von {} Punkten aus {}",
            self.cal.as_ref().map_or(0, |c| c.points().len()),
            self.rows.len(),
            self.source
        ));
        let mut busy = false;
        egui::ScrollArea::vertical()
            .max_height(300.0)
            .show(ui, |ui| {
                egui::Grid::new("calibration points")
                    .striped(true)
                    .show(ui, |ui| {
                        let use_ = ui.strong(t("Use", "Verwenden"));
                        ui.strong(t("Substance", "Substanz"));
                        let rt = ui.strong(t("RT (min)", "RT (min)"));
                        let ri = ui.strong("RI");
                        ui.end_row();
                        for r in &mut self.rows {
                            ui.checkbox(&mut r.used, "")
                                .labelled_by(use_.id)
                                .on_hover_text(t(
                                    "Use this point for the calibration",
                                    "Diesen Punkt für die Kalibrierung verwenden",
                                ));
                            // A mix substance that was not found has no time yet: red until it gets one.
                            let name = egui::RichText::new(&r.name);
                            let found = r.used || r.rt_min > 0.0;
                            let text = match found {
                                true => name,
                                false => name.color(ui.visuals().error_fg_color),
                            };
                            let label = ui.label(text);
                            if !found {
                                label.on_hover_text(t(
                                    "Not found in the run: enter its time or leave it unused",
                                    "Im Lauf nicht gefunden: Zeit eingeben oder ungenutzt lassen",
                                ));
                            }
                            let time = egui::DragValue::new(&mut r.rt_min)
                                .speed(0.001)
                                .range(0.0..=f64::MAX)
                                .fixed_decimals(3);
                            let time = ui
                                .add(time)
                                .labelled_by(rt.id)
                                .on_hover_text(t(
                                    "Retention time of the substance in this run, in minutes",
                                    "Retentionszeit der Substanz in diesem Lauf, in Minuten",
                                ));
                            let n = ui
                                .add(egui::DragValue::new(&mut r.ri).fixed_decimals(1))
                                .labelled_by(ri.id)
                                .on_hover_text(t(
                                    "Retention index of the substance (100 per carbon for n-alkanes)",
                                    "Retentionsindex der Substanz (100 je Kohlenstoff bei n-Alkanen)",
                                ));
                            busy |= [time, n].iter().any(|x| x.dragged() || x.has_focus());
                            ui.end_row();
                        }
                    });
            });
        self.add_point_button(ui);
        if !busy {
            self.apply_edits();
        }
    }

    fn add_point_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .button(t("Add point", "Punkt hinzufügen"))
            .on_hover_text(t(
                "A new line for a marker you know (e.g. FAMEs): its time in this run and its RI. \
                 At least two points with rising times and RIs make a calibration.",
                "Eine neue Zeile für einen bekannten Marker (z. B. FAMEs): seine Zeit in diesem Lauf \
                 und sein RI. Mindestens zwei Punkte mit steigenden Zeiten und RIs ergeben eine Kalibrierung.",
            ))
            .clicked()
        {
            self.add_point();
        }
    }

    /// Calibration per run (old "Calibration Setter"): which runs use their own, and its `.calib`.
    fn run_contents(&mut self, ui: &mut egui::Ui, runs: &[Run]) {
        if self.rt_only || runs.is_empty() || (self.cal.is_none() && self.own.is_empty()) {
            return;
        }
        ui.separator();
        ui.strong(t("Calibration per run", "Kalibrierung je Lauf"));
        egui::Grid::new("run calibrations")
            .striped(true)
            .show(ui, |ui| {
                for run in runs {
                    ui.label(&run.name);
                    let own = self.own.get(&run.path).cloned();
                    ui.label(match &own {
                        Some(c) => tf!("own, {} points", "eigene, {} Punkte", c.points().len()),
                        None => t("table above", "Tabelle oben").into(),
                    });
                    ui.horizontal(|ui| {
                        if let Some(cal) = &self.cal
                            && ui
                                .small_button(t("Use table", "Tabelle verwenden"))
                                .clicked()
                        {
                            self.own.insert(run.path.clone(), cal.clone());
                        }
                        if own.is_some() && ui.small_button(t("Reset", "Zurücksetzen")).clicked() {
                            self.own.remove(&run.path);
                        }
                        if let Some(cal) = own.as_ref().or(self.cal.as_ref())
                            && ui
                                .small_button(t("Save .calib…", ".calib speichern…"))
                                .clicked()
                            && let Some(path) = calib_dialog()
                                .set_file_name(format!("{}.calib", run.name))
                                .save_file()
                            && let Err(e) = std::fs::write(&path, cal.to_calib())
                        {
                            self.error = Some(format!("{}: {e}", path.display()).into());
                        }
                    });
                    ui.end_row();
                }
            });
    }

    /// RT shift correction (§6.5): reference choice and the shift found in each run.
    fn shift_contents(&mut self, ui: &mut egui::Ui, runs: &[Run], lib: Option<&Library>) {
        if self.cal.is_none() && self.own.is_empty() {
            return;
        }
        ui.separator();
        let Some(lib) = lib else {
            ui.label(t(
                "Open a library to correct RT shifts with a reference substance.",
                "Eine Bibliothek öffnen, um RT-Verschiebungen mit einer Referenzsubstanz zu korrigieren.",
            ));
            return;
        };
        let rt_only = self.rt_only;
        ui.horizontal(|ui| {
            let reference = ui.label(t("Reference substance", "Referenzsubstanz"));
            let shown = self.reference.as_deref().unwrap_or(t("None", "Keine"));
            egui::ComboBox::from_id_salt("reference")
                .selected_text(shown)
                .show_ui(ui, |ui| {
                    ui.selectable_value(&mut self.reference, None, t("None", "Keine"));
                    // Only compounds with a library RI (a time in RT-only mode) can anchor the shift.
                    for c in lib.compounds.iter().filter(|c| anchor(rt_only, c) >= 0.0) {
                        let label = format!("{} ({:.1})", c.name, anchor(rt_only, c));
                        ui.selectable_value(&mut self.reference, Some(c.name.clone()), label);
                    }
                })
            .response
            .labelled_by(reference.id)
            .on_hover_text(t(
                "A substance found in every run; its position there corrects the retention-time shift of the run",
                "Eine Substanz, die in jedem Lauf gefunden wird; ihre Lage dort korrigiert die Retentionszeit-Verschiebung des Laufs",
            ));
            let label = ui.label(t("Cutoff", "Grenzwert"));
            let edit = ui
                .add(
                    egui::DragValue::new(&mut self.cutoff_edit)
                        .range(0.0..=1.0)
                        .speed(0.01),
                )
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Lowest match score (0–1) at which the reference counts as found in a run",
                    "Niedrigster Treffer-Score (0–1), ab dem die Referenz in einem Lauf als gefunden gilt",
                ));
            if !edit.dragged() {
                self.ref_cutoff = self.cutoff_edit;
            }
        });
        if self.reference.is_none() {
            return;
        }
        egui::Grid::new("shifts").striped(true).show(ui, |ui| {
            ui.strong(t("Run", "Lauf"));
            ui.strong(t("RT shift (s)", "RT-Verschiebung (s)"));
            ui.end_row();
            for run in runs {
                ui.label(&run.name);
                let shift = match (&run.ident, &run.ident_job) {
                    (_, Some(_)) | (None, None) => t("searching…", "wird gesucht…").into(),
                    (Some((_, _, Some(c))), None) => format!("{:+.2}", c.shift_s),
                    (Some(_), None) => t(
                        "reference not found, unshifted",
                        "Referenz nicht gefunden, unverschoben",
                    )
                    .into(),
                };
                ui.label(shift);
                ui.end_row();
            }
        });
    }
}

/// Mix substances (sorted by RI) found in a run: the library search with spectra only, then the old
/// assignment in time order (`match_mix`). Not found: unused row.
fn find_mix(
    mut mix: Vec<LibraryCompound>,
    comps: &[Component],
    scan_times_s: &[f64],
    s: &SearchSettings,
) -> Vec<Row> {
    mix.sort_by(|a, b| a.spectrum.ri.total_cmp(&b.spectrum.ri));
    // The mix is found by its spectra alone, whatever the search settings say about RI.
    let s = &SearchSettings {
        use_ri: false,
        use_spec: true,
        ..s.clone()
    };
    let id = map_compounds_to_library(comps, None, &PreparedLibrary::new(&mix, s), s);
    let rts: Vec<f64> = (comps.iter())
        .map(|c| scan_times_s[c.apex_scan as usize])
        .collect();
    let ris: Vec<f64> = mix.iter().map(|c| c.spectrum.ri).collect();
    let score = |c: usize, m: usize| {
        (id.hits[c].iter())
            .find(|h| h.compound == m)
            .map_or(0.0, |h| h.overall)
    };
    let found = match_mix(&ris, &rts, 0.0, score);
    mix.iter()
        .zip(found)
        .map(|(c, k)| Row {
            name: c.name.clone(),
            ri: c.spectrum.ri,
            rt_min: k.map_or(0.0, |k| rts[k] / 60.0),
            used: k.is_some(),
        })
        .collect()
}

fn anchor(rt_only: bool, c: &LibraryCompound) -> f64 {
    match rt_only {
        true => c.spectrum.rt_ms / 1000.0,
        false => c.spectrum.ri,
    }
}

fn calib_dialog() -> rfd::FileDialog {
    rfd::FileDialog::new().add_filter(t("RI calibration", "RI-Kalibrierung"), &["calib"])
}

#[cfg(test)]
mod tests {
    use super::*;

    fn substance(name: &str, ri: f64, mass: f32) -> LibraryCompound {
        let mut c = LibraryCompound::new(name);
        c.spectrum.ri = ri;
        (c.spectrum.masses, c.spectrum.values) = (vec![mass], vec![1000]);
        c
    }

    fn component(apex: u32, mass: u32) -> Component {
        Component {
            first_scan: apex,
            apex_scan: apex,
            last_scan: apex,
            model: Vec::new(),
            spectrum: vec![(mass, 1000)],
            peaks: Vec::new(),
            exact_spectrum: Vec::new(),
        }
    }

    #[test]
    fn rt_only_off_brings_back_the_previous_calibration() {
        let cal = RiCalibration::new([(60.0, 1000.0), (120.0, 1100.0)]).unwrap();
        let mut c = Calibration::default();
        c.set(cal.clone(), "x.calib".into());
        let before = (c.cal.clone(), c.rows.clone());
        c.set_rt_only(true);
        assert_eq!(&*c.source, "retention time only");
        c.set_rt_only(false);
        assert_eq!((c.cal.clone(), c.rows.clone()), before);
        assert_eq!(&*c.source, "x.calib");
        assert_eq!(c.cal.as_deref(), Some(&cal));
    }

    #[test]
    fn points_entered_by_hand_make_a_calibration() {
        // FAME markers of the example run 7235eg08 (RT in s, RI).
        let fames = [
            (252.46, 262_320.0),
            (311.36, 323_120.0),
            (369.61, 381_020.0),
        ];
        let mut c = Calibration::default();
        for (rt, ri) in fames {
            c.add_point();
            // The empty line the button adds (time 0) does not disturb the calibration so far.
            let before = c.cal.clone();
            c.apply_edits();
            assert_eq!(c.cal, before);
            let r = c.rows.last_mut().unwrap();
            (r.rt_min, r.ri) = (rt / 60.0, ri);
            // One point is not a calibration yet; the table keeps its error and no `cal`.
            c.apply_edits();
            assert_eq!(c.cal.is_some(), c.rows.len() > 1);
        }
        assert_eq!(&*c.source, "entered by hand");
        assert!(c.error.is_none());
        let cal = c.cal.as_ref().unwrap();
        for (rt, ri) in fames {
            assert!((cal.ri(rt) - ri).abs() < 1e-6);
        }
        // Switching a point off recalibrates from the rest.
        c.rows[2].used = false;
        c.apply_edits();
        assert_eq!(c.cal.as_ref().unwrap().points().len(), 2);
    }

    #[test]
    fn source_follows_a_language_switch() {
        let mut c = Calibration::default();
        c.add_point();
        c.rows[0] = Row {
            name: String::new(),
            ri: 1000.0,
            rt_min: 1.0,
            used: true,
        };
        c.apply_edits();
        let source = c.source.clone();
        crate::lang::set(crate::lang::Lang::De);
        assert_eq!(&*c.source, "von Hand eingegeben");
        crate::lang::set(crate::lang::Lang::En);
        assert_eq!(&*source, "entered by hand");
    }

    #[test]
    fn rt_only_toggle_keeps_the_error_of_an_unfinished_table() {
        let mut c = Calibration::default();
        c.add_point();
        c.rows[0] = Row {
            name: String::new(),
            ri: 1000.0,
            rt_min: 1.0,
            used: true,
        };
        c.apply_edits();
        let error = c.error.clone();
        assert!(error.is_some());
        c.set_rt_only(true);
        c.set_rt_only(false);
        assert_eq!(c.error, error);
        assert!(c.cal.is_none());
    }

    #[test]
    fn rt_only_off_without_calibration_stays_uncalibrated() {
        let mut c = Calibration::default();
        c.set_rt_only(true);
        assert!(c.cal.is_some());
        c.set_rt_only(false);
        assert!(c.cal.is_none() && c.rows.is_empty() && c.source.is_empty() && !c.rt_only);
    }

    #[test]
    fn mix_rows_follow_ri_and_leave_missing_substances_unused() {
        // C12 (RI 1200) and C10 (RI 1000) are in the run, C11 (RI 1100) is not; the library is unordered.
        let mix = vec![
            substance("C12", 1200.0, 60.0),
            substance("C11", 1100.0, 50.0),
            substance("C10", 1000.0, 40.0),
        ];
        let comps = [component(1, 40), component(3, 60)];
        let times = [0.0, 60.0, 90.0, 120.0];
        // Whatever the settings say about RI and spectra, the mix is found by its spectra.
        let odd = SearchSettings {
            use_ri: true,
            use_spec: false,
            ..SearchSettings::default()
        };
        let rows = find_mix(mix, &comps, &times, &odd);
        let rows: Vec<_> = rows
            .iter()
            .map(|r| (&*r.name, r.ri, r.rt_min, r.used))
            .collect();
        assert_eq!(
            rows,
            [
                ("C10", 1000.0, 1.0, true),
                ("C11", 1100.0, 0.0, false),
                ("C12", 1200.0, 2.0, true)
            ]
        );
    }
}
