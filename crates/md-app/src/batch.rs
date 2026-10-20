//! Batch quantification window: scans many runs on a worker thread against the open library and
//! shows the substance × run table (`md_core::batch`). Filters and groups apply without a rescan.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use crate::layout::ToolWindow;
use crate::settings::row;
use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::batch::{
    Completed, Found, RunResult, ScanSettings, ShiftReference, Table, TableSettings, complete_scan,
    scan_runs, table,
};
use md_core::batch::{cell_ions, hit_list, ions_path, settings_path};
use md_core::library::LibraryCompound;
use md_core::mid::frag_export;
use md_core::normalize::{Normalization, POOL_GROUP, Pool, apply};
use md_core::ri::RiCalibration;
use md_core::settings::GcmsSettings;
use md_core::stats::{group_values, student, summary};

use crate::identify::Library;
use crate::spawn;

type Scanned = Option<(Vec<RunResult>, Vec<String>, Completed)>;

struct Job {
    total: usize,
    done: Arc<AtomicUsize>,
    cancel: Arc<AtomicBool>,
    started: Instant,
    lib: Arc<Library>,
    rx: mpsc::Receiver<Scanned>,
}

struct Scan {
    lib: Arc<Library>,
    runs: Arc<Vec<RunResult>>,
    /// Unknown substances pooled by an untargeted scan; they follow the library's entries.
    extra: Arc<Vec<LibraryCompound>>,
    rescanned: usize,
    errors: Vec<String>,
    took: Duration,
}

/// Which normalisation (F13) the table shows; the table is derived again when this changes.
#[derive(Clone, PartialEq)]
struct Norm {
    kind: Kind,
    /// Library index of the internal standard.
    standard: Option<usize>,
    pool_group: String,
    pool_nearest: usize,
    pool_median: bool,
    z_score: bool,
    /// Runs and substances the user removed; part of the key because they change the shown table.
    edits: Edits,
}

/// Columns and rows removed or moved in the shown table. Runs are edited before the
/// normalisation (a removed run takes no part in z-scores), substances after it.
#[derive(Clone, Default, PartialEq)]
struct Edits {
    /// Columns of the table `table` built that are kept, in display order; `None` = all.
    cols: Option<Vec<usize>>,
    /// Library indices of the removed substances.
    rows: Vec<usize>,
}

impl Edits {
    fn cols(&self, n: usize) -> Vec<usize> {
        self.cols.clone().unwrap_or_else(|| (0..n).collect())
    }

    fn is_empty(&self) -> bool {
        self.cols.is_none() && self.rows.is_empty()
    }

    /// "Restore removed": brings back all runs and substances and the original order.
    fn restore(&mut self) {
        *self = Self::default();
    }

    /// Applies what the user asked for in the shown table with `n` columns in the base table.
    fn apply(&mut self, act: Act, n: usize) {
        let mut cols = self.cols(n);
        match act {
            Act::DropRun(j) if j < cols.len() => {
                cols.remove(j);
            }
            Act::MoveRun(j, to) if j < cols.len() && to < cols.len() => cols.swap(j, to),
            Act::DropRow(compound) => return self.rows.push(compound),
            _ => return,
        }
        self.cols = Some(cols);
    }
}

/// What a click in the shown table asks for; columns are numbered as shown.
#[derive(Clone, Copy)]
enum Act {
    DropRun(usize),
    MoveRun(usize, usize),
    DropRow(usize),
}

/// What the shown table lists.
#[derive(Clone, Copy, PartialEq, Default)]
enum Mode {
    /// Summed quant ion integrals per run.
    #[default]
    Values,
    /// Retention time of the hit per run, in minutes.
    Times,
    /// Mean, relative standard error, hits and p-value against the blanks, per replicate group.
    Groups,
}

/// How the shown table is displayed; changing it does not derive the table again.
struct View {
    mode: Mode,
    totals: bool,
    /// A line per quant ion under each substance (the old tree), with the raw ion integrals.
    ion_rows: bool,
    /// Also write the integral of every quant ion when exporting.
    export_ions: bool,
    /// Replicate group of the blank runs (the old program's fixed name).
    blanks: String,
}

impl Default for View {
    fn default() -> Self {
        Self {
            mode: Mode::default(),
            totals: false,
            ion_rows: false,
            export_ions: false,
            blanks: "BLANK".into(),
        }
    }
}

#[derive(Clone, Copy, PartialEq)]
enum Kind {
    None,
    InternalStandard,
    TotalSignal,
    Pool,
}

impl Default for Norm {
    fn default() -> Self {
        Self {
            kind: Kind::None,
            standard: None,
            pool_group: POOL_GROUP.into(),
            pool_nearest: 2,
            pool_median: false,
            z_score: false,
            edits: Edits::default(),
        }
    }
}

/// Options of the scan itself: changing them needs a new run of the batch.
#[derive(Clone, PartialEq)]
struct ScanOptions {
    untargeted: bool,
    rescan: bool,
}

pub struct Batch {
    pub open: bool,
    paths: Vec<PathBuf>,
    options: ScanOptions,
    /// Runs in memory at once; each run is analysed in parallel anyway.
    parallel: usize,
    pub settings: TableSettings,
    /// The applied settings (Analysis → Settings): scan settings and, when they change, the quant ions.
    applied: GcmsSettings,
    /// Whether the search is by retention time only (Analysis → RI calibration).
    rt_only: bool,
    job: Option<Job>,
    scan: Option<Scan>,
    /// Table of `scan` for these settings, or the pattern error; rebuilt when the settings change.
    table: Option<(TableSettings, Result<Table, String>)>,
    norm: Norm,
    view: View,
    /// Total intensities of the shown table (`generation`), when they are displayed.
    totals: (u64, Vec<f64>),
    /// `table` normalised as `norm` says.
    shown: Option<(Norm, Table)>,
    /// Counts the tables `shown` has held, so windows reading it know when it changed.
    generation: u64,
    /// File dates of the table's runs and how many could not be read; read once per table.
    times: Option<(Vec<i64>, usize)>,
    note: String,
    /// FRAG export being written; its note replaces `note` when done.
    frag_job: Option<mpsc::Receiver<String>>,
}

impl Default for Batch {
    fn default() -> Self {
        Self {
            open: false,
            paths: Vec::new(),
            // 15 example runs: 1.1 s with 1, 0.17 s with 4 (90 MiB), 0.10 s with 15 (189 MiB).
            parallel: 4,
            options: ScanOptions {
                untargeted: false,
                rescan: false,
            },
            settings: TableSettings::default(),
            applied: GcmsSettings::default(),
            rt_only: false,
            job: None,
            scan: None,
            table: None,
            norm: Norm::default(),
            view: View::default(),
            totals: (0, Vec::new()),
            shown: None,
            generation: 0,
            times: None,
            note: String::new(),
            frag_job: None,
        }
    }
}

impl Batch {
    /// Takes the applied settings; changed `QUANT_*` values replace those of the window.
    pub fn use_settings(&mut self, applied: &GcmsSettings, rt_only: bool) {
        self.rt_only = rt_only;
        if self.applied != *applied {
            let quant = applied.quant();
            if quant != self.applied.quant() {
                self.settings.quant = quant;
            }
            self.applied = applied.clone();
        }
    }

    /// The table as shown (normalised) with a number that changes whenever it is rebuilt.
    pub fn shown_table(&self) -> Option<(u64, &Table)> {
        self.shown.as_ref().map(|s| (self.generation, &s.1))
    }

    /// A finished scan, as if the batch had run; for tests that draw the result table.
    #[cfg(test)]
    pub fn finished(&mut self, lib: Arc<Library>, runs: Vec<RunResult>) {
        self.scan = Some(Scan {
            lib,
            runs: Arc::new(runs),
            extra: Arc::default(),
            rescanned: 0,
            errors: Vec::new(),
            took: Duration::ZERO,
        });
    }

    /// Whether the shown table is z-scored; a missing value cannot count as 0 there.
    pub fn z_scored(&self) -> bool {
        self.shown.as_ref().is_some_and(|s| s.0.z_score)
    }

    pub fn window(
        &mut self,
        ctx: &egui::Context,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        self.poll(ctx);
        let mut open = self.open;
        egui::Window::new("Batch quantification")
            .tool_window(ctx, 2, [900.0, 600.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui, lib, cal, reference));
        self.open = open;
    }

    fn poll(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.frag_job
            && let Some(r) = crate::poll_job(rx, "MID export")
        {
            (self.note, self.frag_job) = (r.unwrap_or_else(|e| e), None);
        }
        let Some(job) = &self.job else { return };
        match job.rx.try_recv() {
            Ok(Some((runs, mut errors, done))) => {
                errors.extend(done.errors);
                self.scan = Some(Scan {
                    lib: job.lib.clone(),
                    runs: Arc::new(runs),
                    extra: Arc::new(done.extra),
                    rescanned: done.rescanned,
                    errors,
                    took: job.started.elapsed(),
                });
                (self.table, self.note) = (None, String::new());
            }
            Ok(None) => self.note = "Cancelled.".into(),
            Err(mpsc::TryRecvError::Empty) => {
                // Keeps the progress bar moving.
                ctx.request_repaint_after(Duration::from_millis(200));
                return;
            }
            Err(mpsc::TryRecvError::Disconnected) => self.note = "Batch crashed.".into(),
        }
        self.job = None;
    }

    fn start(
        &mut self,
        ctx: &egui::Context,
        lib: &Arc<Library>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        let (done, cancel) = (
            Arc::new(AtomicUsize::new(0)),
            Arc::new(AtomicBool::new(false)),
        );
        let s = ScanSettings {
            // Like the single-run search: only a substance with an RI can be the reference.
            shift: cal
                .and(reference)
                .and_then(|(name, cutoff)| ShiftReference::find(&lib.compounds, name, cutoff)),
            untargeted: self.options.untargeted,
            ..self.applied.scan(cal.cloned())
        };
        let rescan = self.options.rescan;
        let quant = self.settings.quant.clone();
        let (paths, l, d, c, parallel) = (
            self.paths.clone(),
            lib.clone(),
            done.clone(),
            cancel.clone(),
            self.parallel,
        );
        let rx = spawn(ctx, move || {
            let (mut runs, errors) = scan_runs(&paths, &l.prepared, &s, parallel, &d, &c)?;
            let more = complete_scan(&mut runs, &l.compounds, &s, &quant, rescan, parallel, &c)?;
            Some((runs, errors, more))
        });
        self.job = Some(Job {
            total: self.paths.len(),
            done,
            cancel,
            started: Instant::now(),
            lib: lib.clone(),
            rx,
        });
    }

    /// Scan options, quant ions and the filters of the table.
    fn settings_grid(&mut self, ui: &mut egui::Ui) {
        let s = &mut self.settings;
        egui::Grid::new("batch settings").show(ui, |ui| {
            row(
                ui,
                "Runs in memory at once",
                "",
                "How many files are read at the same time; each run is analysed in parallel anyway. \
                 Fewer need less memory.",
                egui::DragValue::new(&mut self.parallel).range(1..=64),
            );
            row(
                ui,
                "Quant ions per substance",
                "QUANT_ION_NUMBER",
                "How many ions of each substance are integrated and summed. Only for substances whose library entry names no quant ions.",
                egui::DragValue::new(&mut s.quant.ion_number).range(1..=20),
            );
            row(
                ui,
                "Min. distance of quant ions",
                "QUANT_MIN_DIST_IONS",
                "Quant ions of one substance are preferably at least this many m/z apart. Only for substances whose library entry names no quant ions.",
                egui::DragValue::new(&mut s.quant.min_dist)
                    .range(0.0..=100.0)
                    .speed(0.1),
            );
            row(
                ui,
                "Quant ion quality cut-off",
                "QUANT_QUALITY_CUTOFF",
                "Ions whose peak quality is below this (lower is better) are picked first. Only for substances whose library entry names no quant ions.",
                egui::DragValue::new(&mut s.quant.quality_cutoff)
                    .range(0.0..=100.0)
                    .speed(0.1),
            );
            row(
                ui,
                "Replicate groups (regex on file name)",
                "",
                "Group = first capture group, else the whole match. A file without a match is a \
                 group of its own, so it passes any reproducibility filter. Empty: one group.",
                egui::TextEdit::singleline(&mut s.replicates).hint_text(r"e.g. ^(.+)_\d+\.cdf$"),
            );
            row(
                ui,
                "Min. reproducibility in a group",
                "",
                "Keep a substance only if some replicate group quantifies it in at least this share \
                 of its runs.",
                egui::Slider::new(&mut s.min_reproducibility, 0.0..=1.0),
            );
            row(
                ui,
                "Min. quant ions found",
                "",
                "A value counts only if at least this many of the substance's quant ions were found.",
                egui::DragValue::new(&mut s.min_ions).range(1..=20),
            );
            row(
                ui,
                "Min. signal-to-noise of an ion",
                "",
                "An ion counts only if its peak has at least this signal-to-noise ratio (old \
                 \"Req. S/N\").",
                egui::DragValue::new(&mut s.min_sn)
                    .range(0.0..=1e6)
                    .speed(0.1),
            );
            let hint = "A value counts only if the mean peak quality of its ions is at most this. \
                        Off: no limit.";
            let name = ui.label("Max. discrepancy index").on_hover_text(hint);
            let mut limit = s.max_discrepancy.is_finite();
            ui.horizontal(|ui| {
                ui.add(egui::Checkbox::without_text(&mut limit))
                    .labelled_by(name.id)
                    .on_hover_text(hint);
                let mut v = if limit { s.max_discrepancy } else { 0.0 };
                ui.add_enabled(limit, egui::DragValue::new(&mut v).range(0.0..=1e6).speed(0.1))
                    .labelled_by(name.id)
                    .on_hover_text(hint)
                    .on_disabled_hover_text("Tick the box to set a limit.");
                s.max_discrepancy = if limit { v } else { f32::INFINITY };
            });
            ui.end_row();
            row(
                ui,
                "Untargeted",
                "",
                "Also list components without a library hit: they are pooled over all runs into \
                 unknown substances. Applies to the next scan.",
                egui::Checkbox::without_text(&mut self.options.untargeted),
            );
            row(
                ui,
                "Rescan missing ions",
                "",
                "Search the ion traces again in the window of the RI difference around the expected \
                 positions. Applies to the next scan, with the quant ions chosen now.",
                egui::Checkbox::without_text(&mut self.options.rescan),
            );
        });
        crate::settings::ranges(
            ui,
            "Ions avoided as quant ions",
            "Ion ranges (m/z from–to) not picked as quant ions while other ions are left. Only for substances whose library entry names no quant ions.\nOld key: QUANT_EXCL_IONS",
            &mut s.quant.excluded,
        );
    }

    /// Whether `lib` was prepared for the search settings and the RT-only mode that apply now.
    fn library_ready(&self, lib: &Library) -> bool {
        lib.settings == self.applied.search() && lib.prepared.is_rt_only() == self.rt_only
    }

    /// Progress and Cancel while a scan runs, else Start and what it will use.
    fn job_row(
        &mut self,
        ui: &mut egui::Ui,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        ui.horizontal(|ui| match &self.job {
            Some(job) => {
                let d = job.done.load(Relaxed);
                let text = match d >= job.total {
                    // Pooling and the rescan run after the last file.
                    true => "pooling and rescanning…".into(),
                    false => format!("{d} / {} runs", job.total),
                };
                ui.add(
                    egui::ProgressBar::new(d as f32 / job.total.max(1) as f32)
                        .text(text)
                        .desired_width(300.0),
                );
                if ui.button("Cancel").clicked() {
                    job.cancel.store(true, Relaxed);
                }
            }
            None => {
                // Not while the library is rebuilt for new search settings or a switch of RT-only: it
                // would be scaled differently, and the calibration is hidden meanwhile.
                let ready = lib.filter(|l| !self.paths.is_empty() && self.library_ready(l));
                if ui
                    .add_enabled(ready.is_some(), egui::Button::new("Start"))
                    .clicked()
                    && let Some(l) = ready
                {
                    self.start(ui.ctx(), l, cal, reference);
                }
                match lib {
                    Some(l) => ui.label(format!(
                        "Library {}, {}",
                        l.name,
                        if cal.is_some() {
                            "RI calibrated"
                        } else {
                            "no RI calibration"
                        }
                    )),
                    None => ui.label("Open a library first (File → Open library)."),
                };
                if lib.is_some_and(|l| !self.library_ready(l)) {
                    ui.label("Library is being prepared for the new search or RT-only settings.");
                }
                ui.label(&self.note);
            }
        });
    }

    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        ui.horizontal(|ui| {
            if ui.button("Add runs…").clicked()
                && let Some(p) = rfd::FileDialog::new()
                    .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                    .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
                    .add_filter("FastFlight raw data (MsData.FFT)", &["FFT", "fft"])
                    .add_filter("Saved components (.cmp)", &["cmp", "CMP"])
                    .pick_files()
            {
                self.paths.extend(p);
                self.paths.sort();
                self.paths.dedup();
            }
            if ui.button("Clear").clicked() {
                self.paths.clear();
            }
            ui.label(format!("{} runs", self.paths.len()));
        });
        self.settings_grid(ui);
        self.job_row(ui, lib, cal, reference);
        let Some(scan) = &self.scan else { return };
        if self.table.as_ref().is_none_or(|t| t.0 != self.settings) {
            let t = table(&scan.runs, &scan.lib.compounds, &scan.extra, &self.settings)
                .map_err(|e| e.to_string());
            self.table = Some((self.settings.clone(), t));
            (self.shown, self.times, self.norm.edits) = (None, None, Edits::default());
        }
        ui.separator();
        ui.label(format!(
            "{} runs scanned in {} against {}; {} unknown substances, {} ions added by the rescan",
            scan.runs.len(),
            crate::view::took(scan.took),
            scan.lib.name,
            scan.extra.len(),
            scan.rescanned
        ));
        for e in &scan.errors {
            ui.colored_label(ui.visuals().error_fg_color, e);
        }
        match self.table.as_ref().map(|t| &t.1) {
            Some(Ok(t)) => {
                norm_controls(ui, &mut self.norm, t);
                pool_notes(ui, &self.norm, t, &scan.runs, &mut self.times);
                if self.shown.as_ref().is_none_or(|s| s.0 != self.norm) {
                    let times = self.times.as_ref().map(|t| t.0.as_slice());
                    self.shown = Some((self.norm.clone(), shown_table(t, &self.norm, times)));
                    self.generation += 1;
                }
                if let Some((_, shown)) = &self.shown {
                    ui.horizontal(|ui| {
                        let ex = Export {
                            scan,
                            settings: &self.settings,
                            extra: shown_settings(&self.norm, shown, t.runs.len()),
                            ions: self.view.export_ions,
                        };
                        if let Some(note) = export_button(ui, shown, &ex) {
                            self.note = note;
                        }
                        if self.frag_job.is_none()
                            && let Some(rx) = frag_button(ui, shown, scan)
                        {
                            (self.note, self.frag_job) = ("Exporting MIDs…".into(), Some(rx));
                        }
                    });
                    view_controls(
                        ui,
                        &mut self.view,
                        &mut self.norm.edits,
                        &group_names(shown),
                        self.norm.z_score,
                    );
                    if self.view.totals && self.totals.0 != self.generation {
                        self.totals = (self.generation, shown.total_intensities());
                    }
                    let totals = self.view.totals.then_some(self.totals.1.as_slice());
                    let d = Detail {
                        scan,
                        settings: &self.settings,
                    };
                    if let Some(act) = table_view(ui, shown, &self.norm, &self.view, totals, &d) {
                        self.norm.edits.apply(act, t.runs.len());
                    }
                }
            }
            Some(Err(e)) => {
                ui.colored_label(
                    ui.visuals().error_fg_color,
                    format!("Replicate pattern: {e}"),
                );
            }
            None => {}
        }
    }
}

/// "Export MIDs" of the old batch table: isotopomer shares of the fragments in each shown
/// substance's `FRAG` feature, written on a worker thread that sends back a note.
fn frag_button(ui: &mut egui::Ui, t: &Table, scan: &Scan) -> Option<mpsc::Receiver<String>> {
    if !ui
        .button("Export MIDs (FRAG)…")
        .on_hover_text("Total ion counts and isotopomer shares of the fragments in each substance's FRAG feature, e.g. (174/8)(261/5).")
        .clicked()
    {
        return None;
    }
    let p = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name("mids.csv")
        .save_file()?;
    let rows: Vec<usize> = t.rows.iter().map(|r| r.compound).collect();
    let sources = t.sources.clone();
    let (lib, extra, runs) = (scan.lib.clone(), scan.extra.clone(), scan.runs.clone());
    Some(spawn(ui.ctx(), move || {
        let n = lib.compounds.len();
        // Unknowns of an untargeted scan follow the library's entries.
        let compounds: Vec<(usize, &LibraryCompound)> = (rows.into_iter())
            .filter_map(|i| {
                Some((
                    i,
                    lib.compounds
                        .get(i)
                        .or_else(|| extra.get(i.checked_sub(n)?))?,
                ))
            })
            .collect();
        let e = frag_export(&shown_runs(&runs, &sources), &compounds);
        let mut note = match std::fs::write(&p, &e.csv) {
            Ok(()) => format!("Saved {}", p.display()),
            Err(err) => format!("Export failed: {err}"),
        };
        for m in &e.errors {
            note += &format!("\n{m}");
        }
        match e.warnings.len() {
            0 => {}
            1 => note += &format!("\n{}", e.warnings[0]),
            n => note += &format!("\n{n} substances have no fragment (FRAG) defined."),
        }
        note
    }))
}

/// The scanned runs of the shown table's columns, in its order: removed runs stay out of exports.
fn shown_runs(runs: &[RunResult], sources: &[usize]) -> Vec<RunResult> {
    sources
        .iter()
        .filter_map(|i| runs.get(*i).cloned())
        .collect()
}

/// What an export writes besides the table itself.
struct Export<'a> {
    scan: &'a Scan,
    settings: &'a TableSettings,
    /// Settings of the shown table that `TableSettings` does not hold (normalisation, removals).
    extra: Vec<(&'static str, String)>,
    /// Also write the integral of every quant ion.
    ions: bool,
}

/// Normalisation and edits of the shown table, for the settings file next to an export.
fn shown_settings(n: &Norm, t: &Table, base_runs: usize) -> Vec<(&'static str, String)> {
    let name = |c: Option<usize>| {
        (t.rows.iter().find(|r| Some(r.compound) == c))
            .map_or_else(|| "none chosen".into(), |r| r.name.clone())
    };
    let normalisation = match n.kind {
        Kind::None => "none".into(),
        Kind::InternalStandard => format!("internal standard {}", name(n.standard)),
        Kind::TotalSignal => "total signal".into(),
        Kind::Pool => format!(
            "pool group {}, nearest {} runs, {}",
            n.pool_group,
            n.pool_nearest,
            if n.pool_median { "median" } else { "mean" }
        ),
    };
    vec![
        ("NORMALISATION", normalisation),
        ("Z_SCORE", n.z_score.to_string()),
        (
            "REMOVED_RUNS",
            (base_runs - t.runs.len().min(base_runs)).to_string(),
        ),
        ("REMOVED_SUBSTANCES", n.edits.rows.len().to_string()),
    ]
}

/// The export has the values the table shows, plus the replicate groups and the settings used
/// next to it, and with `ex.ions` the raw integral of every quant ion. Returns a note after a
/// save attempt.
fn export_button(ui: &mut egui::Ui, t: &Table, ex: &Export) -> Option<String> {
    if !ui.button("Export CSV…").clicked() {
        return None;
    }
    let p = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name("results.csv")
        .save_file()?;
    Some(match save_results(&p, t, ex) {
        Ok(note) => note,
        Err(e) => format!("Export failed: {e}"),
    })
}

/// Writes the results CSV `p` and the files next to it; the note names what was saved and what
/// could not be.
fn save_results(p: &Path, t: &Table, ex: &Export) -> std::io::Result<String> {
    let create = |p: &Path| std::fs::File::create(p).map(std::io::BufWriter::new);
    t.write_csv(create(p)?)?;
    let sp = settings_path(p);
    let mut results = vec![
        ("replicate groups", t.write_groups_file(p).transpose()),
        (
            "settings",
            Some(
                create(&sp)
                    .and_then(|f| ex.settings.write_csv(&ex.extra, f))
                    .map(|()| sp),
            ),
        ),
    ];
    if !ex.ions {
        // An older ions file would not match this results file.
        std::fs::remove_file(ions_path(p)).ok();
    } else {
        let ip = ions_path(p);
        let runs: Vec<&RunResult> = (t.sources.iter())
            .filter_map(|i| ex.scan.runs.get(*i))
            .collect();
        let r = create(&ip).and_then(|f| t.write_ions_csv(&runs, ex.settings, f));
        results.push(("quant ions", Some(r.map(|()| ip))));
    }
    let mut saved = vec![p.display().to_string()];
    let mut failed = Vec::new();
    for (what, r) in results {
        match r {
            Some(Ok(path)) => saved.push(path.display().to_string()),
            Some(Err(e)) => failed.push(format!("{what}: {e}")),
            None => {}
        }
    }
    let mut text = format!("Saved {}", saved.join(", "));
    if !failed.is_empty() {
        text += &format!("; not saved: {}", failed.join("; "));
    }
    Ok(text)
}

/// The table as the user's edits and the normalisation leave it: runs are removed and moved
/// first, substances are removed last.
fn shown_table(t: &Table, n: &Norm, times: Option<&[i64]>) -> Table {
    let cols = n.edits.cols.as_deref();
    let edited = cols.map(|c| t.select_runs(c));
    let base = edited.as_ref().unwrap_or(t);
    // The file dates follow the runs that are left.
    let times: Vec<i64> = match (times, cols) {
        (Some(ts), Some(cols)) => cols
            .iter()
            .map(|c| ts.get(*c).copied().unwrap_or(0))
            .collect(),
        (Some(ts), None) => ts.to_vec(),
        (None, _) => Vec::new(),
    };
    apply(base, &mode(n, base, &times), n.z_score).without_rows(&n.edits.rows)
}

/// With `z_scored` there are no group statistics: the old program computes them on the normalised
/// values only, and a mean, hits above 0 or a missing value as 0 say nothing about z-scores.
fn view_controls(
    ui: &mut egui::Ui,
    v: &mut View,
    edits: &mut Edits,
    groups: &[&str],
    z_scored: bool,
) {
    if z_scored && v.mode == Mode::Groups {
        v.mode = Mode::Values;
    }
    // Two rows: the window must not need more than the room beside the run list.
    ui.horizontal(|ui| {
        for (mode, title) in [
            (Mode::Values, "Values"),
            (Mode::Times, "Retention times (min)"),
            (Mode::Groups, "Group statistics"),
        ] {
            let locked = mode == Mode::Groups && z_scored;
            let r = ui
                .add_enabled_ui(!locked, |ui| ui.selectable_value(&mut v.mode, mode, title))
                .inner;
            if mode == Mode::Groups {
                r.on_hover_text("Per replicate group: mean, relative standard error (standard error / mean), hits (#) and the p-value of Student's t-test against the blank group.")
                    .on_disabled_hover_text("Not with z-scores: a missing value cannot count as 0 there, and \"#\" (values above 0) means nothing. Turn the z-score off; the group statistics then use the normalised values, like the old program.");
            }
        }
    });
    ui.horizontal(|ui| {
        ui.add_enabled(
            v.mode == Mode::Values,
            egui::Checkbox::new(&mut v.totals, "Total sample intensities"),
        )
        .on_hover_text(
            "Adds a first row with the sum of each run's values (a missing value counts as 0).",
        )
        .on_disabled_hover_text("Only in the Values view.");
        ui.add_enabled(
            v.mode == Mode::Values,
            egui::Checkbox::new(&mut v.ion_rows, "Quant ion lines"),
        )
        .on_hover_text("Adds a line under each substance for every quant ion, with its raw integral per run.")
        .on_disabled_hover_text("Only in the Values view.");
        ui.checkbox(&mut v.export_ions, "Export quant ions")
            .on_hover_text("Also writes results.ions.csv with the raw integral of every quant ion.");
        if v.mode == Mode::Groups {
            let blank = ui.label("Blank group");
            ui.add(egui::TextEdit::singleline(&mut v.blanks).desired_width(80.0))
                .labelled_by(blank.id)
                .on_hover_text(
                    "Replicate group of the blank runs; every other group gets a p-value against it.",
                );
            if !groups.contains(&v.blanks.as_str()) {
                ui.weak("no such group, so no p-values");
            }
        }
        if !edits.is_empty()
            && ui
                .button("Restore removed")
                .on_hover_text(
                    "Brings back the removed runs and substances and the original order of the runs.",
                )
                .clicked()
        {
            edits.restore();
        }
    });
}

/// Replicate groups in name order, like the old statistics table.
fn group_names(t: &Table) -> Vec<&str> {
    let mut g: Vec<&str> = t.groups.iter().map(String::as_str).collect();
    g.sort_unstable();
    g.dedup();
    g
}

/// A run without replicate group (empty pattern) belongs to the whole batch.
fn group_label(g: &str) -> String {
    if g.is_empty() {
        "all runs".into()
    } else {
        g.into()
    }
}

/// A statistic that is undefined (NaN) shows as a dash.
fn number(v: f64, decimals: usize) -> String {
    if v.is_finite() {
        format!("{v:.decimals$}")
    } else {
        "–".into()
    }
}

/// A mean with 6 significant digits, like the old table (Qt shows a double as `'g'` with 6), so
/// normalised fractions stay readable; without exponent.
fn significant(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return number(v, 0);
    }
    if v.abs() < 1e-10 {
        return format!("{v:.5e}");
    }
    let decimals = (5 - v.abs().log10().floor() as i32).clamp(0, 15) as usize;
    let s = format!("{v:.decimals$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').into()
    } else {
        s
    }
}

/// Titles and sub-titles of the columns after "Quant ions".
fn heads(t: &Table, v: &View, groups: &[&str], blank: bool) -> Vec<(String, String)> {
    if v.mode != Mode::Groups {
        return (t.runs.iter().zip(&t.groups))
            .map(|(r, g)| (r.clone(), group_label(g)))
            .collect();
    }
    let stats = ["Mean", "Rel. Std. Err.", "#", "P (blank)"];
    let per_group = if blank { 4 } else { 3 };
    (groups.iter())
        .flat_map(|g| {
            stats[..per_group]
                .iter()
                .map(move |s| (s.to_string(), group_label(g)))
        })
        .collect()
}

/// The cells of one substance in "Group statistics": per group mean, relative standard error,
/// hits and (with a blank group) the p-value against the blanks.
fn group_cells(t: &Table, row: usize, groups: &[&str], blank: Option<&str>) -> Vec<String> {
    let values = |g: &str| group_values(t, row, &[g.to_string()], true);
    let blanks = blank.map(values);
    let mut out = Vec::new();
    for g in groups {
        let x = values(g);
        let s = summary(&x);
        out.push(significant(s.mean));
        out.push(number(s.rel_std_err, 3));
        out.push(s.hits.to_string());
        if let Some(b) = &blanks {
            let p = (blank != Some(g)).then(|| student(&x, b)).flatten();
            out.push(p.map_or("–".into(), |p| format!("{p:.4}")));
        }
    }
    out
}

/// Cells of substance `row` after "Quant ions".
fn cells(t: &Table, v: &View, groups: &[&str], blank: bool, row: usize) -> Vec<String> {
    let r = &t.rows[row];
    match v.mode {
        Mode::Values => (r.values.iter())
            .map(|v| v.map_or("–".into(), |v| format!("{v:.0}")))
            .collect(),
        Mode::Times => (r.rts_s.iter())
            .map(|rt| rt.map_or("–".into(), |s| format!("{:.2}", s / 60.0)))
            .collect(),
        Mode::Groups => group_cells(t, row, groups, blank.then_some(v.blanks.as_str())),
    }
}

/// Buttons in the header of run `j` of `n`: move left, move right, remove.
fn run_buttons(ui: &mut egui::Ui, j: usize, n: usize) -> Option<Act> {
    let mut act = None;
    ui.horizontal(|ui| {
        if j > 0
            && ui
                .small_button("‹")
                .on_hover_text("Move run left")
                .clicked()
        {
            act = Some(Act::MoveRun(j, j - 1));
        }
        if j + 1 < n
            && ui
                .small_button("›")
                .on_hover_text("Move run right")
                .clicked()
        {
            act = Some(Act::MoveRun(j, j + 1));
        }
        if ui
            .small_button("×")
            .on_hover_text("Remove this run from the results (Restore removed brings it back)")
            .clicked()
        {
            act = Some(Act::DropRun(j));
        }
    });
    act
}

/// Where the cells of the shown table come from, for the details behind a cell.
struct Detail<'a> {
    scan: &'a Scan,
    settings: &'a TableSettings,
}

impl Detail<'_> {
    /// The hit of `compound` in the run behind column `col` of `t`.
    fn found(&self, t: &Table, col: usize, compound: usize) -> Option<&Found> {
        let run = self.scan.runs.get(*t.sources.get(col)?)?;
        run.found.iter().find(|f| f.compound == compound)
    }

    /// The library hits of the component behind a cell: the assigned one first (red), then the
    /// alternatives, like the old tooltip.
    fn hits_ui(&self, ui: &mut egui::Ui, t: &Table, col: usize, compound: usize) {
        let Some(f) = self.found(t, col, compound) else {
            ui.label("Not found in this run.");
            return;
        };
        let lib = &self.scan.lib.compounds;
        for (k, (name, score)) in hit_list(f, lib, &self.scan.extra).into_iter().enumerate() {
            let text = format!("{name} ({score:.3})");
            match k {
                0 => ui.colored_label(ui.visuals().error_fg_color, text),
                _ => ui.label(text),
            };
        }
    }

    /// The raw integral of quant ion `k` of substance `row` in column `col`, where the cell counts
    /// as quantified (whatever the normalisation made of its value).
    fn ion(&self, t: &Table, row: usize, col: usize, k: usize) -> Option<f64> {
        let r = &t.rows[row];
        let f = self.found(t, col, r.compound)?;
        cell_ions(f, &r.ions, self.settings)[k]
    }
}

/// One line of the table body.
#[derive(Clone, Copy)]
enum Line {
    Total,
    Substance(usize),
    /// Quant ion `k` of substance `row`.
    Ion(usize, usize),
}

fn body_lines(t: &Table, v: &View, with_total: bool) -> Vec<Line> {
    let ions = v.ion_rows && v.mode == Mode::Values;
    let mut lines = Vec::with_capacity(t.rows.len() + 1);
    lines.extend(with_total.then_some(Line::Total));
    for (i, r) in t.rows.iter().enumerate() {
        lines.push(Line::Substance(i));
        if ions {
            lines.extend((0..r.ions.len()).map(|k| Line::Ion(i, k)));
        }
    }
    lines
}

/// What the table body needs to draw a line.
struct Body<'a> {
    t: &'a Table,
    v: &'a View,
    d: &'a Detail<'a>,
    groups: &'a [&'a str],
    blank: bool,
    totals: &'a [f64],
}

impl Body<'_> {
    /// Draws one line; returns what the user clicked in it.
    fn draw(&self, row: &mut egui_extras::TableRow, line: Line) -> Option<Act> {
        match line {
            Line::Total => total_row(row, self.totals),
            Line::Substance(i) => return self.substance(row, i),
            Line::Ion(i, k) => self.ion(row, i, k),
        }
        None
    }

    fn substance(&self, row: &mut egui_extras::TableRow, i: usize) -> Option<Act> {
        let (t, r) = (self.t, &self.t.rows[i]);
        let mut act = None;
        row.col(|ui| {
            if ui
                .small_button("×")
                .on_hover_text(
                    "Remove this substance from the results (Restore removed brings it back)",
                )
                .clicked()
            {
                act = Some(Act::DropRow(r.compound));
            }
        });
        let ions: Vec<String> = r.ions.iter().map(u32::to_string).collect();
        let lead = [
            r.name.clone(),
            format!("{:.2}", r.rt_s / 60.0),
            ions.join(" "),
        ];
        for text in lead {
            row.col(|ui| {
                ui.label(text);
            });
        }
        for (col, text) in cells(t, self.v, self.groups, self.blank, i)
            .into_iter()
            .enumerate()
        {
            row.col(|ui| {
                let label = ui.label(text);
                if self.v.mode == Mode::Values {
                    label.on_hover_ui(|ui| self.d.hits_ui(ui, t, col, r.compound));
                }
            });
        }
        act
    }

    fn ion(&self, row: &mut egui_extras::TableRow, i: usize, k: usize) {
        let mass = self.t.rows[i].ions[k];
        row.col(|_| {});
        row.col(|ui| {
            ui.weak(format!("    ion {mass}"));
        });
        row.col(|_| {});
        row.col(|_| {});
        for col in 0..self.t.runs.len() {
            row.col(|ui| {
                let v = self.d.ion(self.t, i, col, k);
                ui.weak(v.map_or("–".into(), |v| format!("{v:.0}")));
            });
        }
    }
}

/// Rows are virtualised, so only the visible ones are laid out each frame. Returns what the user
/// clicked to change the table.
/// shortcut: all columns are laid out per visible row; virtualise columns too if studies with
/// thousands of runs feel slow. The total row scrolls with the rows instead of sticking.
fn table_view(
    ui: &mut egui::Ui,
    t: &Table,
    norm: &Norm,
    v: &View,
    totals: Option<&[f64]>,
    d: &Detail,
) -> Option<Act> {
    let what = match (v.mode, norm.kind, norm.z_score) {
        (Mode::Times, ..) => "retention times of the hits in minutes",
        (Mode::Groups, ..) => "per replicate group; a missing value counts as 0",
        (_, Kind::None, false) => "summed quant ion integrals",
        (_, _, false) => "normalised quant ion integrals",
        (_, Kind::None, true) => "z-scores of the summed integrals",
        (_, _, true) => "z-scores of the normalised integrals",
    };
    ui.label(format!(
        "{} substances × {} runs; values are {what}",
        t.rows.len(),
        t.runs.len()
    ));
    let groups = group_names(t);
    let blank = v.mode == Mode::Groups && groups.contains(&v.blanks.as_str());
    let heads = heads(t, v, &groups, blank);
    let totals = totals.filter(|_| v.mode == Mode::Values);
    let lines = body_lines(t, v, totals.is_some());
    let body = Body {
        t,
        v,
        d,
        groups: &groups,
        blank,
        totals: totals.unwrap_or_default(),
    };
    let mut act = None;
    egui::ScrollArea::horizontal().show(ui, |ui| {
        // Rows have a fixed height, so cells must not wrap.
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        TableBuilder::new(ui)
            .id_salt("batch table")
            .striped(true)
            .column(Column::exact(28.0))
            .column(Column::auto().at_least(160.0))
            .columns(Column::auto(), 2)
            .columns(Column::auto().at_least(70.0), heads.len())
            .header(54.0, |mut h| {
                h.col(|_| {});
                for title in ["Substance", "RT (min)", "Quant ions"] {
                    h.col(|ui| {
                        ui.strong(title);
                    });
                }
                for (j, (title, sub)) in heads.iter().enumerate() {
                    h.col(|ui| {
                        ui.vertical(|ui| {
                            ui.strong(title);
                            ui.label(sub);
                            if v.mode != Mode::Groups {
                                act = run_buttons(ui, j, heads.len()).or(act);
                            }
                        });
                    });
                }
            })
            .body(|rows| {
                rows.rows(18.0, lines.len(), |mut row| {
                    let line = lines[row.index()];
                    act = body.draw(&mut row, line).or(act);
                });
            });
    });
    act
}

fn total_row(row: &mut egui_extras::TableRow, totals: &[f64]) {
    row.col(|_| {});
    row.col(|ui| {
        ui.strong("Total sample intensity");
    });
    row.col(|_| {});
    row.col(|_| {});
    for t in totals {
        row.col(|ui| {
            ui.strong(format!("{t:.0}"));
        });
    }
}

fn norm_controls(ui: &mut egui::Ui, n: &mut Norm, t: &Table) {
    ui.horizontal(|ui| {
        ui.label("Normalise by");
        for (kind, title) in [
            (Kind::None, "nothing"),
            (Kind::InternalStandard, "internal standard"),
            (Kind::TotalSignal, "total signal"),
            (Kind::Pool, "pool"),
        ] {
            ui.selectable_value(&mut n.kind, kind, title);
        }
        ui.checkbox(&mut n.z_score, "z-score");
    });
    match n.kind {
        Kind::InternalStandard => {
            ui.horizontal(|ui| {
                let label = ui.label("Internal standard");
                let name = |c: Option<usize>| {
                    t.rows
                        .iter()
                        .find(|r| Some(r.compound) == c)
                        .map_or("choose…", |r| r.name.as_str())
                };
                egui::ComboBox::from_id_salt("norm standard")
                    .selected_text(name(n.standard))
                    .show_ui(ui, |ui| {
                        for r in &t.rows {
                            ui.selectable_value(&mut n.standard, Some(r.compound), &r.name);
                        }
                    })
                    .response
                    .labelled_by(label.id);
            });
        }
        Kind::Pool => {
            ui.horizontal(|ui| {
                let group = ui.label("Pool runs are the replicate group");
                ui.add(egui::TextEdit::singleline(&mut n.pool_group).desired_width(80.0))
                    .labelled_by(group.id);
                let nearest = ui.add(egui::DragValue::new(&mut n.pool_nearest).range(1..=99));
                let label = ui.label("nearest pool runs");
                nearest.labelled_by(label.id);
                ui.checkbox(&mut n.pool_median, "median");
            })
            .response
            .on_hover_text("Each run is divided by the mean or median of the nearest pool runs in time (file date). Pool runs keep their raw values.");
        }
        Kind::None | Kind::TotalSignal => {}
    }
}

fn mode(n: &Norm, t: &Table, times: &[i64]) -> Normalization {
    match n.kind {
        Kind::None => Normalization::None,
        Kind::TotalSignal => Normalization::TotalSignal,
        Kind::InternalStandard => Normalization::InternalStandard(
            t.rows
                .iter()
                .position(|r| Some(r.compound) == n.standard)
                .unwrap_or(usize::MAX),
        ),
        Kind::Pool => Normalization::Pool(Pool {
            group: n.pool_group.clone(),
            times: times.to_vec(),
            nearest: n.pool_nearest,
            median: n.pool_median,
        }),
    }
}

/// The old program orders runs by acquisition date; the file date stands in for it.
fn run_times(t: &Table, runs: &[RunResult]) -> (Vec<i64>, usize) {
    let mut missing = 0;
    let times = t
        .sources
        .iter()
        .map(|i| {
            runs.get(*i)
                .and_then(|r| std::fs::metadata(&r.path).ok()?.modified().ok())
                .and_then(|m| m.duration_since(std::time::UNIX_EPOCH).ok())
                .map_or_else(
                    || {
                        missing += 1;
                        0
                    },
                    |d| d.as_secs() as i64,
                )
        })
        .collect();
    (times, missing)
}

/// Reads the run dates once and warns when the pool settings cannot work.
fn pool_notes(
    ui: &mut egui::Ui,
    n: &Norm,
    t: &Table,
    runs: &[RunResult],
    times: &mut Option<(Vec<i64>, usize)>,
) {
    if n.kind != Kind::Pool {
        return;
    }
    let (_, missing) = times.get_or_insert_with(|| run_times(t, runs));
    if *missing > 0 {
        ui.label(format!(
            "No file date for {missing} runs; they count as time 0."
        ));
    }
    if !t.groups.contains(&n.pool_group) {
        ui.label("No run is in the pool group: every normalised value is empty.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quant_settings_follow_the_applied_ones_but_keep_local_edits() {
        let mut b = Batch::default();
        let mut applied = GcmsSettings {
            quant_ion_number: 5,
            ..Default::default()
        };
        b.use_settings(&applied, false);
        assert_eq!(b.settings.quant.ion_number, 5);
        // An edit in the window survives settings that leave QUANT_* alone ...
        b.settings.quant.ion_number = 2;
        applied.ls_ri_diff += 1.0;
        b.use_settings(&applied, false);
        assert_eq!(b.settings.quant.ion_number, 2);
        // ... and the next change of QUANT_* replaces it.
        applied.quant_ion_number = 4;
        b.use_settings(&applied, false);
        assert_eq!(b.settings.quant.ion_number, 4);
    }

    #[test]
    fn edits_drop_move_and_restore_columns_and_rows() {
        let mut e = Edits::default();
        e.apply(Act::MoveRun(0, 1), 3);
        assert_eq!(e.cols, Some(vec![1, 0, 2]));
        e.apply(Act::DropRun(1), 3);
        assert_eq!(e.cols, Some(vec![1, 2]));
        // Out of range clicks (a stale frame) change nothing.
        e.apply(Act::DropRun(5), 3);
        e.apply(Act::MoveRun(1, 2), 3);
        assert_eq!(e.cols, Some(vec![1, 2]));
        e.apply(Act::DropRow(7), 3);
        assert_eq!(e.rows, [7]);
        assert!(!e.is_empty());
    }

    #[test]
    fn group_statistics_are_readable_after_normalising_and_locked_for_z_scores() {
        let t = Table {
            runs: vec!["a".into(), "b".into(), "c".into()],
            groups: vec!["A".into(); 3],
            rows: vec![md_core::batch::Row {
                compound: 0,
                name: "s0".into(),
                rt_s: 60.0,
                ions: vec![73],
                values: vec![Some(1.0), Some(2.0), Some(10.0)],
                rts_s: vec![Some(60.0); 3],
            }],
            ..Default::default()
        };
        assert_eq!(
            group_cells(&t, 0, &["A"], None)[..3],
            ["4.33333", "0.657", "3"]
        );
        // By total signal every value is 1: the old code rounded the mean to "1" and fractions to 0.
        let mut n = Norm {
            kind: Kind::TotalSignal,
            ..Default::default()
        };
        let shown = shown_table(&t, &n, None);
        assert_eq!(group_cells(&shown, 0, &["A"], None)[0], "1");
        let shown = Table {
            rows: vec![md_core::batch::Row {
                values: vec![Some(0.0123), Some(0.0124), Some(0.0128)],
                ..t.rows[0].clone()
            }],
            ..t.clone()
        };
        assert_eq!(group_cells(&shown, 0, &["A"], None)[0], "0.0125");
        assert_eq!(significant(1_234_567.0), "1234567");
        assert_eq!(significant(-0.0123456789), "-0.0123457");
        assert_eq!(significant(f64::NAN), "–");
        assert_eq!(significant(-1e-20), "-1.00000e-20");

        // z-scores of 1, 2, 10: only one is above 0, so "#" would say 1 of 3 hits.
        n = Norm {
            z_score: true,
            ..Default::default()
        };
        let z: Vec<String> = (shown_table(&t, &n, None).rows[0].values.iter())
            .map(|v| format!("{:.2}", v.unwrap()))
            .collect();
        assert_eq!(z, ["-0.68", "-0.47", "1.15"]);
        let mut v = View {
            mode: Mode::Groups,
            ..Default::default()
        };
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(Default::default(), |ui| {
            view_controls(ui, &mut v, &mut Edits::default(), &["A"], true);
        });
        out.textures_delta.clear();
        assert!(v.mode == Mode::Values);
    }

    #[test]
    fn the_shown_table_has_no_removed_runs_or_substances() {
        let row = |compound, values: [f64; 3]| md_core::batch::Row {
            compound,
            name: format!("s{compound}"),
            rt_s: 60.0,
            ions: vec![73],
            values: values.map(Some).to_vec(),
            rts_s: vec![Some(60.0); 3],
        };
        let t = Table {
            runs: vec!["a".into(), "b".into(), "c".into()],
            groups: vec![String::new(); 3],
            rows: vec![row(0, [1.0, 2.0, 3.0]), row(1, [4.0, 5.0, 6.0])],
            ..Default::default()
        };
        let mut n = Norm::default();
        n.edits.apply(Act::DropRun(0), 3);
        n.edits.apply(Act::DropRow(0), 3);
        let shown = shown_table(&t, &n, None);
        assert_eq!(shown.runs, ["b", "c"]);
        assert_eq!(shown.rows.len(), 1);
        assert_eq!(shown.rows[0].values, [Some(5.0), Some(6.0)]);
        // Normalised by total signal: the removed run no longer counts, the removed row still does.
        n.kind = Kind::TotalSignal;
        let shown = shown_table(&t, &n, None);
        assert_eq!(shown.rows[0].values, [Some(5.0 / 7.0), Some(6.0 / 9.0)]);
        // "Restore removed" brings back the runs, the substances and the original order.
        n.kind = Kind::None;
        n.edits.apply(Act::MoveRun(0, 1), 3);
        n.edits.restore();
        let shown = shown_table(&t, &n, None);
        assert_eq!(shown.runs, t.runs);
        assert_eq!(shown.rows, t.rows);
    }

    #[test]
    fn exports_use_the_runs_of_the_shown_table_in_its_order() {
        let run = |n: &str| RunResult {
            path: format!("dir/{n}.cdf").into(),
            ..Default::default()
        };
        let runs = [run("a"), run("b"), run("c")];
        let shown = shown_runs(&runs, &[2, 0]);
        let names: Vec<_> = shown
            .iter()
            .map(|r| r.path.file_name().unwrap().to_owned())
            .collect();
        assert_eq!(names, ["c.cdf", "a.cdf"]);
    }

    #[test]
    fn runs_with_the_same_name_in_two_folders_are_told_apart() {
        let root = std::env::temp_dir().join(format!("md-same-name-{}", std::process::id()));
        let at = |dir: &str, secs: u64| {
            std::fs::create_dir_all(root.join(dir)).unwrap();
            let path = root.join(dir).join("a.cdf");
            let file = std::fs::File::create(&path).unwrap();
            file.set_modified(std::time::UNIX_EPOCH + Duration::from_secs(secs))
                .unwrap();
            RunResult {
                path,
                ..Default::default()
            }
        };
        let runs = [at("d1", 100), at("d2", 200)];
        let t = md_core::batch::table(&runs, &[], &[], &Default::default()).unwrap();
        assert_eq!(run_times(&t, &runs), (vec![100, 200], 0));
        let t = t.select_runs(&[1, 0]);
        let shown: Vec<_> = shown_runs(&runs, &t.sources)
            .into_iter()
            .map(|r| r.path)
            .collect();
        assert_eq!(shown, [runs[1].path.clone(), runs[0].path.clone()]);
        std::fs::remove_dir_all(&root).unwrap();
    }

    #[test]
    fn start_waits_for_the_library_to_match_search_settings_and_rt_only() {
        let search = GcmsSettings::default().search();
        let lib = |rt_only: bool| Library {
            path: PathBuf::new(),
            name: String::new(),
            prepared: if rt_only {
                md_core::search::PreparedLibrary::new_rt_only(&[], &search)
            } else {
                md_core::search::PreparedLibrary::new(&[], &search)
            },
            compounds: Vec::new(),
            settings: search.clone(),
        };
        let mut b = Batch::default();
        assert!(b.library_ready(&lib(false)));
        assert!(!b.library_ready(&lib(true)));
        b.use_settings(&GcmsSettings::default(), true);
        assert!(b.library_ready(&lib(true)));
        assert!(!b.library_ready(&lib(false)));
    }

    /// One run with one hit of two quant ions and an alternative library entry.
    fn one_hit_scan() -> (Scan, Table) {
        use md_core::batch::{IonPeak, Row};
        let search = GcmsSettings::default().search();
        let compound = |n: &str| LibraryCompound::new(n);
        let ion = |mass, integral| IonPeak {
            mass,
            integral,
            quality: 0.5,
            signal_to_noise: 10.0,
        };
        let run = RunResult {
            path: "dir/a.cdf".into(),
            found: vec![Found {
                compound: 0,
                score: 0.9,
                rt_s: 60.0,
                ions: vec![ion(73, 500.0), ion(147, 100.0)],
                rescanned: false,
                rescued: Vec::new(),
                alternatives: vec![(1, 0.8)],
            }],
            ..Default::default()
        };
        let lib = Library {
            path: PathBuf::new(),
            name: String::new(),
            prepared: md_core::search::PreparedLibrary::new(&[], &search),
            compounds: vec![compound("Glucose"), compound("Fructose")],
            settings: search,
        };
        let scan = Scan {
            lib: Arc::new(lib),
            runs: Arc::new(vec![run]),
            extra: Arc::default(),
            rescanned: 0,
            errors: Vec::new(),
            took: Duration::ZERO,
        };
        let table = Table {
            runs: vec!["a.cdf".into()],
            sources: vec![0],
            groups: vec![String::new()],
            rows: vec![Row {
                compound: 0,
                name: "Glucose".into(),
                rt_s: 60.0,
                ions: vec![73, 147],
                values: vec![Some(600.0)],
                rts_s: vec![Some(60.0)],
            }],
        };
        (scan, table)
    }

    #[test]
    fn a_cell_shows_its_ion_integrals_and_the_other_hits_of_its_component() {
        let (scan, t) = one_hit_scan();
        let settings = TableSettings::default();
        let d = Detail {
            scan: &scan,
            settings: &settings,
        };
        assert_eq!(d.ion(&t, 0, 0, 0), Some(500.0));
        assert_eq!(d.ion(&t, 0, 0, 1), Some(100.0));
        let f = d.found(&t, 0, 0).unwrap();
        let hits = hit_list(f, &scan.lib.compounds, &scan.extra);
        assert_eq!(hits, [("Glucose".into(), 0.9), ("Fructose".into(), 0.8)]);
        // The ion lines follow their substance; without the checkbox there are none.
        let mut v = View::default();
        assert_eq!(body_lines(&t, &v, false).len(), 1);
        v.ion_rows = true;
        assert_eq!(body_lines(&t, &v, true).len(), 1 + 1 + 2);
        v.mode = Mode::Times;
        assert_eq!(body_lines(&t, &v, false).len(), 1);
    }

    #[test]
    fn an_export_writes_the_settings_and_on_request_the_ion_integrals() {
        let (scan, t) = one_hit_scan();
        let dir = std::env::temp_dir().join(format!("md-app-export-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let csv = dir.join("results.csv");
        let settings = TableSettings::default();
        let mut ex = Export {
            scan: &scan,
            settings: &settings,
            extra: vec![("NORMALISATION", "none".into())],
            ions: false,
        };
        let note = save_results(&csv, &t, &ex).unwrap();
        assert!(note.contains("results.settings.csv"), "{note}");
        assert!(!dir.join("results.ions.csv").exists());
        let text = std::fs::read_to_string(dir.join("results.settings.csv")).unwrap();
        assert!(text.ends_with("NORMALISATION,none\n"), "{text}");

        ex.ions = true;
        save_results(&csv, &t, &ex).unwrap();
        let ions = std::fs::read_to_string(dir.join("results.ions.csv")).unwrap();
        assert_eq!(
            ions,
            "\u{feff}Substance,Quant_ion,a.cdf\nGlucose,73,500\nGlucose,147,100\n"
        );
        std::fs::remove_dir_all(&dir).ok();
    }
}
