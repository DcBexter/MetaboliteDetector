//! Batch quantification window: scans many runs on a worker thread against the open library and
//! shows the substance × run table (`md_core::batch`). Filters and groups apply without a rescan.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use crate::layout::ToolWindow;
use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::batch::{
    Completed, RunResult, ScanSettings, ShiftReference, Table, TableSettings, complete_scan,
    scan_runs, table,
};
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
    /// Replicate group of the blank runs (the old program's fixed name).
    blanks: String,
}

impl Default for View {
    fn default() -> Self {
        Self {
            mode: Mode::default(),
            totals: false,
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
    settings: TableSettings,
    /// The applied settings (Analysis → Settings): scan settings and, when they change, the quant ions.
    applied: GcmsSettings,
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
    pub fn use_settings(&mut self, applied: &GcmsSettings) {
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
            ui.label("Runs in memory at once");
            ui.add(egui::DragValue::new(&mut self.parallel).range(1..=64));
            ui.end_row();
            ui.label("Quant ions per substance");
            ui.add(egui::DragValue::new(&mut s.quant.ion_number).range(1..=20));
            ui.end_row();
            ui.label("Min. distance of quant ions");
            ui.add(
                egui::DragValue::new(&mut s.quant.min_dist)
                    .range(0.0..=100.0)
                    .speed(0.1),
            );
            ui.end_row();
            ui.label("Quant ion quality cut-off");
            ui.add(
                egui::DragValue::new(&mut s.quant.quality_cutoff)
                    .range(0.0..=100.0)
                    .speed(0.1),
            );
            ui.end_row();
            ui.label("Replicate groups (regex on file name)");
            ui.add(
                egui::TextEdit::singleline(&mut s.replicates).hint_text(r"e.g. ^(.+)_\d+\.cdf$"),
            )
            .on_hover_text("Group = first capture group, else the whole match. A file without a match is a group of its own, so it passes any reproducibility filter. Empty: one group.");
            ui.end_row();
            ui.label("Min. reproducibility in a group");
            ui.add(egui::Slider::new(&mut s.min_reproducibility, 0.0..=1.0));
            ui.end_row();
            ui.label("Min. quant ions found");
            ui.add(egui::DragValue::new(&mut s.min_ions).range(1..=20));
            ui.end_row();
            ui.label("Min. signal-to-noise of an ion");
            ui.add(
                egui::DragValue::new(&mut s.min_sn)
                    .range(0.0..=1e6)
                    .speed(0.1),
            );
            ui.end_row();
            ui.label("Max. discrepancy index");
            let mut limit = s.max_discrepancy.is_finite();
            ui.horizontal(|ui| {
                ui.checkbox(&mut limit, "");
                let mut v = if limit { s.max_discrepancy } else { 0.0 };
                ui.add_enabled(limit, egui::DragValue::new(&mut v).range(0.0..=1e6).speed(0.1))
                    .on_hover_text("A value counts only if the mean peak quality of its ions is at most this.");
                s.max_discrepancy = if limit { v } else { f32::INFINITY };
            });
            ui.end_row();
            ui.label("Untargeted");
            ui.checkbox(&mut self.options.untargeted, "also list unidentified components")
                .on_hover_text("Components without library hit are pooled over all runs into unknown substances. Applies to the next scan.");
            ui.end_row();
            ui.label("Rescan missing ions");
            ui.checkbox(&mut self.options.rescan, "search the ion traces")
                .on_hover_text("Searches the window of the RI difference around the expected positions. Applies to the next scan, with the quant ions chosen now.");
            ui.end_row();
        });
        crate::settings::ranges(
            ui,
            "QUANT_EXCL_IONS",
            "Ion ranges (m/z from–to) never picked as quant ions.",
            &mut s.quant.excluded,
        );
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
                // Not while the library is rebuilt for new search settings: it would be scaled differently.
                let ready =
                    lib.filter(|l| !self.paths.is_empty() && l.settings == self.applied.search());
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
                if lib.is_some_and(|l| l.settings != self.applied.search()) {
                    ui.label("Library is being prepared for the new search settings.");
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
            "{} runs scanned in {:.2?} against {}; {} unknown substances, {} ions added by the rescan",
            scan.runs.len(),
            scan.took,
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
                        if let Some(note) = export_button(ui, shown) {
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
                    );
                    if self.view.totals && self.totals.0 != self.generation {
                        self.totals = (self.generation, shown.total_intensities());
                    }
                    let totals = self.view.totals.then_some(self.totals.1.as_slice());
                    if let Some(act) = table_view(ui, shown, &self.norm, &self.view, totals) {
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
        let e = frag_export(&runs, &compounds);
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

/// The export has the values the table shows. Returns a note after a save attempt.
fn export_button(ui: &mut egui::Ui, t: &Table) -> Option<String> {
    if !ui.button("Export CSV…").clicked() {
        return None;
    }
    let p = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name("results.csv")
        .save_file()?;
    let saved = std::fs::File::create(&p)
        .map(std::io::BufWriter::new)
        .and_then(|f| t.write_csv(f));
    Some(match saved.map(|()| t.write_groups_file(&p)) {
        Ok(Ok(None)) => format!("Saved {}", p.display()),
        Ok(Ok(Some(g))) => format!("Saved {} and {}", p.display(), g.display()),
        Ok(Err(e)) => format!("Saved {}, but not the replicate groups: {e}", p.display()),
        Err(e) => format!("Export failed: {e}"),
    })
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

fn view_controls(ui: &mut egui::Ui, v: &mut View, edits: &mut Edits, groups: &[&str]) {
    ui.horizontal(|ui| {
        for (mode, title) in [
            (Mode::Values, "Values"),
            (Mode::Times, "Retention times (min)"),
            (Mode::Groups, "Group statistics"),
        ] {
            let r = ui.selectable_value(&mut v.mode, mode, title);
            if mode == Mode::Groups {
                r.on_hover_text("Per replicate group: mean, relative standard error (standard error / mean), hits (#) and the p-value of Student's t-test against the blank group.");
            }
        }
        ui.add_enabled(
            v.mode == Mode::Values,
            egui::Checkbox::new(&mut v.totals, "Total sample intensities"),
        )
        .on_hover_text(
            "Adds a first row with the sum of each run's values (a missing value counts as 0).",
        )
        .on_disabled_hover_text("Only in the Values view.");
        if v.mode == Mode::Groups {
            ui.label("Blank group");
            ui.add(egui::TextEdit::singleline(&mut v.blanks).desired_width(80.0))
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
            *edits = Edits::default();
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
        out.push(number(s.mean, 0));
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
            .body(|body| {
                let lead = totals.is_some() as usize;
                body.rows(18.0, t.rows.len() + lead, |mut row| {
                    let Some(i) = row.index().checked_sub(lead) else {
                        total_row(&mut row, totals.unwrap_or_default());
                        return;
                    };
                    let r = &t.rows[i];
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
                    let lead = [r.name.clone(), format!("{:.2}", r.rt_s / 60.0), ions.join(" ")];
                    for text in lead.into_iter().chain(cells(t, v, &groups, blank, i)) {
                        row.col(|ui| {
                            ui.label(text);
                        });
                    }
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
                ui.label("Internal standard");
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
                    });
            });
        }
        Kind::Pool => {
            ui.horizontal(|ui| {
                ui.label("Pool runs are the replicate group");
                ui.add(egui::TextEdit::singleline(&mut n.pool_group).desired_width(80.0));
                ui.label("nearest pool runs");
                ui.add(egui::DragValue::new(&mut n.pool_nearest).range(1..=99));
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
/// shortcut: files are found by name, so equal names in different folders are not told apart.
fn run_times(t: &Table, runs: &[RunResult]) -> (Vec<i64>, usize) {
    let by_name: std::collections::HashMap<_, _> = runs
        .iter()
        .filter_map(|r| Some((r.path.file_name()?.to_string_lossy().into_owned(), &r.path)))
        .collect();
    let mut missing = 0;
    let times = t
        .runs
        .iter()
        .map(|name| {
            by_name
                .get(name)
                .and_then(|p| std::fs::metadata(p).ok()?.modified().ok())
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
        b.use_settings(&applied);
        assert_eq!(b.settings.quant.ion_number, 5);
        // An edit in the window survives settings that leave QUANT_* alone ...
        b.settings.quant.ion_number = 2;
        applied.ls_ri_diff += 1.0;
        b.use_settings(&applied);
        assert_eq!(b.settings.quant.ion_number, 2);
        // ... and the next change of QUANT_* replaces it.
        applied.quant_ion_number = 4;
        b.use_settings(&applied);
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
    }
}
