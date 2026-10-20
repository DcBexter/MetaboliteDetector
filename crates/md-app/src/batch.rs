//! Batch quantification window: scans many runs on a worker thread against the open library and
//! shows the substance × run table (`md_core::batch`). Filters and groups apply without a rescan.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

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
            .open(&mut open)
            .default_size([900.0, 600.0])
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
            self.shown = None;
            self.times = None;
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
                    let n = &self.norm;
                    let times = self.times.as_ref().map(|t| t.0.as_slice());
                    let shown = apply(t, &mode(n, t, times.unwrap_or_default()), n.z_score);
                    self.shown = Some((n.clone(), shown));
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
                    table_view(ui, shown, &self.norm);
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
    Some(
        match std::fs::File::create(&p)
            .map(std::io::BufWriter::new)
            .and_then(|f| t.write_csv(f))
        {
            Ok(()) => format!("Saved {}", p.display()),
            Err(e) => format!("Export failed: {e}"),
        },
    )
}

/// Rows are virtualised, so only the visible ones are laid out each frame.
/// shortcut: all columns are laid out per visible row; virtualise columns too if studies with
/// thousands of runs feel slow.
fn table_view(ui: &mut egui::Ui, t: &Table, norm: &Norm) {
    let what = match (norm.kind, norm.z_score) {
        (Kind::None, false) => "summed quant ion integrals",
        (_, false) => "normalised quant ion integrals",
        (Kind::None, true) => "z-scores of the summed integrals",
        (_, true) => "z-scores of the normalised integrals",
    };
    ui.label(format!(
        "{} substances × {} runs; values are {what}",
        t.rows.len(),
        t.runs.len()
    ));
    egui::ScrollArea::horizontal().show(ui, |ui| {
        // Rows have a fixed height, so cells must not wrap.
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        TableBuilder::new(ui)
            .id_salt("batch table")
            .striped(true)
            .column(Column::auto().at_least(160.0))
            .columns(Column::auto(), 2)
            .columns(Column::auto().at_least(70.0), t.runs.len())
            .header(36.0, |mut h| {
                for title in ["Substance", "RT (min)", "Quant ions"] {
                    h.col(|ui| {
                        ui.strong(title);
                    });
                }
                for (run, group) in t.runs.iter().zip(&t.groups) {
                    h.col(|ui| {
                        ui.vertical(|ui| {
                            ui.strong(run);
                            ui.label(group);
                        });
                    });
                }
            })
            .body(|body| {
                body.rows(18.0, t.rows.len(), |mut row| {
                    let r = &t.rows[row.index()];
                    let ions: Vec<String> = r.ions.iter().map(u32::to_string).collect();
                    row.col(|ui| {
                        ui.label(&r.name);
                    });
                    row.col(|ui| {
                        ui.label(format!("{:.2}", r.rt_s / 60.0));
                    });
                    row.col(|ui| {
                        ui.label(ions.join(" "));
                    });
                    for v in &r.values {
                        row.col(|ui| {
                            ui.label(v.map_or("–".into(), |v| format!("{v:.0}")));
                        });
                    }
                });
            });
    });
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
}
