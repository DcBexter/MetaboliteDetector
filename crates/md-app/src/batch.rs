//! Batch quantification window: scans many runs on a worker thread against the open library and
//! shows the substance × run table (`md_core::batch`). Filters and groups apply without a rescan.

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::batch::{RunResult, ScanSettings, Table, TableSettings, scan_runs, table};
use md_core::normalize::{Normalization, POOL_GROUP, Pool, apply};
use md_core::ri::RiCalibration;
use md_core::search::SearchSettings;

use crate::identify::Library;
use crate::spawn;

type Scanned = Option<(Vec<RunResult>, Vec<String>)>;

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
    runs: Vec<RunResult>,
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

pub struct Batch {
    pub open: bool,
    paths: Vec<PathBuf>,
    /// Runs in memory at once; each run is analysed in parallel anyway.
    parallel: usize,
    settings: TableSettings,
    /// Peak detection and deconvolution settings of the app (Analysis → Settings), set every frame.
    pub analysis: crate::settings::Key,
    job: Option<Job>,
    scan: Option<Scan>,
    /// Table of `scan` for these settings, or the pattern error; rebuilt when the settings change.
    table: Option<(TableSettings, Result<Table, String>)>,
    norm: Norm,
    /// `table` normalised as `norm` says.
    shown: Option<(Norm, Table)>,
    note: String,
}

impl Default for Batch {
    fn default() -> Self {
        Self {
            open: false,
            paths: Vec::new(),
            // 15 example runs: 1.1 s with 1, 0.17 s with 4 (90 MiB), 0.10 s with 15 (189 MiB).
            parallel: 4,
            settings: TableSettings::default(),
            analysis: Default::default(),
            job: None,
            scan: None,
            table: None,
            norm: Norm::default(),
            shown: None,
            note: String::new(),
        }
    }
}

impl Batch {
    pub fn window(
        &mut self,
        ctx: &egui::Context,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
    ) {
        self.poll(ctx);
        let mut open = self.open;
        egui::Window::new("Batch quantification")
            .open(&mut open)
            .default_size([900.0, 600.0])
            .show(ctx, |ui| self.contents(ui, lib, cal));
        self.open = open;
    }

    fn poll(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        match job.rx.try_recv() {
            Ok(Some((runs, errors))) => {
                self.scan = Some(Scan {
                    lib: job.lib.clone(),
                    runs,
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

    fn start(&mut self, ctx: &egui::Context, lib: &Arc<Library>, cal: Option<&RiCalibration>) {
        let (done, cancel) = (
            Arc::new(AtomicUsize::new(0)),
            Arc::new(AtomicBool::new(false)),
        );
        // The library was prepared with the app's search settings, which differ only in the cutoff.
        let s = ScanSettings {
            peaks: self.analysis.0,
            deconv: self.analysis.1.clone(),
            search: SearchSettings {
                cutoff: SearchSettings::default().cutoff,
                ..lib.settings.clone()
            },
            calibration: cal.cloned(),
        };
        let (paths, l, d, c, parallel) = (
            self.paths.clone(),
            lib.clone(),
            done.clone(),
            cancel.clone(),
            self.parallel,
        );
        let rx = spawn(ctx, move || {
            scan_runs(&paths, &l.prepared, &s, parallel, &d, &c)
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

    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
    ) {
        ui.horizontal(|ui| {
            if ui.button("Add runs…").clicked()
                && let Some(p) = rfd::FileDialog::new()
                    .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                    .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
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
        });
        ui.horizontal(|ui| match &self.job {
            Some(job) => {
                let d = job.done.load(Relaxed);
                ui.add(
                    egui::ProgressBar::new(d as f32 / job.total.max(1) as f32)
                        .text(format!("{d} / {} runs", job.total))
                        .desired_width(300.0),
                );
                if ui.button("Cancel").clicked() {
                    job.cancel.store(true, Relaxed);
                }
            }
            None => {
                let ready = lib.filter(|_| !self.paths.is_empty());
                if ui
                    .add_enabled(ready.is_some(), egui::Button::new("Start"))
                    .clicked()
                    && let Some(l) = ready
                {
                    self.start(ui.ctx(), l, cal);
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
                ui.label(&self.note);
            }
        });
        let Some(scan) = &self.scan else { return };
        if self.table.as_ref().is_none_or(|t| t.0 != self.settings) {
            let t =
                table(&scan.runs, &scan.lib.compounds, &self.settings).map_err(|e| e.to_string());
            self.table = Some((self.settings.clone(), t));
            self.shown = None;
        }
        ui.separator();
        ui.label(format!(
            "{} runs scanned in {:.2?} against {}",
            scan.runs.len(),
            scan.took,
            scan.lib.name
        ));
        for e in &scan.errors {
            ui.colored_label(ui.visuals().error_fg_color, e);
        }
        match self.table.as_ref().map(|t| &t.1) {
            Some(Ok(t)) => {
                norm_controls(ui, &mut self.norm, t);
                if self.shown.as_ref().is_none_or(|s| s.0 != self.norm) {
                    let n = &self.norm;
                    let shown = apply(t, &mode(n, t, &scan.runs), n.z_score);
                    self.shown = Some((n.clone(), shown));
                }
                if let Some((_, shown)) = &self.shown {
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

/// The old program orders runs by acquisition date; the file date stands in for it.
/// shortcut: files are found by name, so equal names in different folders are not told apart.
fn mode(n: &Norm, t: &Table, runs: &[RunResult]) -> Normalization {
    match n.kind {
        Kind::None => Normalization::None,
        Kind::TotalSignal => Normalization::TotalSignal,
        Kind::InternalStandard => Normalization::InternalStandard(
            t.rows
                .iter()
                .position(|r| Some(r.compound) == n.standard)
                .unwrap_or(usize::MAX),
        ),
        Kind::Pool => {
            let times = t
                .runs
                .iter()
                .map(|name| {
                    runs.iter()
                        .find(|r| {
                            r.path
                                .file_name()
                                .is_some_and(|f| f.to_string_lossy() == *name)
                        })
                        .and_then(|r| std::fs::metadata(&r.path).ok()?.modified().ok())
                        .and_then(|m| m.duration_since(std::time::UNIX_EPOCH).ok())
                        .map_or(0, |d| d.as_secs() as i64)
                })
                .collect();
            Normalization::Pool(Pool {
                group: n.pool_group.clone(),
                times,
                nearest: n.pool_nearest,
                median: n.pool_median,
            })
        }
    }
}
