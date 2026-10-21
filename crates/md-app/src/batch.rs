//! Batch quantification window: scans many runs on a worker thread against the open library and
//! shows the substance × run table (`md_core::batch`). Filters and groups apply without a rescan.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use crate::errors;
use crate::lang::{self, tf};
use crate::layout::ToolWindow;
use crate::settings::row;
use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::batch::{
    Completed, Found, RunResult, SampleGroups, SampleTable, ScanSettings, ShiftReference, Table,
    TableSettings, complete_scan, path_key, run_key, scan_runs, table,
};
use md_core::batch::{
    cell_ions, group_label, group_names, group_stats_path, hit_list, ions_path, settings_path,
};
use md_core::library::LibraryCompound;
use md_core::mid::frag_export;
use md_core::normalize::{Normalization, POOL_GROUP, Pool, apply, pool_mapping_text};
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
    /// Per run: library index of a substance -> position in `RunResult::found`, so the details
    /// behind a cell are found without searching every frame.
    index: Vec<HashMap<usize, usize>>,
}

/// See `Scan::index`.
fn hit_index(runs: &[RunResult]) -> Vec<HashMap<usize, usize>> {
    (runs.iter())
        .map(|r| {
            r.found
                .iter()
                .enumerate()
                .map(|(i, f)| (f.compound, i))
                .collect()
        })
        .collect()
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
    /// Write the export for a spreadsheet with German locale (`;`, decimal comma).
    excel_de: bool,
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
            excel_de: false,
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
    /// Save the components of each run as `.cmp` (old Compound Detection Wizard).
    save_cmp: Option<md_core::batch::CmpTarget>,
}

pub struct Batch {
    pub open: bool,
    paths: Vec<PathBuf>,
    options: ScanOptions,
    /// Runs in memory at once; each run is analysed in parallel anyway.
    parallel: usize,
    pub settings: TableSettings,
    /// The sample table behind `settings.sample_groups`.
    samples: Option<Samples>,
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

/// A loaded sample table and the two columns that give run and group.
struct Samples {
    /// File name of the table, for the settings file.
    name: String,
    table: SampleTable,
    file: usize,
    group: usize,
    error: Option<String>,
}

impl Samples {
    /// File and chosen columns, for the settings file; none while the table gives no groups.
    fn source(&self) -> Option<String> {
        let col = |i: usize| &self.table.columns[i];
        self.error.is_none().then(|| {
            format!(
                "{} (files {}, groups {})",
                self.name,
                col(self.file),
                col(self.group)
            )
        })
    }

    /// Reads the chosen columns into `groups`; a bad table leaves no groups.
    fn apply(&mut self, groups: &mut SampleGroups) {
        match self.table.groups(self.file, self.group) {
            Ok(g) => (*groups, self.error) = (g, None),
            Err(e) => (*groups, self.error) = (SampleGroups::new(), Some(e)),
        }
    }
}

/// "Groups from table": replicate groups from a CSV/TSV table of file name → group, for runs whose
/// names do not tell them. Runs the table lacks keep the pattern's group.
fn samples_row(
    ui: &mut egui::Ui,
    samples: &mut Option<Samples>,
    groups: &mut SampleGroups,
    paths: &[PathBuf],
    note: &mut String,
) {
    let hint = lang::t(
        "Load a CSV, TSV or semicolon table with a column of run file names and a column \
         with the replicate group (e.g. samples.txt: CDF_FILE, TIME_POINT). It wins over \
         the pattern for the runs it names; other runs keep the pattern's group.",
        "Eine CSV-, TSV- oder Semikolon-Tabelle mit einer Spalte der Lauf-Dateinamen und einer \
         Spalte der Replikatgruppe laden (z. B. samples.txt: CDF_FILE, TIME_POINT). Sie hat Vorrang \
         vor dem Muster für die Läufe, die sie nennt; andere Läufe behalten die Gruppe des Musters.",
    );
    ui.label(lang::t(
        "Replicate groups from table",
        "Replikatgruppen aus Tabelle",
    ))
    .on_hover_text(hint);
    ui.horizontal_wrapped(|ui| {
        if ui
            .button(lang::t("Table…", "Tabelle…"))
            .on_hover_text(hint)
            .clicked()
            && let Some(p) = rfd::FileDialog::new()
                .add_filter(
                    lang::t("Sample table", "Probentabelle"),
                    &["txt", "tsv", "csv"],
                )
                .pick_file()
        {
            match std::fs::read_to_string(&p)
                .map_err(|e| e.to_string())
                .and_then(|t| SampleTable::parse(&t))
            {
                Ok(table) => {
                    let (file, group) = table.default_columns();
                    let mut s = Samples {
                        name: p.file_name().unwrap_or_default().to_string_lossy().into(),
                        table,
                        file,
                        group,
                        error: None,
                    };
                    s.apply(groups);
                    *samples = Some(s);
                }
                Err(e) => *note = format!("{}: {}", p.display(), errors::text(&e)),
            }
        }
        let Some(s) = samples else { return };
        let mut changed = false;
        let cols = [
            (
                "File column",
                lang::t("File column", "Dateispalte"),
                lang::t(
                    "Column with the run file names.",
                    "Spalte mit den Lauf-Dateinamen.",
                ),
                &mut s.file,
            ),
            (
                "Group column",
                lang::t("Group column", "Gruppenspalte"),
                lang::t(
                    "Column with the replicate group.",
                    "Spalte mit der Replikatgruppe.",
                ),
                &mut s.group,
            ),
        ];
        for (id, label, hover, col) in cols {
            egui::ComboBox::new(id, label)
                .selected_text(&s.table.columns[*col])
                .show_ui(ui, |ui| {
                    for (i, c) in s.table.columns.iter().enumerate() {
                        changed |= ui.selectable_value(col, i, c).changed();
                    }
                })
                .response
                .on_hover_text(hover);
        }
        if changed {
            s.apply(groups);
        }
        match &s.error {
            Some(e) => {
                ui.colored_label(ui.visuals().error_fg_color, errors::text(e));
            }
            None => {
                let found = (paths.iter())
                    .filter(|p| groups.contains_key(&run_key(&p.to_string_lossy())))
                    .count();
                let n = paths.len();
                ui.label(tf!(
                    "{found} of {n} runs in the table",
                    "{found} von {n} Läufen in der Tabelle"
                ));
            }
        }
        if ui
            .button(lang::t("Remove", "Entfernen"))
            .on_hover_text(lang::t(
                "Forget the table; the groups come from the pattern again.",
                "Die Tabelle verwerfen; die Gruppen kommen wieder aus dem Muster.",
            ))
            .clicked()
        {
            (*samples, *groups) = (None, SampleGroups::new());
        }
    });
    ui.end_row();
}

/// What a run without a group by hand gets: the table's group, else the pattern's, else (like
/// `table`) a group of its own or "all runs".
fn hint_text(from_table: Option<&String>, has_pattern: bool, has_table: bool) -> &str {
    match (from_table, has_pattern, has_table) {
        (Some(g), ..) => g,
        (None, true, _) => lang::t("pattern", "Muster"),
        (None, false, true) => lang::t("own group", "eigene Gruppe"),
        (None, false, false) => lang::t("all runs", "alle Läufe"),
    }
}

/// "Group per run": type a replicate group for a run; it wins over the table and the pattern. An
/// empty field takes the run back to them. `table` only fills the field's hint.
fn manual_groups_row(
    ui: &mut egui::Ui,
    manual: &mut SampleGroups,
    table: &SampleGroups,
    has_pattern: bool,
    paths: &[PathBuf],
) {
    let hint = lang::t(
        "Set the replicate group of single runs by hand. It wins over the table and the \
         pattern; an empty field leaves the run to them.",
        "Die Replikatgruppe einzelner Läufe von Hand setzen. Sie hat Vorrang vor Tabelle und \
         Muster; ein leeres Feld überlässt den Lauf ihnen.",
    );
    ui.label(lang::t("Group per run", "Gruppe pro Lauf"))
        .on_hover_text(hint);
    let set = (paths.iter())
        .filter(|p| manual.contains_key(&path_key(p)))
        .count();
    let n = paths.len();
    egui::CollapsingHeader::new(tf!("{set} of {n} runs set", "{set} von {n} Läufen gesetzt"))
        .id_salt("manual groups")
        .show(ui, |ui| {
            if ui
                .add_enabled(
                    !manual.is_empty(),
                    egui::Button::new(lang::t("Clear all", "Alle leeren")),
                )
                .clicked()
            {
                manual.clear();
            }
            let height = ui.text_style_height(&egui::TextStyle::Body) + 8.0;
            egui::ScrollArea::vertical()
                .id_salt("manual groups list")
                .max_height(200.0)
                .show_rows(ui, height, paths.len(), |ui, rows| {
                    for path in &paths[rows] {
                        let name = path.file_name().unwrap_or_default().to_string_lossy();
                        let key = path_key(path);
                        let mut text = manual.get(&key).cloned().unwrap_or_default();
                        ui.horizontal(|ui| {
                            let label = ui.label(name.as_ref());
                            let edit = egui::TextEdit::singleline(&mut text)
                                .desired_width(100.0)
                                .hint_text(hint_text(
                                    table.get(&run_key(&name)),
                                    has_pattern,
                                    !table.is_empty(),
                                ));
                            if ui.add(edit).labelled_by(label.id).changed() {
                                if text.trim().is_empty() {
                                    manual.remove(&key);
                                } else {
                                    manual.insert(key, text);
                                }
                            }
                        });
                    }
                });
        });
    ui.end_row();
}

/// "Save components as .cmp": off, next to each run, or into a chosen folder.
fn save_cmp_row(ui: &mut egui::Ui, target: &mut Option<md_core::batch::CmpTarget>) {
    use md_core::batch::CmpTarget;
    let hint = lang::t(
        "Also save the components of each run as run.cmp, with the RI of the search, \
         like the old Compound Detection Wizard. Saved runs load in \"Results of the \
         active run\" and scan faster as batch input. Applies to the next scan.",
        "Die Komponenten jedes Laufs zusätzlich als run.cmp speichern, mit dem RI der Suche, \
         wie der alte Compound Detection Wizard. Gespeicherte Läufe laden in „Ergebnisse des \
         aktiven Laufs“ und werden als Batch-Eingabe schneller durchsucht. Gilt ab dem nächsten Scan.",
    );
    let name = ui
        .label(lang::t(
            "Save components as .cmp",
            "Komponenten als .cmp speichern",
        ))
        .on_hover_text(hint);
    ui.horizontal(|ui| {
        let mut on = target.is_some();
        if ui
            .add(egui::Checkbox::without_text(&mut on))
            .labelled_by(name.id)
            .on_hover_text(hint)
            .changed()
        {
            *target = on.then_some(CmpTarget::BesideRun);
        }
        let Some(t) = target else { return };
        ui.label(match t {
            CmpTarget::BesideRun => lang::t("next to each run", "neben jedem Lauf").into(),
            CmpTarget::Dir(d) => d.display().to_string(),
        });
        if ui.button(lang::t("Folder…", "Ordner…")).clicked()
            && let Some(d) = rfd::FileDialog::new().pick_folder()
        {
            *t = CmpTarget::Dir(d);
        }
        if matches!(t, CmpTarget::Dir(_))
            && ui
                .button(lang::t("Next to runs", "Neben die Läufe"))
                .clicked()
        {
            *t = CmpTarget::BesideRun;
        }
    });
    ui.end_row();
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
                save_cmp: None,
            },
            settings: TableSettings::default(),
            samples: None,
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
            index: hit_index(&runs),
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
        // As tall as the room below its spot (`tool_window` caps it), so the result table gets the
        // room the settings leave.
        egui::Window::new(lang::t("Batch quantification", "Batch-Quantifizierung"))
            .id(egui::Id::unique("Batch quantification"))
            .tool_window(ctx, 2, [900.0, 2000.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui, lib, cal, reference));
        self.open = open;
    }

    fn poll(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.frag_job
            && let Some(r) = crate::poll_job(rx, lang::t("MID export", "MID-Export"))
        {
            (self.note, self.frag_job) = (r.unwrap_or_else(|e| e), None);
        }
        let Some(job) = &self.job else { return };
        match job.rx.try_recv() {
            Ok(Some((runs, mut errors, done))) => {
                errors.extend(done.errors);
                self.scan = Some(Scan {
                    lib: job.lib.clone(),
                    index: hit_index(&runs),
                    runs: Arc::new(runs),
                    extra: Arc::new(done.extra),
                    rescanned: done.rescanned,
                    errors,
                    took: job.started.elapsed(),
                });
                (self.table, self.note) = (None, String::new());
            }
            Ok(None) => self.note = lang::t("Cancelled.", "Abgebrochen.").into(),
            Err(mpsc::TryRecvError::Empty) => {
                // Keeps the progress bar moving.
                ctx.request_repaint_after(Duration::from_millis(200));
                return;
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.note = lang::t("Batch crashed.", "Batch abgestürzt.").into()
            }
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
            save_cmp: self.options.save_cmp.clone(),
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
    pub fn settings_grid(&mut self, ui: &mut egui::Ui) {
        let s = &mut self.settings;
        egui::Grid::new("batch settings").show(ui, |ui| {
            row(
                ui,
                lang::t(
                    "Replicate groups (regex on file name)",
                    "Replikatgruppen (Regex auf Dateinamen)",
                ),
                "",
                lang::t(
                    "Group = first capture group, else the whole match. A file without a match is a \
                     group of its own, so it passes any reproducibility filter. Empty: one group.",
                    "Gruppe = erste Erfassungsgruppe, sonst der ganze Treffer. Eine Datei ohne Treffer \
                     ist eine eigene Gruppe und besteht so jeden Reproduzierbarkeitsfilter. Leer: eine Gruppe.",
                ),
                egui::TextEdit::singleline(&mut s.replicates)
                    .hint_text(lang::t(r"e.g. ^(.+)_\d+\.cdf$", r"z. B. ^(.+)_\d+\.cdf$")),
            );
            samples_row(
                ui,
                &mut self.samples,
                &mut s.sample_groups,
                &self.paths,
                &mut self.note,
            );
            manual_groups_row(
                ui,
                &mut s.manual_groups,
                &s.sample_groups,
                !s.replicates.is_empty(),
                &self.paths,
            );
            row(
                ui,
                lang::t(
                    "Min. reproducibility in a group",
                    "Min. Reproduzierbarkeit in einer Gruppe",
                ),
                "",
                lang::t(
                    "Keep a substance only if some replicate group quantifies it in at least this share \
                     of its runs.",
                    "Eine Substanz nur behalten, wenn eine Replikatgruppe sie in mindestens diesem Anteil \
                     ihrer Läufe quantifiziert.",
                ),
                egui::Slider::new(&mut s.min_reproducibility, 0.0..=1.0),
            );
            row(
                ui,
                lang::t("Untargeted", "Ungezielt"),
                "",
                lang::t(
                    "Also list components without a library hit: they are pooled over all runs into \
                     unknown substances. Applies to the next scan.",
                    "Auch Komponenten ohne Bibliothekstreffer aufführen: Sie werden über alle Läufe zu \
                     unbekannten Substanzen zusammengefasst. Gilt ab dem nächsten Scan.",
                ),
                egui::Checkbox::without_text(&mut self.options.untargeted),
            );
            save_cmp_row(ui, &mut self.options.save_cmp);
        });
        crate::guided::advanced(ui, "batch settings", |ui| {
            egui::Grid::new("batch advanced").show(ui, |ui| {
                row(
                    ui,
                    lang::t("Runs in memory at once", "Läufe gleichzeitig im Speicher"),
                    "",
                    lang::t(
                        "How many files are read at the same time; each run is analysed in parallel anyway. \
                         Fewer need less memory.",
                        "Wie viele Dateien gleichzeitig gelesen werden; jeder Lauf wird ohnehin parallel analysiert. \
                         Weniger braucht weniger Speicher.",
                    ),
                    egui::DragValue::new(&mut self.parallel).range(1..=64),
                );
                row(
                    ui,
                    lang::t("Quant ions per substance", "Quant-Ionen pro Substanz"),
                    "QUANT_ION_NUMBER",
                    lang::t(
                        "How many ions of each substance are integrated and summed. Only for substances whose library entry names no quant ions.",
                        "Wie viele Ionen jeder Substanz integriert und summiert werden. Nur für Substanzen, deren Bibliothekseintrag keine Quant-Ionen nennt.",
                    ),
                    egui::DragValue::new(&mut s.quant.ion_number).range(1..=20),
                );
                row(
                    ui,
                    lang::t("Min. distance of quant ions", "Min. Abstand der Quant-Ionen"),
                    "QUANT_MIN_DIST_IONS",
                    lang::t(
                        "Quant ions of one substance are preferably at least this many m/z apart. Only for substances whose library entry names no quant ions.",
                        "Quant-Ionen einer Substanz liegen möglichst mindestens so viele m/z auseinander. Nur für Substanzen, deren Bibliothekseintrag keine Quant-Ionen nennt.",
                    ),
                    egui::DragValue::new(&mut s.quant.min_dist)
                        .range(0.0..=100.0)
                        .speed(0.1),
                );
                row(
                    ui,
                    lang::t("Quant ion quality cut-off", "Qualitätsschwelle der Quant-Ionen"),
                    "QUANT_QUALITY_CUTOFF",
                    lang::t(
                        "Ions whose peak quality is below this (lower is better) are picked first. Only for substances whose library entry names no quant ions.",
                        "Ionen, deren Peak-Qualität darunter liegt (niedriger ist besser), werden zuerst gewählt. Nur für Substanzen, deren Bibliothekseintrag keine Quant-Ionen nennt.",
                    ),
                    egui::DragValue::new(&mut s.quant.quality_cutoff)
                        .range(0.0..=100.0)
                        .speed(0.1),
                );
                row(
                    ui,
                    lang::t("Min. quant ions found", "Min. gefundene Quant-Ionen"),
                    "",
                    lang::t(
                        "A value counts only if at least this many of the substance's quant ions were found.",
                        "Ein Wert zählt nur, wenn mindestens so viele Quant-Ionen der Substanz gefunden wurden.",
                    ),
                    egui::DragValue::new(&mut s.min_ions).range(1..=20),
                );
                row(
                    ui,
                    lang::t(
                        "Min. signal-to-noise of an ion",
                        "Min. Signal-Rausch-Verhältnis eines Ions",
                    ),
                    "",
                    lang::t(
                        "An ion counts only if its peak has at least this signal-to-noise ratio (old \
                         \"Req. S/N\").",
                        "Ein Ion zählt nur, wenn sein Peak mindestens dieses Signal-Rausch-Verhältnis hat (alt \
                         „Req. S/N“).",
                    ),
                    egui::DragValue::new(&mut s.min_sn)
                        .range(0.0..=1e6)
                        .speed(0.1),
                );
                let hint = lang::t(
                    "A value counts only if the mean peak quality of its ions is at most this. \
                     Off: no limit.",
                    "Ein Wert zählt nur, wenn die mittlere Peak-Qualität seiner Ionen höchstens so groß ist. \
                     Aus: keine Grenze.",
                );
                let name = ui
                    .label(lang::t("Max. discrepancy index", "Max. Diskrepanzindex"))
                    .on_hover_text(hint);
                let mut limit = s.max_discrepancy.is_finite();
                ui.horizontal(|ui| {
                    ui.add(egui::Checkbox::without_text(&mut limit))
                        .labelled_by(name.id)
                        .on_hover_text(hint);
                    let mut v = if limit { s.max_discrepancy } else { 0.0 };
                    ui.add_enabled(limit, egui::DragValue::new(&mut v).range(0.0..=1e6).speed(0.1))
                        .labelled_by(name.id)
                        .on_hover_text(hint)
                        .on_disabled_hover_text(lang::t(
                            "Tick the box to set a limit.",
                            "Das Kästchen ankreuzen, um eine Grenze zu setzen.",
                        ));
                    s.max_discrepancy = if limit { v } else { f32::INFINITY };
                });
                ui.end_row();
                row(
                    ui,
                    lang::t("Rescan missing ions", "Fehlende Ionen erneut suchen"),
                    "",
                    lang::t(
                        "Search the ion traces again in the window of the RI difference around the expected \
                         positions. Applies to the next scan, with the quant ions chosen now.",
                        "Die Ionenspuren im Fenster der RI-Differenz um die erwarteten Positionen erneut \
                         durchsuchen. Gilt ab dem nächsten Scan, mit den jetzt gewählten Quant-Ionen.",
                    ),
                    egui::Checkbox::without_text(&mut self.options.rescan),
                );
            });
            crate::settings::ranges(
                ui,
                lang::t(
                    "Ions avoided as quant ions",
                    "Als Quant-Ionen gemiedene Ionen",
                ),
                lang::t(
                    "Ion ranges (m/z from–to) not picked as quant ions while other ions are left. Only for substances whose library entry names no quant ions.\nOld key: QUANT_EXCL_IONS",
                    "Ionenbereiche (m/z von–bis), die nicht als Quant-Ionen gewählt werden, solange andere Ionen übrig sind. Nur für Substanzen, deren Bibliothekseintrag keine Quant-Ionen nennt.\nAlter Schlüssel: QUANT_EXCL_IONS",
                ),
                &mut s.quant.excluded,
            );
        });
    }

    /// Whether `lib` was prepared for the search settings and the RT-only mode that apply now.
    fn library_ready(&self, lib: &Library) -> bool {
        lib.settings == self.applied.search() && lib.prepared.is_rt_only() == self.rt_only
    }

    /// Progress and Cancel while a scan runs, else Start and what it will use.
    pub fn job_row(
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
                    true => lang::t("pooling and rescanning…", "zusammenfassen und erneut suchen…")
                        .into(),
                    false => tf!("{d} / {} runs", "{d} / {} Läufe", job.total),
                };
                ui.add(
                    egui::ProgressBar::new(d as f32 / job.total.max(1) as f32)
                        .text(text)
                        .desired_width(300.0),
                );
                if ui.button(lang::t("Cancel", "Abbrechen")).clicked() {
                    job.cancel.store(true, Relaxed);
                }
            }
            None => {
                // Not while the library is rebuilt for new search settings or a switch of RT-only: it
                // would be scaled differently, and the calibration is hidden meanwhile.
                let ready = lib.filter(|l| !self.paths.is_empty() && self.library_ready(l));
                if ui
                    .add_enabled(ready.is_some(), egui::Button::new(lang::t("Start", "Starten")))
                    .clicked()
                    && let Some(l) = ready
                {
                    self.start(ui.ctx(), l, cal, reference);
                }
                match lib {
                    Some(l) => ui.label(match cal.is_some() {
                        true => tf!("Library {}, RI calibrated", "Bibliothek {}, RI kalibriert", l.name),
                        false => tf!(
                            "Library {}, no RI calibration",
                            "Bibliothek {}, keine RI-Kalibrierung",
                            l.name
                        ),
                    }),
                    None => ui.label(lang::t(
                        "Open a library first (File → Open library).",
                        "Erst eine Bibliothek öffnen (Datei → Bibliothek öffnen).",
                    )),
                };
                if lib.is_some_and(|l| !self.library_ready(l)) {
                    ui.label(lang::t(
                        "Library is being prepared for the new search or RT-only settings.",
                        "Die Bibliothek wird für die neuen Such- oder Nur-RT-Einstellungen vorbereitet.",
                    ));
                }
                if self.paths.is_empty() {
                    ui.label(lang::t("Add runs first.", "Läufe fehlen."));
                }
                ui.label(errors::text(&self.note));
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
        self.runs_row(ui);
        // A short screen scrolls the settings, so the result table below stays on it. Measured on
        // the screen: the window's own height follows its content and would shrink with the cap.
        egui::ScrollArea::vertical()
            .id_salt("batch settings")
            .max_height(ui.ctx().content_rect().height() * 0.4)
            .show(ui, |ui| self.settings_grid(ui));
        self.job_row(ui, lib, cal, reference);
        self.results(ui);
    }

    /// Runs scanned and errors of the last scan, for the guided Start step.
    pub fn summary(&self) -> Option<(usize, usize)> {
        self.scan.as_ref().map(|s| (s.runs.len(), s.errors.len()))
    }

    /// Whether runs are chosen; the guided steps go on once they are.
    pub fn has_runs(&self) -> bool {
        !self.paths.is_empty()
    }

    /// The guided component detection is the batch that saves the components; it starts with that on.
    /// Returns whether it was off, so a cancelled detection can switch it off again.
    pub fn save_cmp_by_default(&mut self) -> bool {
        let off = self.options.save_cmp.is_none();
        (self.options.save_cmp).get_or_insert(md_core::batch::CmpTarget::BesideRun);
        off
    }

    pub fn stop_saving_cmp(&mut self) {
        self.options.save_cmp = None;
    }

    #[cfg(test)]
    pub fn saves_cmp(&self) -> bool {
        self.options.save_cmp.is_some()
    }

    /// The `.cmp` choice alone, for the guided component detection.
    pub fn cmp_grid(&mut self, ui: &mut egui::Ui) {
        egui::Grid::new("batch cmp").show(ui, |ui| save_cmp_row(ui, &mut self.options.save_cmp));
    }

    pub fn runs_row(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui
                .button(lang::t("Add runs…", "Läufe hinzufügen…"))
                .clicked()
                && let Some(p) = rfd::FileDialog::new()
                    .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                    .add_filter(
                        lang::t(
                            "Old MetaboliteDetector (.bin + .idx)",
                            "Alter MetaboliteDetector (.bin + .idx)",
                        ),
                        &["bin"],
                    )
                    .add_filter(
                        lang::t(
                            "FastFlight raw data (MsData.FFT)",
                            "FastFlight-Rohdaten (MsData.FFT)",
                        ),
                        &["FFT", "fft"],
                    )
                    .add_filter(
                        lang::t("Saved components (.cmp)", "Gespeicherte Komponenten (.cmp)"),
                        &["cmp", "CMP"],
                    )
                    .pick_files()
            {
                self.paths.extend(p);
                self.paths.sort();
                self.paths.dedup();
            }
            if ui.button(lang::t("Clear", "Leeren")).clicked() {
                self.paths.clear();
            }
            ui.label(tf!("{} runs", "{} Läufe", self.paths.len()));
        });
    }

    /// The last scan: its summary, normalisation, exports and the table.
    fn results(&mut self, ui: &mut egui::Ui) {
        let Some(scan) = &self.scan else { return };
        if self.table.as_ref().is_none_or(|t| t.0 != self.settings) {
            let t = table(&scan.runs, &scan.lib.compounds, &scan.extra, &self.settings)
                .map_err(|e| e.to_string());
            self.table = Some((self.settings.clone(), t));
            (self.shown, self.times, self.norm.edits) = (None, None, Edits::default());
        }
        ui.separator();
        ui.label(tf!(
            "{} runs scanned in {} against {}; {} unknown substances, {} ions added by the rescan",
            "{} Läufe in {} gegen {} durchsucht; {} unbekannte Substanzen, {} Ionen durch die erneute Suche ergänzt",
            scan.runs.len(),
            crate::view::took(scan.took),
            scan.lib.name,
            scan.extra.len(),
            scan.rescanned
        ));
        for e in &scan.errors {
            ui.colored_label(ui.visuals().error_fg_color, errors::text(e));
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
                    ui.horizontal_wrapped(|ui| {
                        let ex = Export {
                            scan,
                            settings: &self.settings,
                            extra: [
                                vec![(
                                    "GROUP_SOURCE",
                                    self.settings.group_source(
                                        self.samples.as_ref().and_then(Samples::source).as_deref(),
                                        &scan
                                            .runs
                                            .iter()
                                            .map(|r| r.path.clone())
                                            .collect::<Vec<_>>(),
                                    ),
                                )],
                                shown_settings(&self.norm, shown, t.runs.len()),
                                scan_settings(scan),
                            ]
                            .concat(),
                            ions: self.view.export_ions,
                            excel_de: self.view.excel_de,
                            blank: &self.view.blanks,
                            group_stats: !self.norm.z_score,
                        };
                        if let Some(note) = export_button(ui, shown, &ex) {
                            self.note = note;
                        }
                        let times = self.times.as_ref().map(|t| t.0.as_slice());
                        pool_mapping_button(ui, t, &self.norm, times, &mut self.note);
                        if self.frag_job.is_none()
                            && let Some(rx) = frag_button(ui, shown, scan)
                        {
                            (self.note, self.frag_job) = (
                                lang::t("Exporting MIDs…", "MIDs werden exportiert…").into(),
                                Some(rx),
                            );
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
                let e = errors::text(e);
                ui.colored_label(
                    ui.visuals().error_fg_color,
                    tf!("Replicate pattern: {e}", "Replikatmuster: {e}"),
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
        .button(lang::t("Export MIDs (FRAG)…", "MIDs exportieren (FRAG)…"))
        .on_hover_text(lang::t(
            "Total ion counts and isotopomer shares of the fragments in each substance's FRAG feature, e.g. (174/8)(261/5).",
            "Gesamtionenzahlen und Isotopomer-Anteile der Fragmente im FRAG-Merkmal jeder Substanz, z. B. (174/8)(261/5).",
        ))
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
            Ok(()) => tf!("Saved {}", "Gespeichert: {}", p.display()),
            Err(err) => tf!("Export failed: {err}", "Export fehlgeschlagen: {err}"),
        };
        for m in &e.errors {
            note += &format!("\n{}", errors::text(m));
        }
        match e.warnings.len() {
            0 => {}
            1 => note += &format!("\n{}", errors::text(&e.warnings[0])),
            n => {
                note += &tf!(
                    "\n{n} substances have no fragment (FRAG) defined.",
                    "\nFür {n} Substanzen ist kein Fragment (FRAG) definiert."
                )
            }
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
    /// Convert the files for a spreadsheet with German locale (`md_core::csv::to_excel_de`).
    excel_de: bool,
    /// The replicate group the group statistics test against (if the table has it).
    blank: &'a str,
    /// Write the group statistics; not for z-scores, which the view does not show them for either.
    group_stats: bool,
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

/// How the runs were scanned, for the settings file next to an export.
fn scan_settings(scan: &Scan) -> Vec<(&'static str, String)> {
    vec![
        ("LIBRARY", scan.lib.name.clone()),
        ("RT_ONLY", scan.lib.prepared.is_rt_only().to_string()),
        ("SEARCH_CUTOFF", scan.lib.settings.cutoff.to_string()),
        ("RUNS_SCANNED", scan.runs.len().to_string()),
        ("UNKNOWN_SUBSTANCES", scan.extra.len().to_string()),
        ("IONS_ADDED_BY_RESCAN", scan.rescanned.to_string()),
    ]
}

/// The export has the values the table shows, plus the replicate groups and the settings used
/// next to it, and with `ex.ions` the raw integral of every quant ion. Returns a note after a
/// save attempt.
fn export_button(ui: &mut egui::Ui, t: &Table, ex: &Export) -> Option<String> {
    if !ui
        .button(lang::t("Export CSV…", "CSV exportieren…"))
        .clicked()
    {
        return None;
    }
    let p = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name("results.csv")
        .save_file()?;
    Some(match save_results(&p, t, ex) {
        Ok(note) => note,
        Err(e) => tf!("Export failed: {e}", "Export fehlgeschlagen: {e}"),
    })
}

/// Removes an older file next to an export that would not match the new results file; only a
/// file that is there but cannot be removed is an error (a script would read its old numbers).
fn remove_stale(path: &Path) -> std::io::Result<Option<PathBuf>> {
    match std::fs::remove_file(path) {
        Err(e) if e.kind() != std::io::ErrorKind::NotFound => Err(e),
        _ => Ok(None),
    }
}

/// Files are written in the standard format first; the German spreadsheet variant is a rewrite.
fn german(path: &Path, text_cols: usize, ex: &Export) -> std::io::Result<()> {
    if ex.excel_de {
        md_core::csv::excel_de_file(path, text_cols)?;
    }
    Ok(())
}

/// Writes the results CSV `p` and the files next to it; the note names what was saved and what
/// could not be.
fn save_results(p: &Path, t: &Table, ex: &Export) -> std::io::Result<String> {
    let create = |p: &Path| std::fs::File::create(p).map(std::io::BufWriter::new);
    t.write_csv(create(p)?)?;
    german(p, 3, ex)?;
    let sp = settings_path(p);
    let group_stats = if ex.group_stats {
        t.write_group_stats_file(p, Some(ex.blank))
    } else {
        remove_stale(&group_stats_path(p))
    };
    // (text columns for the German variant, name in the note, file)
    let quant_ions = lang::t("quant ions", "Quant-Ionen");
    let mut results = vec![
        (
            2,
            lang::t("replicate groups", "Replikatgruppen"),
            t.write_groups_file(p).transpose(),
        ),
        (
            3,
            lang::t("group statistics", "Gruppenstatistik"),
            group_stats.transpose(),
        ),
        (
            2,
            lang::t("settings", "Einstellungen"),
            Some(
                create(&sp)
                    .and_then(|f| ex.settings.write_csv(&ex.extra, f))
                    .map(|()| sp),
            ),
        ),
    ];
    if !ex.ions {
        results.push((2, quant_ions, remove_stale(&ions_path(p)).transpose()));
    } else {
        let ip = ions_path(p);
        let runs: Vec<&RunResult> = (t.sources.iter())
            .filter_map(|i| ex.scan.runs.get(*i))
            .collect();
        let r = create(&ip).and_then(|f| t.write_ions_csv(&runs, ex.settings, f));
        results.push((2, quant_ions, Some(r.map(|()| ip))));
    }
    let mut saved = vec![p.display().to_string()];
    let mut failed = Vec::new();
    for (cols, what, r) in results {
        let r = r.map(|r| r.and_then(|f| german(&f, cols, ex).map(|()| f)));
        match r {
            Some(Ok(path)) => saved.push(path.display().to_string()),
            Some(Err(e)) => failed.push(format!("{what}: {e}")),
            None => {}
        }
    }
    let mut text = tf!("Saved {}", "Gespeichert: {}", saved.join(", "));
    if !failed.is_empty() {
        text += &tf!(
            "; not saved: {}",
            "; nicht gespeichert: {}",
            failed.join("; ")
        );
    }
    Ok(text)
}

/// The runs as the user left them (removed and moved), with their file dates.
fn edited_runs<'a>(
    t: &'a Table,
    n: &Norm,
    times: Option<&[i64]>,
) -> (std::borrow::Cow<'a, Table>, Vec<i64>) {
    use std::borrow::Cow;
    let cols = n.edits.cols.as_deref();
    let base = cols.map_or(Cow::Borrowed(t), |c| Cow::Owned(t.select_runs(c)));
    // The file dates follow the runs that are left.
    let times = match (times, cols) {
        (Some(ts), Some(cols)) => cols
            .iter()
            .map(|c| ts.get(*c).copied().unwrap_or(0))
            .collect(),
        (Some(ts), None) => ts.to_vec(),
        (None, _) => Vec::new(),
    };
    (base, times)
}

/// The table as the user's edits and the normalisation leave it: runs are removed and moved
/// first, substances are removed last.
fn shown_table(t: &Table, n: &Norm, times: Option<&[i64]>) -> Table {
    let (base, times) = edited_runs(t, n, times);
    apply(&base, &mode(n, &base, &times), n.z_score).without_rows(&n.edits.rows)
}

/// The old program's "Export pool mapping": which pool runs each run was divided by, for the runs
/// the table shows. A save attempt leaves a note.
fn pool_mapping_button(
    ui: &mut egui::Ui,
    t: &Table,
    n: &Norm,
    times: Option<&[i64]>,
    note: &mut String,
) {
    if n.kind != Kind::Pool
        || !ui
            .button(lang::t(
                "Export pool mapping…",
                "Pool-Zuordnung exportieren…",
            ))
            .clicked()
    {
        return;
    }
    let (base, times) = edited_runs(t, n, times);
    let Normalization::Pool(p) = mode(n, &base, &times) else {
        return;
    };
    let Some(mut path) = rfd::FileDialog::new()
        .add_filter(lang::t("Text", "Text"), &["txt"])
        .set_file_name("pool_mapping.txt")
        .save_file()
    else {
        return;
    };
    if path.extension().is_none_or(|e| e != "txt") {
        path.as_mut_os_string().push(".txt");
    }
    *note = match std::fs::write(&path, pool_mapping_text(&base, &p)) {
        Ok(()) => tf!("Saved {}", "Gespeichert: {}", path.display()),
        Err(e) => tf!("Export failed: {e}", "Export fehlgeschlagen: {e}"),
    };
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
            (Mode::Values, lang::t("Values", "Werte")),
            (
                Mode::Times,
                lang::t("Retention times (min)", "Retentionszeiten (min)"),
            ),
            (Mode::Groups, lang::t("Group statistics", "Gruppenstatistik")),
        ] {
            let locked = mode == Mode::Groups && z_scored;
            let r = ui
                .add_enabled_ui(!locked, |ui| ui.selectable_value(&mut v.mode, mode, title))
                .inner;
            if mode == Mode::Groups {
                r.on_hover_text(lang::t(
                    "Per replicate group: mean, relative standard error (standard error / mean), hits (#) and the p-value of Student's t-test against the blank group.",
                    "Pro Replikatgruppe: Mittelwert, relativer Standardfehler (Standardfehler / Mittelwert), Treffer (#) und der p-Wert des Student-t-Tests gegen die Blindwertgruppe.",
                ))
                .on_disabled_hover_text(lang::t(
                    "Not with z-scores: a missing value cannot count as 0 there, and \"#\" (values above 0) means nothing. Turn the z-score off; the group statistics then use the normalised values, like the old program.",
                    "Nicht mit z-Scores: Dort kann ein fehlender Wert nicht als 0 zählen, und „#“ (Werte über 0) sagt nichts aus. Den z-Score ausschalten; die Gruppenstatistik nutzt dann die normalisierten Werte, wie das alte Programm.",
                ));
            }
        }
    });
    ui.horizontal(|ui| {
        ui.add_enabled(
            v.mode == Mode::Values,
            egui::Checkbox::new(
                &mut v.totals,
                lang::t("Total sample intensities", "Probensummen"),
            ),
        )
        .on_hover_text(lang::t(
            "Adds a first row with the sum of each run's values (a missing value counts as 0).",
            "Fügt eine erste Zeile mit der Summe der Werte jedes Laufs hinzu (ein fehlender Wert zählt als 0).",
        ))
        .on_disabled_hover_text(lang::t(
            "Only in the Values view.",
            "Nur in der Ansicht Werte.",
        ));
        ui.add_enabled(
            v.mode == Mode::Values,
            egui::Checkbox::new(
                &mut v.ion_rows,
                lang::t("Quant ion lines", "Quant-Ion-Zeilen"),
            ),
        )
        .on_hover_text(lang::t(
            "Adds a line under each substance for every quant ion, with its raw integral per run.",
            "Fügt unter jeder Substanz eine Zeile pro Quant-Ion hinzu, mit seinem Rohintegral pro Lauf.",
        ))
        .on_disabled_hover_text(lang::t(
            "Only in the Values view.",
            "Nur in der Ansicht Werte.",
        ));
        ui.checkbox(
            &mut v.export_ions,
            lang::t("Export quant ions", "Quant-Ionen-Export"),
        )
        .on_hover_text(lang::t(
            "Also writes results.ions.csv with the raw integral of every quant ion.",
            "Schreibt zusätzlich results.ions.csv mit dem Rohintegral jedes Quant-Ions.",
        ));
        ui.checkbox(&mut v.excel_de, lang::t("Excel (German)", "Excel (deutsch)"))
            .on_hover_text(lang::t(
                "Exports with semicolons and decimal commas, for a spreadsheet with German locale. Off: commas and decimal points; open those with locale English (USA) in the import dialog.",
                "Exportiert mit Semikolons und Dezimalkommas, für eine Tabellenkalkulation mit deutschem Gebietsschema. Aus: Kommas und Dezimalpunkte; diese im Importdialog mit dem Gebietsschema Englisch (USA) öffnen.",
            ));
        if v.mode == Mode::Groups {
            let blank = ui.label(lang::t("Blank group", "Blindwertgruppe"));
            ui.add(egui::TextEdit::singleline(&mut v.blanks).desired_width(80.0))
                .labelled_by(blank.id)
                .on_hover_text(lang::t(
                    "Replicate group of the blank runs; every other group gets a p-value against it.",
                    "Replikatgruppe der Blindwertläufe; jede andere Gruppe bekommt einen p-Wert dagegen.",
                ));
            if !groups.contains(&v.blanks.as_str()) {
                ui.weak(lang::t(
                    "no such group, so no p-values",
                    "keine solche Gruppe, daher keine p-Werte",
                ));
            }
        }
        if !edits.is_empty()
            && ui
                .button(lang::t("Restore removed", "Entferntes wiederherstellen"))
                .on_hover_text(lang::t(
                    "Brings back the removed runs and substances and the original order of the runs.",
                    "Holt die entfernten Läufe und Substanzen und die ursprüngliche Reihenfolge der Läufe zurück.",
                ))
                .clicked()
        {
            edits.restore();
        }
    });
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
pub(crate) fn significant(v: f64) -> String {
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
            .map(|(r, g)| (r.clone(), group_label(g).to_string()))
            .collect();
    }
    let stats = [
        lang::t("Mean", "Mittelwert"),
        lang::t("Rel. Std. Err.", "Rel. Std.-Fehler"),
        "#",
        lang::t("P (blank)", "P (Blindwert)"),
    ];
    let per_group = if blank { 4 } else { 3 };
    (groups.iter())
        .flat_map(|g| {
            stats[..per_group]
                .iter()
                .map(move |s| (s.to_string(), group_label(g).to_string()))
        })
        .chain([(
            "P (ANOVA)".to_string(),
            lang::t("all groups", "alle Gruppen").to_string(),
        )])
        .collect()
}

/// The cells of one substance in "Group statistics": per group mean, relative standard error,
/// hits and (with a blank group) the p-value against the blanks; last the ANOVA p-value over all
/// the groups, blanks included, like the old "p-value" column.
fn group_cells(t: &Table, row: usize, groups: &[&str], blank: Option<&str>) -> Vec<String> {
    let p = |p: Option<f64>| p.map_or("–".into(), |p| format!("{p:.4}"));
    let stats = t.group_stats(row, groups, blank);
    let mut out = Vec::new();
    for (k, s) in stats.groups.iter().enumerate() {
        out.push(significant(s.mean));
        out.push(number(s.rel_std_err, 3));
        out.push(s.hits.to_string());
        if blank.is_some() {
            out.push(p(stats.p_blank[k]));
        }
    }
    out.push(p(stats.anova));
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
                .on_hover_text(lang::t("Move run left", "Lauf nach links"))
                .clicked()
        {
            act = Some(Act::MoveRun(j, j - 1));
        }
        if j + 1 < n
            && ui
                .small_button("›")
                .on_hover_text(lang::t("Move run right", "Lauf nach rechts"))
                .clicked()
        {
            act = Some(Act::MoveRun(j, j + 1));
        }
        if ui
            .small_button("×")
            .on_hover_text(lang::t(
                "Remove this run from the results (Restore removed brings it back)",
                "Diesen Lauf aus den Ergebnissen entfernen (Entferntes wiederherstellen holt ihn zurück)",
            ))
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
        let run = *t.sources.get(col)?;
        let hit = *self.scan.index.get(run)?.get(&compound)?;
        self.scan.runs.get(run)?.found.get(hit)
    }

    /// The library hits of the component behind a cell: the assigned one first (red), then the
    /// alternatives, like the old tooltip.
    fn hits_ui(&self, ui: &mut egui::Ui, t: &Table, col: usize, compound: usize) {
        let Some(f) = self.found(t, col, compound) else {
            ui.label(lang::t(
                "Not found in this run.",
                "In diesem Lauf nicht gefunden.",
            ));
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
                .on_hover_text(lang::t(
                    "Remove this substance from the results (Restore removed brings it back)",
                    "Diese Substanz aus den Ergebnissen entfernen (Entferntes wiederherstellen holt sie zurück)",
                ))
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
            ui.weak(tf!("    ion {mass}", "    Ion {mass}"));
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
        (Mode::Times, ..) => lang::t(
            "retention times of the hits in minutes",
            "Retentionszeiten der Treffer in Minuten",
        ),
        (Mode::Groups, ..) => lang::t(
            "per replicate group; a missing value counts as 0",
            "pro Replikatgruppe; ein fehlender Wert zählt als 0",
        ),
        (_, Kind::None, false) => lang::t(
            "summed quant ion integrals",
            "summierte Integrale der Quant-Ionen",
        ),
        (_, _, false) => lang::t(
            "normalised quant ion integrals",
            "normalisierte Integrale der Quant-Ionen",
        ),
        (_, Kind::None, true) => lang::t(
            "z-scores of the summed integrals",
            "z-Scores der summierten Integrale",
        ),
        (_, _, true) => lang::t(
            "z-scores of the normalised integrals",
            "z-Scores der normalisierten Integrale",
        ),
    };
    let (rows, runs) = (t.rows.len(), t.runs.len());
    ui.label(tf!(
        "{rows} substances × {runs} runs; values are {what}",
        "{rows} Substanzen × {runs} Läufe; Werte sind {what}"
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
    // Both bars belong to this area, so they sit at the window edges. The table's own vertical bar
    // would be at the far end of the wide table, out of sight; the header then scrolls with the rows.
    egui::ScrollArea::both()
        .auto_shrink([false, true])
        .show(ui, |ui| {
            // Rows have a fixed height, so cells must not wrap.
            ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
            TableBuilder::new(ui)
                .id_salt("batch table")
                .vscroll(false)
                .striped(true)
                .column(Column::exact(28.0))
                .column(Column::auto().at_least(160.0))
                .columns(Column::auto(), 2)
                .columns(Column::auto().at_least(70.0), heads.len())
                .header(54.0, |mut h| {
                    h.col(|_| {});
                    for title in [
                        lang::t("Substance", "Substanz"),
                        "RT (min)",
                        lang::t("Quant ions", "Quant-Ionen"),
                    ] {
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
        ui.strong(lang::t(
            "Total sample intensity",
            "Gesamtintensität der Probe",
        ));
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
        ui.label(lang::t("Normalise by", "Normieren nach"));
        for (kind, title) in [
            (Kind::None, lang::t("nothing", "nichts")),
            (
                Kind::InternalStandard,
                lang::t("internal standard", "internem Standard"),
            ),
            (Kind::TotalSignal, lang::t("total signal", "Gesamtsignal")),
            (Kind::Pool, lang::t("pool", "Pool")),
        ] {
            ui.selectable_value(&mut n.kind, kind, title);
        }
        ui.checkbox(&mut n.z_score, lang::t("z-score", "z-Score"));
    });
    match n.kind {
        Kind::InternalStandard => {
            ui.horizontal(|ui| {
                let label = ui.label(lang::t("Internal standard", "Interner Standard"));
                let name = |c: Option<usize>| {
                    t.rows
                        .iter()
                        .find(|r| Some(r.compound) == c)
                        .map_or(lang::t("choose…", "wählen…"), |r| r.name.as_str())
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
                let group = ui.label(lang::t(
                    "Pool runs are the replicate group",
                    "Pool-Läufe sind die Replikatgruppe",
                ));
                ui.add(egui::TextEdit::singleline(&mut n.pool_group).desired_width(80.0))
                    .labelled_by(group.id);
                let nearest = ui.add(egui::DragValue::new(&mut n.pool_nearest).range(1..=99));
                let label = ui.label(lang::t("nearest pool runs", "nächste Pool-Läufe"));
                nearest.labelled_by(label.id);
                ui.checkbox(&mut n.pool_median, lang::t("median", "Median"));
            })
            .response
            .on_hover_text(lang::t(
                "Each run is divided by the mean or median of the nearest pool runs in time (file date). Pool runs keep their raw values.",
                "Jeder Lauf wird durch Mittelwert oder Median der zeitlich nächsten Pool-Läufe (Dateidatum) geteilt. Pool-Läufe behalten ihre Rohwerte.",
            ));
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
        ui.label(tf!(
            "No file date for {missing} runs; they count as time 0.",
            "Kein Dateidatum für {missing} Läufe; sie zählen als Zeitpunkt 0."
        ));
    }
    if !t.groups.contains(&n.pool_group) {
        ui.label(lang::t(
            "No run is in the pool group: every normalised value is empty.",
            "Kein Lauf ist in der Pool-Gruppe: Jeder normalisierte Wert ist leer.",
        ));
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
    fn anova_column_covers_all_groups_and_stands_last() {
        // Hand calculation (md-core stats test): groups of 1 2 3 / 2 3 4 / 5 6 7 give p = 27 / 4096.
        let v = [1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let t = Table {
            runs: (0..9).map(|i| i.to_string()).collect(),
            groups: ["A", "A", "A", "B", "B", "B", "C", "C", "C"]
                .map(String::from)
                .to_vec(),
            rows: vec![md_core::batch::Row {
                compound: 0,
                name: "s0".into(),
                rt_s: 60.0,
                ions: vec![73],
                values: v.map(Some).to_vec(),
                rts_s: vec![Some(60.0); 9],
            }],
            ..Default::default()
        };
        let groups = ["A", "B", "C"];
        let cells = group_cells(&t, 0, &groups, None);
        assert_eq!(cells.len(), 3 * 3 + 1);
        assert_eq!(cells.last().unwrap(), "0.0066");
        // With a blank group its column comes between the group's, the ANOVA still last and over all.
        let cells = group_cells(&t, 0, &groups, Some("C"));
        assert_eq!(cells.len(), 3 * 4 + 1);
        assert_eq!(cells.last().unwrap(), "0.0066");
        let view = View {
            mode: Mode::Groups,
            ..Default::default()
        };
        assert_eq!(heads(&t, &view, &groups, true).len(), cells.len());
        assert_eq!(heads(&t, &view, &groups, false).len(), 3 * 3 + 1);
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
        // One group: nothing to compare, p = 1 like the old code.
        assert_eq!(group_cells(&t, 0, &["A"], None)[3], "1.0000");
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
            index: hit_index(std::slice::from_ref(&run)),
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
        std::fs::remove_dir_all(&dir).ok(); // leftovers of an aborted run
        std::fs::create_dir_all(&dir).unwrap();
        let csv = dir.join("results.csv");
        let settings = TableSettings::default();
        let mut ex = Export {
            scan: &scan,
            settings: &settings,
            extra: [vec![("NORMALISATION", "none".into())], scan_settings(&scan)].concat(),
            ions: false,
            excel_de: false,
            blank: "BLANK",
            group_stats: true,
        };
        let note = save_results(&csv, &t, &ex).unwrap();
        assert!(note.contains("results.settings.csv"), "{note}");
        assert!(!dir.join("results.ions.csv").exists());
        let text = std::fs::read_to_string(dir.join("results.settings.csv")).unwrap();
        assert!(text.contains("\nNORMALISATION,none\n"), "{text}");
        assert!(text.contains("\nRUNS_SCANNED,1\n"), "{text}");
        assert!(text.contains("\nIONS_ADDED_BY_RESCAN,0\n"), "{text}");

        ex.ions = true;
        save_results(&csv, &t, &ex).unwrap();
        let ions = std::fs::read_to_string(dir.join("results.ions.csv")).unwrap();
        assert_eq!(
            ions,
            "\u{feff}Substance,Quant_ion,a.cdf\nGlucose,73,500\nGlucose,147,100\n"
        );

        // The German spreadsheet variant converts every file of the export.
        ex.excel_de = true;
        save_results(&csv, &t, &ex).unwrap();
        let ions = std::fs::read_to_string(dir.join("results.ions.csv")).unwrap();
        assert_eq!(
            ions,
            "\u{feff}Substance;Quant_ion;a.cdf\nGlucose;73;500\nGlucose;147;100\n"
        );
        let main = std::fs::read_to_string(&csv).unwrap();
        assert!(
            main.starts_with("\u{feff}Substance;RT_s;Quant_ions;"),
            "{main}"
        );
        let text = std::fs::read_to_string(dir.join("results.settings.csv")).unwrap();
        assert!(text.contains("\nNORMALISATION;none\n"), "{text}");
        ex.excel_de = false;

        // Group statistics only for a table that has groups, and not for z-scores (the view has none).
        let grouped = Table {
            groups: vec!["ctrl".into()],
            ..t.clone()
        };
        save_results(&csv, &grouped, &ex).unwrap();
        assert!(dir.join("results.groupstats.csv").exists());
        ex.group_stats = false;
        save_results(&csv, &grouped, &ex).unwrap();
        assert!(!dir.join("results.groupstats.csv").exists());
        // Without groups there is no file either, and an old one goes.
        std::fs::write(dir.join("results.groupstats.csv"), "old").unwrap();
        ex.group_stats = true;
        save_results(&csv, &t, &ex).unwrap();
        assert!(!dir.join("results.groupstats.csv").exists());

        // An old file that cannot be removed (here a directory of that name) is reported.
        std::fs::create_dir_all(dir.join("results.groupstats.csv").join("x")).unwrap();
        ex.group_stats = false;
        let note = save_results(&csv, &grouped, &ex).unwrap();
        assert!(note.contains("not saved: group statistics"), "{note}");
        // The same for an old ions file.
        std::fs::remove_file(dir.join("results.ions.csv")).unwrap();
        std::fs::create_dir_all(dir.join("results.ions.csv").join("x")).unwrap();
        let note = save_results(&csv, &grouped, &ex).unwrap();
        assert!(note.contains("quant ions:"), "{note}");
        std::fs::remove_dir_all(&dir).ok();
    }
}
