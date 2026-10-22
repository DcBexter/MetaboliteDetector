//! MID window (F16, old `MIDWizard`): scans labelled (and unlabelled) runs on a worker thread
//! against the library's substances with fragments and shows their mass isotopomer distributions
//! (`md_core::mid`).

use std::path::PathBuf;
use std::sync::Arc;
use std::time::{Duration, Instant};

use crate::errors;
use crate::lang::{Msg, msg, t, tf};
use crate::layout::ToolWindow;
use eframe::egui;
use egui_plot::{Bar, BarChart, BoxElem, BoxPlot, BoxSpread, Legend, LegendGrouping, Plot};
use md_core::batch::{
    DEFAULT_PARALLEL, QuantSettings, ShiftReference, Study, StudyResult, run_study,
};
use md_core::isotope::IsotopeTable;
use md_core::library::LibraryCompound;
use md_core::mid::{FragmentMid, Method, Mids, determine, group_stat, mid_library};
use md_core::ri::RiCalibration;
use md_core::search::PreparedLibrary;
use md_core::settings::GcmsSettings;

use crate::Job;
use crate::identify::Library;

/// The old method page's choices, in its order.
#[derive(Clone, Copy, PartialEq)]
enum Kind {
    Uncorrected,
    Reference,
    Library,
    Formula,
}

struct Done {
    method: Method,
    /// The substances with fragments the MIDs refer to.
    lib: Vec<LibraryCompound>,
    /// One per fit: the labelled runs' names, or the group.
    columns: Vec<String>,
    /// Full path of each per-run column, for its hover; empty for the group methods.
    paths: Vec<String>,
    /// Position of each column when determined, before the user moved runs.
    runs: Vec<usize>,
    mids: Mids,
    /// Per-run methods: mean and standard error of the labelled runs per substance, fragment and
    /// isotopomer (old group `GL`); empty for the group methods.
    stats: Vec<Vec<Vec<(f64, f64)>>>,
    /// The old program's text exports: the MIDs, and for the per-run methods the statistics.
    text: String,
    stats_text: Option<String>,
    errors: Vec<Msg>,
    /// Runs that loaded and were scanned.
    scanned: usize,
    took: Duration,
}

/// The one column of the group methods, kept in English as the plot key; [`shown_name`] translates it.
const GROUP_COLUMN: &str = "labelled";

/// Column `i` as drawn: the group methods have no paths, and their column follows the language.
fn shown_name<'a>(columns: &'a [String], paths: &[String], i: usize) -> &'a str {
    match paths.is_empty() {
        true => t("labelled", "markiert"),
        false => &columns[i],
    }
}

/// A move in the result view: run columns (old "Change column order…") or substances with their
/// diagrams (old "Change order of diagrams…"), from one position to another.
#[derive(Clone, Copy)]
enum Move {
    Run(usize, usize),
    Compound(usize, usize),
}

impl Done {
    /// Applies a move to the table, the diagrams and the exports. The values stay as determined:
    /// the statistics keep summing the runs in their first order, so not even the last bit moves.
    fn apply(&mut self, m: Move) {
        match m {
            Move::Run(a, b) => {
                self.columns.swap(a, b);
                self.runs.swap(a, b);
                if !self.paths.is_empty() {
                    self.paths.swap(a, b);
                }
                (self.mids.compounds.iter_mut())
                    .flat_map(|c| &mut c.fragments)
                    .for_each(|f| f.fits.swap(a, b));
            }
            Move::Compound(a, b) => {
                self.mids.compounds.swap(a, b);
                if !self.stats.is_empty() {
                    self.stats.swap(a, b);
                }
                if self.stats_text.is_some() {
                    let mut first = self.mids.clone();
                    for f in first.compounds.iter_mut().flat_map(|c| &mut c.fragments) {
                        let mut fits = vec![None; f.fits.len()];
                        for (j, fit) in std::mem::take(&mut f.fits).into_iter().enumerate() {
                            fits[self.runs[j]] = fit;
                        }
                        f.fits = fits;
                    }
                    self.stats_text = Some(first.stats_text(&self.lib));
                }
            }
        }
        self.text = self.mids_text();
    }

    /// Statistics and text exports of the MIDs as determined.
    fn refresh(&mut self) {
        let per_run = self.method.per_run();
        self.runs = (0..self.columns.len()).collect();
        self.stats = (self.mids.compounds.iter())
            .filter(|_| per_run)
            .map(|c| {
                (c.fragments.iter())
                    .map(|f| (0..f.isotopomers).map(|k| group_stat(f, k)).collect())
                    .collect()
            })
            .collect();
        self.text = self.mids_text();
        self.stats_text = per_run.then(|| self.mids.stats_text(&self.lib));
    }

    fn mids_text(&self) -> String {
        match self.method.per_run() {
            true => (self.mids).run_text(&self.lib, self.method, &self.columns),
            false => self.mids.group_text(&self.lib),
        }
    }
}

struct Running {
    started: Instant,
    job: Job<Result<Done, Msg>>,
}

pub struct MidWindow {
    pub open: bool,
    labelled: Vec<PathBuf>,
    unlabelled: Vec<PathBuf>,
    kind: Kind,
    tracer: String,
    m1: bool,
    /// Natural abundance of the tracer's M+1 isotope in percent, for the M1 correction.
    abundance: f64,
    /// Tracer `abundance` was last filled in for.
    abundance_for: String,
    rescan: bool,
    /// The applied settings (Analysis → Settings), set every frame.
    pub applied: GcmsSettings,
    table: IsotopeTable,
    job: Option<Running>,
    result: Option<Done>,
    note: Msg,
}

impl Default for MidWindow {
    fn default() -> Self {
        Self {
            open: false,
            labelled: Vec::new(),
            unlabelled: Vec::new(),
            kind: Kind::Reference,
            tracer: "C".into(),
            m1: false,
            abundance: 0.0,
            abundance_for: String::new(),
            rescan: true,
            applied: GcmsSettings::default(),
            table: IsotopeTable::standard(),
            job: None,
            result: None,
            note: Msg::default(),
        }
    }
}

impl MidWindow {
    pub fn window(
        &mut self,
        ctx: &egui::Context,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        self.poll(ctx);
        let mut open = self.open;
        egui::Window::new(t(
            "Mass isotopomer distributions",
            "Massenisotopomerverteilungen",
        ))
        .id(egui::Id::unique("Mass isotopomer distributions"))
        .tool_window(ctx, 3, [800.0, 600.0])
        .open(&mut open)
        .show(ctx, |ui| self.contents(ui, lib, cal, reference));
        self.open = open;
    }

    fn poll(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        match job
            .job
            .poll(|| t("MID determination", "Die MID-Bestimmung").into())
        {
            Some(Ok(Ok(done))) => (self.result, self.note) = (Some(done), Msg::default()),
            Some(Ok(Err(e)) | Err(e)) => self.note = e,
            None => {
                ctx.request_repaint_after(Duration::from_millis(200));
                return;
            }
        }
        self.job = None;
    }

    fn method(&self) -> Method {
        match self.kind {
            Kind::Uncorrected => Method::Uncorrected,
            // The old wizard divides the abundance spin box by 100.
            Kind::Reference => Method::Reference {
                m1: if self.m1 { self.abundance / 100.0 } else { 0.0 },
            },
            Kind::Library => Method::Library,
            Kind::Formula => Method::Formula,
        }
    }

    fn start(
        &mut self,
        ctx: &egui::Context,
        lib: &Arc<Library>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        let mut s = self.applied.scan(cal.cloned());
        // Like the single-run search: only with an RI calibration.
        let shift_ref = cal.and(reference).map(|(n, c)| (n.to_owned(), c));
        let method = self.method();
        // Unlabelled runs only take part in the reference correction (old group `GUL`).
        let unlabelled = match self.kind {
            Kind::Reference => self.unlabelled.clone(),
            _ => Vec::new(),
        };
        let paths: Vec<PathBuf> = self.labelled.iter().chain(&unlabelled).cloned().collect();
        let columns: Vec<String> = match method.per_run() {
            true => (self.labelled.iter())
                // Old `QFileInfo::baseName`: the file name up to the first dot.
                .map(|p| {
                    let name = p.file_name().unwrap_or_default().to_string_lossy();
                    name.split('.').next().unwrap_or_default().to_owned()
                })
                .collect(),
            false => vec![GROUP_COLUMN.into()],
        };
        let column_paths: Vec<String> = match method.per_run() {
            true => (self.labelled.iter())
                .map(|p| p.display().to_string())
                .collect(),
            false => Vec::new(),
        };
        let (nl, nu) = (self.labelled.len(), unlabelled.len());
        let rescan = self.rescan;
        let (l, tracer, table) = (
            lib.clone(),
            self.tracer.trim().to_owned(),
            self.table.clone(),
        );
        let job = Job::spawn(ctx, move |progress, cancel| {
            let started = Instant::now();
            let mids_lib = mid_library(&l.compounds);
            if mids_lib.is_empty() {
                return Err(msg!(
                    "No substance of the library has a fragment with sum formula (MID fragments in the library editor).",
                    "Keine Substanz der Bibliothek hat ein Fragment mit Summenformel (MID-Fragmente im Bibliothekseditor)."
                ));
            }
            // The shift reference is an index into the scanned library, which keeps only
            // substances with fragments.
            s.shift = (shift_ref.as_ref())
                .and_then(|(n, cutoff)| ShiftReference::find(&mids_lib, n, *cutoff));
            let prepared = match l.prepared.is_rt_only() {
                true => PreparedLibrary::new_rt_only(&mids_lib, &l.settings),
                false => PreparedLibrary::new(&mids_lib, &l.settings),
            };
            let cancelled = || msg!("Cancelled.", "Abgebrochen.");
            // The quant ions are the fragment ions, so the rescan fills those that were not deconvoluted.
            let study = Study {
                runs: &paths,
                library: &mids_lib,
                prepared: &prepared,
                scan: &s,
                quant: &QuantSettings::default(),
                rescan,
                parallel: DEFAULT_PARALLEL,
            };
            let StudyResult { runs, errors, .. } =
                run_study(&study, progress, cancel).ok_or_else(cancelled)?;
            let mut errors: Vec<Msg> = errors.into_iter().map(Msg::from).collect();
            if let Some((n, _)) = shift_ref.filter(|_| s.shift.is_none()) {
                errors.push(msg!(
                    "RT shift reference {n} is not a substance with MID fragments; scanned without RT shift.",
                    "RT-Verschiebungsreferenz {n} ist keine Substanz mit MID-Fragmenten; ohne RT-Verschiebung gescannt."
                ));
            }
            // A run that fails to load is missing from `runs`; `input` is its position in `paths`.
            let position = |range: std::ops::Range<usize>| {
                (runs.iter().enumerate())
                    .filter(|(_, r)| range.contains(&r.input))
                    .map(|(i, _)| i)
                    .collect::<Vec<usize>>()
            };
            let labelled = position(0..nl);
            let unlabelled = position(nl..nl + nu);
            let loaded = |k: usize| runs.iter().any(|r| r.input == k);
            let mut mids = determine(
                &runs,
                &mids_lib,
                &labelled,
                &unlabelled,
                method,
                &tracer,
                &table,
            );
            if method.per_run() {
                // One column per labelled run as chosen, also for runs that did not load.
                for f in mids.compounds.iter_mut().flat_map(|c| &mut c.fragments) {
                    let mut fits = std::mem::take(&mut f.fits).into_iter();
                    f.fits = (0..nl)
                        .map(|k| loaded(k).then(|| fits.next().flatten()).flatten())
                        .collect();
                }
            }
            let mut done = Done {
                method,
                lib: mids_lib,
                columns,
                paths: column_paths,
                runs: Vec::new(),
                mids,
                stats: Vec::new(),
                text: String::new(),
                stats_text: None,
                errors,
                scanned: runs.len(),
                took: started.elapsed(),
            };
            done.refresh();
            Ok(done)
        });
        job.progress.set_total(nl + nu);
        self.job = Some(Running {
            started: Instant::now(),
            job,
        });
    }

    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        self.settings(ui);
        self.run_lists(ui);
        self.start_row(ui, lib, cal, reference);
        let Some(r) = &mut self.result else { return };
        ui.separator();
        result_view(ui, r, &mut self.note);
    }

    /// Progress and Cancel while it runs, else Start and what it will use.
    pub fn start_row(
        &mut self,
        ui: &mut egui::Ui,
        lib: Option<&Arc<Library>>,
        cal: Option<&RiCalibration>,
        reference: Option<(&str, f64)>,
    ) {
        ui.horizontal(|ui| match &self.job {
            Some(job) => {
                let (d, total) = (job.job.progress.done(), job.job.progress.total());
                ui.add(
                    egui::ProgressBar::new(d as f32 / total.max(1) as f32)
                        .text(tf!(
                            "{d} / {} runs, {:.0?}",
                            "{d} / {} Läufe, {:.0?}",
                            total,
                            job.started.elapsed()
                        ))
                        .desired_width(300.0),
                );
                if ui.button(t("Cancel", "Abbrechen")).clicked() {
                    job.job.cancel.cancel();
                }
            }
            None => {
                let problem = self.problem();
                // Not while the library is rebuilt for new search settings, as in the batch window.
                let ready =
                    lib.filter(|l| problem.is_none() && l.settings == self.applied.search());
                if ui
                    .add_enabled(ready.is_some(), egui::Button::new(t("Start", "Starten")))
                    .clicked()
                    && let Some(l) = ready
                {
                    self.start(ui.ctx(), l, cal, reference);
                }
                match (lib, problem) {
                    (None, _) => ui.label(t(
                        "Open a library first (File → Open library).",
                        "Erst eine Bibliothek öffnen (Datei → Bibliothek öffnen).",
                    )),
                    (_, Some(p)) => ui.label(p),
                    (Some(l), None) if l.settings != self.applied.search() => ui.label(t(
                        "Library is being prepared for the new search settings.",
                        "Die Bibliothek wird für die neuen Sucheinstellungen vorbereitet.",
                    )),
                    (Some(l), None) => ui.label(match cal.is_some() {
                        true => tf!(
                            "Library {}, RI calibrated",
                            "Bibliothek {}, RI-kalibriert",
                            l.name
                        ),
                        false => tf!(
                            "Library {}, no RI calibration",
                            "Bibliothek {}, keine RI-Kalibrierung",
                            l.name
                        ),
                    }),
                };
            }
        });
        if !self.note.is_empty() {
            ui.label(&*self.note);
        }
    }

    /// Runs scanned and problems of the last determination, for the guided Start step.
    pub fn summary(&self) -> Option<(usize, usize)> {
        let r = self.result.as_ref()?;
        Some((r.scanned, r.errors.len() + r.mids.errors.len()))
    }

    /// Why the MIDs cannot be determined yet.
    pub fn problem(&self) -> Option<&'static str> {
        let tracer = self.tracer.trim();
        if self.labelled.is_empty() {
            Some(t("Add labelled runs.", "Markierte Läufe hinzufügen."))
        } else if self.kind == Kind::Reference && self.unlabelled.is_empty() {
            Some(t(
                "Add unlabelled runs for the reference correction.",
                "Unmarkierte Läufe für die Referenzkorrektur hinzufügen.",
            ))
        } else if !tracer.is_empty() && self.table.isotopes(tracer).is_none() {
            Some(t(
                "The tracer is not an element of the isotope table.",
                "Der Tracer ist kein Element der Isotopentabelle.",
            ))
        } else if tracer.is_empty() && matches!(self.kind, Kind::Library | Kind::Formula) {
            Some(t(
                "This correction needs a tracer element.",
                "Diese Korrektur braucht ein Tracer-Element.",
            ))
        } else {
            None
        }
    }

    pub fn settings(&mut self, ui: &mut egui::Ui) {
        egui::Grid::new("mid settings").show(ui, |ui| {
            ui.label(t("Correction", "Korrektur"));
            ui.vertical(|ui| {
                // A grid column starts as narrow as its first frame; the short titles never wrap.
                ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
                for (kind, title, hint) in [
                    (Kind::Uncorrected, t("Uncorrected MIDs", "Unkorrigierte MIDs"), t("Mean and confidence interval of the labelled runs' normalised fragment spectra.", "Mittelwert und Konfidenzintervall der normierten Fragmentspektren der markierten Läufe.")),
                    (Kind::Reference, t("Corrected using reference", "Korrigiert mit Referenz"), t("Labelled runs against unlabelled runs of the same substances.", "Markierte Läufe gegen unmarkierte Läufe derselben Substanzen.")),
                    (Kind::Library, t("Corrected using library", "Korrigiert mit Bibliothek"), t("Every labelled run against the library spectrum, corrected for the tracer's natural M+1.", "Jeder markierte Lauf gegen das Bibliotheksspektrum, korrigiert um das natürliche M+1 des Tracers.")),
                    (Kind::Formula, t("Corrected using sum formulas", "Korrigiert mit Summenformeln"), t("Every labelled run against the natural isotope pattern of the fragment's sum formula.", "Jeder markierte Lauf gegen das natürliche Isotopenmuster der Summenformel des Fragments.")),
                ] {
                    ui.radio_value(&mut self.kind, kind, title).on_hover_text(hint);
                }
            });
            ui.end_row();
            let label = ui.label(t("Tracer element", "Tracer-Element"));
            ui.add(egui::TextEdit::singleline(&mut self.tracer).desired_width(40.0))
                .labelled_by(label.id);
            ui.end_row();
        });
        // Pre-fill the tracer's natural M+1 abundance whenever the tracer changes, also while
        // "Advanced" is closed.
        if self.kind == Kind::Reference && self.abundance_for != self.tracer {
            self.abundance_for = self.tracer.clone();
            if let Some(i) = self
                .table
                .isotopes(self.tracer.trim())
                .and_then(|i| i.get(1))
            {
                self.abundance = i.abundance * 100.0;
            }
        }
        crate::guided::advanced(ui, "mid settings", |ui| {
            egui::Grid::new("mid advanced").show(ui, |ui| {
                if self.kind == Kind::Reference {
                    ui.label(t("M1 correction", "M1-Korrektur"));
                    ui.horizontal(|ui| {
                        ui.checkbox(&mut self.m1, t("natural abundance", "natürliche Häufigkeit"));
                        ui.add_enabled(
                            self.m1,
                            egui::DragValue::new(&mut self.abundance)
                                .range(0.0..=100.0)
                                .speed(0.01)
                                .suffix(" %"),
                        );
                    });
                    ui.end_row();
                }
                let hint = t(
                    "Search the ion traces again for fragment ions the deconvolution missed",
                    "Die Ionenspuren erneut nach Fragment-Ionen durchsuchen, die die Dekonvolution übersehen hat",
                );
                let rescan = ui.label(t("Rescan missing ions", "Fehlende Ionen neu suchen")).on_hover_text(hint);
                ui.add(egui::Checkbox::without_text(&mut self.rescan))
                    .labelled_by(rescan.id)
                    .on_hover_text(hint);
                ui.end_row();
            });
        });
    }

    pub fn run_lists(&mut self, ui: &mut egui::Ui) {
        runs(
            ui,
            t("Labelled runs", "Markierte Läufe"),
            &mut self.labelled,
        );
        if self.kind == Kind::Reference {
            runs(
                ui,
                t("Unlabelled runs", "Unmarkierte Läufe"),
                &mut self.unlabelled,
            );
        }
    }
}

fn runs(ui: &mut egui::Ui, title: &str, paths: &mut Vec<PathBuf>) {
    ui.horizontal(|ui| {
        ui.label(format!("{title}: {}", paths.len()));
        if ui.button(t("Add…", "Hinzufügen…")).clicked()
            && let Some(p) = rfd::FileDialog::new()
                .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                .add_filter(
                    t(
                        "Old MetaboliteDetector (.bin + .idx)",
                        "Altes MetaboliteDetector (.bin + .idx)",
                    ),
                    &["bin"],
                )
                .pick_files()
        {
            paths.extend(p);
            paths.sort();
            paths.dedup();
        }
        if ui.button(t("Clear", "Leeren")).clicked() {
            paths.clear();
        }
    })
    .response
    .on_hover_text(
        paths
            .iter()
            .map(|p| p.display().to_string())
            .collect::<Vec<_>>()
            .join("\n"),
    );
}

/// Substances, fragments and per column the MID ± CI of each isotopomer and R², with a bar chart.
/// shortcut: every substance header is laid out each frame; fine for hundreds of substances.
fn result_view(ui: &mut egui::Ui, r: &mut Done, note: &mut Msg) {
    ui.horizontal_wrapped(|ui| {
        let took = crate::view::took(r.took);
        ui.label(match r.mids.compounds.len() {
            1 => tf!(
                "1 substance with MIDs in {took}",
                "1 Substanz mit MIDs in {took}"
            ),
            n => tf!(
                "{n} substances with MIDs in {took}",
                "{n} Substanzen mit MIDs in {took}"
            ),
        });
        let mids = [
            t("Export MIDs…", "MIDs exportieren…"),
            t("Copy MIDs", "MIDs kopieren"),
        ];
        let copied = msg!("Copied the MIDs.", "MIDs kopiert.");
        export_buttons(ui, mids, copied, "mids.txt", &r.text, note);
        if let Some(text) = &r.stats_text {
            let stats = [
                t("Export statistics…", "Statistik exportieren…"),
                t("Copy statistics", "Statistik kopieren"),
            ];
            let copied = msg!("Copied the statistics.", "Statistik kopiert.");
            export_buttons(ui, stats, copied, "mid-statistics.txt", text, note);
        }
    });
    let error = ui.visuals().error_fg_color;
    if r.mids.compounds.is_empty() {
        ui.colored_label(
            error,
            t(
                "No compounds found during MID determination. Check library and RI window.",
                "Bei der MID-Bestimmung keine Substanzen gefunden. Bibliothek und RI-Fenster prüfen.",
            ),
        );
    }
    for e in &r.errors {
        ui.colored_label(error, &**e);
    }
    for e in &r.mids.errors {
        ui.colored_label(error, errors::text(e));
    }
    let mut mv = None;
    let n = r.mids.compounds.len();
    egui::ScrollArea::both().show(ui, |ui| {
        for (i, c) in r.mids.compounds.iter().enumerate() {
            let id = ui.make_persistent_id(("mid compound", c.compound));
            let mut state = egui::collapsing_header::CollapsingState::load_with_default_open(
                ui.ctx(),
                id,
                i == 0,
            );
            let head = ui.horizontal(|ui| {
                state.show_toggle_button(ui, egui::collapsing_header::paint_default_icon);
                let name = egui::RichText::new(&r.lib[c.compound].name).strong();
                if ui
                    .add(egui::Label::new(name).sense(egui::Sense::click()))
                    .clicked()
                {
                    state.toggle(ui);
                }
                if n < 2 {
                    return;
                }
                let up = ui.add_enabled(i > 0, egui::Button::new("▲").small());
                if up
                    .on_hover_text(t("Move substance up", "Substanz nach oben"))
                    .clicked()
                {
                    mv = Some(Move::Compound(i, i - 1));
                }
                let down = ui.add_enabled(i + 1 < n, egui::Button::new("▼").small());
                if down
                    .on_hover_text(t("Move substance down", "Substanz nach unten"))
                    .clicked()
                {
                    mv = Some(Move::Compound(i, i + 1));
                }
            });
            state.show_body_indented(&head.response, ui, |ui| {
                for (j, f) in c.fragments.iter().enumerate() {
                    let stats = r.stats.get(i).map_or(&[][..], |s| &s[j]);
                    ui.strong(tf!("Fragment {}–{}", "Fragment {}–{}", f.lo, f.hi));
                    let id = (c.compound, j);
                    mv = fragment_grid(ui, f, (&r.columns, &r.paths), stats, id).or(mv);
                    fragment_plot(ui, f, (&r.columns, &r.paths), stats, id);
                }
            });
        }
    });
    if let Some(m) = mv {
        r.apply(m);
    }
}

/// "Export …" and "Copy …" of one of the old program's text exports; `labels` are `[export button, copy button]`,
/// `copied` the note after copying.
fn export_buttons(
    ui: &mut egui::Ui,
    labels: [&str; 2],
    copied: Msg,
    file: &str,
    text: &str,
    note: &mut Msg,
) {
    let [export, copy] = labels;
    if ui.button(export).clicked()
        && let Some(p) = rfd::FileDialog::new()
            .add_filter(t("Text", "Text"), &["txt", "tsv"])
            .set_file_name(file)
            .save_file()
    {
        *note = match std::fs::write(&p, text) {
            Ok(()) => msg!("Saved {}", "Gespeichert: {}", p.display()),
            Err(e) => msg!("Export failed: {e}", "Export fehlgeschlagen: {e}"),
        };
    }
    if ui
        .button(copy)
        .on_hover_text(t(
            "Tab-separated, as the export file",
            "Tabulatorgetrennt, wie die Exportdatei",
        ))
        .clicked()
    {
        ui.ctx().copy_text(text.to_owned());
        *note = copied;
    }
}

/// Cell text of a value, a dash for none (NaN in the old program).
fn cell(v: Option<f64>) -> String {
    v.filter(|v| v.is_finite())
        .map_or("–".into(), |v| format!("{v:.4}"))
}

/// Header of run column `i`: its name with the path on hover, between ‹ › that move it.
fn run_head(ui: &mut egui::Ui, (columns, paths): (&[String], &[String]), i: usize) -> Option<Move> {
    let mut mv = None;
    ui.horizontal(|ui| {
        if i > 0
            && ui
                .small_button("‹")
                .on_hover_text(t("Move run left", "Lauf nach links"))
                .clicked()
        {
            mv = Some(Move::Run(i, i - 1));
        }
        let head = ui.strong(shown_name(columns, paths, i));
        if let Some(path) = paths.get(i) {
            head.on_hover_text(path);
        }
        if i + 1 < columns.len()
            && ui
                .small_button("›")
                .on_hover_text(t("Move run right", "Lauf nach rechts"))
                .clicked()
        {
            mv = Some(Move::Run(i, i + 1));
        }
    });
    mv
}

/// The MID table of a fragment; ‹ › in a run's header move that run in every table and export.
fn fragment_grid(
    ui: &mut egui::Ui,
    f: &FragmentMid,
    (columns, paths): (&[String], &[String]),
    stats: &[(f64, f64)],
    id: (usize, usize),
) -> Option<Move> {
    let mut mv = None;
    egui::Grid::new(id).striped(true).show(ui, |ui| {
        ui.label("");
        for i in 0..columns.len() {
            mv = run_head(ui, (columns, paths), i).or(mv);
        }
        if !stats.is_empty() {
            ui.strong(t("Mean", "Mittelwert")).on_hover_text(t(
                "Mean of the labelled runs; shares of 0 and runs without MID do not count",
                "Mittelwert der markierten Läufe; Anteile von 0 und Läufe ohne MID zählen nicht",
            ));
            ui.strong(t("StdErr", "Std.-Fehler")).on_hover_text(t(
                "Standard error of the mean (SD / √n), none for fewer than two runs",
                "Standardfehler des Mittelwerts (SD / √n), keiner bei weniger als zwei Läufen",
            ));
        }
        ui.end_row();
        for k in 0..f.isotopomers {
            ui.label(format!("M{k} (m/z {})", f.lo as usize + k));
            for fit in &f.fits {
                let fit = fit.as_ref();
                let v = fit.and_then(|x| x.mids.get(k).copied());
                let m = cell(v);
                let shown = ui.label(match fit.and_then(|x| x.cis.get(k).copied().flatten()) {
                    Some(ci) => format!("{m} ± {ci:.4}"),
                    None => m,
                });
                if v.is_some_and(|v| v < 0.0) {
                    shown.on_hover_text(t(
                        "Below 0 after the correction: noise around 0, not a real share",
                        "Nach der Korrektur unter 0: Rauschen um 0, kein echter Anteil",
                    ));
                }
            }
            if let Some((mean, se)) = stats.get(k) {
                ui.label(cell(Some(*mean)));
                ui.label(cell(Some(*se)));
            }
            ui.end_row();
        }
        ui.label("R²").on_hover_text(t(
            "How well the isotopomers explain the measured spectra; 1 is a perfect fit. It is always 1 \
             when there are no more measured values than isotopomers, as with one library spectrum. \
             – means no fit: uncorrected MIDs are not fitted.",
            "Wie gut die Isotopomere die gemessenen Spektren erklären; 1 ist eine perfekte Anpassung. \
             Es ist immer 1, wenn es nicht mehr Messwerte als Isotopomere gibt, wie bei einem \
             Bibliotheksspektrum. – heißt keine Anpassung: unkorrigierte MIDs werden nicht angepasst.",
        ));
        for fit in &f.fits {
            ui.label(cell(fit.as_ref().and_then(|x| x.r2)));
        }
        ui.end_row();
    });
    mv
}

/// Key, name and per isotopomer the value and its error.
type Series<'a> = (&'a str, &'a str, Vec<(f64, Option<f64>)>);

/// Per series: key for the plot (hidden series stay hidden in another language), name and per
/// isotopomer the value and its error. Runs are keyed by their path, since runs like `s.1.cdf` and
/// `s.2.cdf` share the column name `s`; the group methods have no paths and one column.
fn plot_series<'a>(
    f: &FragmentMid,
    columns: &'a [String],
    paths: &'a [String],
    stats: &[(f64, f64)],
) -> Vec<Series<'a>> {
    let mut series: Vec<_> = (columns.iter().zip(&f.fits).enumerate())
        .map(|(j, (c, fit))| {
            let name = shown_name(columns, paths, j);
            let v = (0..f.isotopomers).map(|k| match fit {
                Some(x) => (
                    x.mids.get(k).copied().unwrap_or(f64::NAN),
                    x.cis.get(k).copied().flatten(),
                ),
                None => (f64::NAN, None),
            });
            (paths.get(j).unwrap_or(c).as_str(), name, v.collect())
        })
        .collect();
    if !stats.is_empty() {
        series.push((
            "Mean ± StdErr",
            t("Mean ± StdErr", "Mittelwert ± Std.-Fehler"),
            stats.iter().map(|s| (s.0, Some(s.1))).collect(),
        ));
    }
    series
}

/// Bars per isotopomer: one per run or group, and for the per-run methods the mean; whiskers show
/// the CI of a group or the standard error of the mean.
fn fragment_plot(
    ui: &mut egui::Ui,
    f: &FragmentMid,
    (columns, paths): (&[String], &[String]),
    stats: &[(f64, f64)],
    id: (usize, usize),
) {
    let n = f.isotopomers;
    let series = plot_series(f, columns, paths, stats);
    let width = 0.8 / series.len().max(1) as f64;
    let whisker = ui.visuals().text_color();
    let top = ui.cursor().min;
    Plot::new(("mid plot", id))
        .height(160.0)
        // By id: the default groups by name, which would merge runs with the same short name again.
        .legend(Legend::default().grouping(LegendGrouping::ById))
        .allow_drag(false)
        .allow_zoom(false)
        .allow_scroll(false)
        .y_axis_label(t("Share", "Anteil"))
        .x_axis_formatter(move |m, _| match md_core::view::axis_index(m.value) {
            Some(k) if k < n => format!("M{k}"),
            _ => String::new(),
        })
        .show(ui, |p| {
            for (j, (key, name, values)) in series.iter().enumerate() {
                let x = |k: usize| k as f64 - 0.4 + width * (j as f64 + 0.5);
                let bars = (values.iter().enumerate())
                    .filter(|(_, v)| v.0.is_finite())
                    .map(|(k, v)| Bar::new(x(k), v.0).width(width * 0.9).name(format!("M{k}")));
                p.bar_chart(
                    BarChart::new(*name, bars.collect())
                        .id(*key)
                        .color(crate::color(j)),
                );
                // Same name and id as the bars: one legend entry hides both.
                let whiskers = (values.iter().enumerate())
                    .filter_map(|(k, (v, e))| {
                        Some((k, *v, e.filter(|e| e.is_finite() && v.is_finite())?))
                    })
                    .map(|(k, v, e)| {
                        BoxElem::new(x(k), BoxSpread::new(v - e, v, v, v, v + e))
                            .box_width(0.0)
                            .whisker_width(width * 0.5)
                            .stroke(egui::Stroke::new(1.0, whisker))
                    });
                p.box_plot(
                    BoxPlot::new(*name, whiskers.collect())
                        .id(*key)
                        .allow_hover(false),
                );
            }
        });
    crate::plot_export::button(
        ui,
        crate::plot_export::Area::since(ui, top),
        &format!("mid-{}-{}", id.0 + 1, id.1 + 1),
    );
}

#[cfg(test)]
impl MidWindow {
    /// Labelled runs and the sum formula method, for the screenshots.
    pub fn demo(&mut self, labelled: Vec<PathBuf>) {
        (self.labelled, self.kind) = (labelled, Kind::Formula);
    }

    pub fn busy(&self) -> bool {
        self.job.is_some()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_core::mid::{CompoundMid, MidFit};

    /// Two substances, runs `a`/`b` with the shares `x` and `1 - x` per run.
    fn done(method: Method, columns: &[&str]) -> Done {
        let fit = |x: f64| {
            Some(MidFit {
                mids: vec![x, 1.0 - x],
                r2: Some(x),
                cis: vec![Some(x / 10.0), None],
            })
        };
        let compound = |c: usize| CompoundMid {
            compound: c,
            fragments: vec![FragmentMid {
                lo: 73,
                hi: 74,
                isotopomers: 2,
                fits: (0..columns.len())
                    .map(|j| fit(0.1 * (c * 3 + j + 1) as f64))
                    .collect(),
            }],
        };
        let mut d = Done {
            method,
            lib: ["Alanine", "Glycine"].map(LibraryCompound::new).to_vec(),
            columns: columns.iter().map(|c| c.to_string()).collect(),
            paths: Vec::new(),
            runs: Vec::new(),
            mids: Mids {
                compounds: vec![compound(0), compound(1)],
                errors: Vec::new(),
            },
            stats: Vec::new(),
            text: String::new(),
            stats_text: None,
            errors: Vec::new(),
            scanned: 0,
            took: Duration::ZERO,
        };
        d.refresh();
        d
    }

    #[test]
    fn exports_follow_the_moved_runs_and_substances() {
        let mut d = done(Method::Formula, &["a", "b"]);
        let before = (d.text.clone(), d.stats_text.clone());
        d.apply(Move::Run(1, 0));
        assert_eq!(d.columns, ["b", "a"]);
        let head = "\tb\ta\t\tCI (b)\tCI (a)\t\tR2 (b)\tR2 (a)\n";
        assert!(d.text.starts_with(head), "{}", d.text);
        // Run b of Alanine is 0.2, run a 0.1: the values move with their column.
        assert!(
            d.text
                .contains("73 (M0)\t0.2\t0.1\t\t0.02\t0.01\t\t0.2\t0.1\n")
        );
        d.apply(Move::Compound(0, 1));
        assert!(d.text.starts_with("\tb\ta"));
        assert!(d.text.find("Glycine").unwrap() < d.text.find("Alanine").unwrap());
        let stats = d.stats_text.as_deref().unwrap();
        assert!(stats.find("Glycine").unwrap() < stats.find("Alanine").unwrap());
        // The table's mean of Glycine M0 (0.4 and 0.5) moved up with it.
        assert!((d.stats[0][0][0].0 - 0.45).abs() < 1e-12);
        d.apply(Move::Compound(1, 0));
        d.apply(Move::Run(0, 1));
        assert_eq!((d.text, d.stats_text), before);
    }

    #[test]
    fn moving_runs_keeps_the_statistics_bit_for_bit() {
        // Alanine M0 is 0.1, 0.2, 0.3: the mean is 0.20000000000000004, summed backwards 0.19999999999999998.
        let mut d = done(Method::Formula, &["a", "b", "c"]);
        let (stats, stats_text) = (d.stats.clone(), d.stats_text.clone());
        d.apply(Move::Run(2, 0));
        d.apply(Move::Compound(0, 1));
        d.apply(Move::Compound(1, 0));
        assert_eq!(d.columns, ["c", "b", "a"]);
        assert_eq!((d.stats, d.stats_text), (stats, stats_text));
    }

    #[test]
    fn runs_with_the_same_short_name_get_distinct_plot_keys() {
        let mut d = done(Method::Formula, &["s", "s"]);
        d.paths = ["x/s.1.cdf", "x/s.2.cdf"].map(String::from).to_vec();
        d.apply(Move::Run(1, 0));
        let f = &d.mids.compounds[0].fragments[0];
        let series = plot_series(f, &d.columns, &d.paths, &[(0.0, 0.0)]);
        let keys: Vec<_> = series.iter().map(|s| s.0).collect();
        assert_eq!(keys, ["x/s.2.cdf", "x/s.1.cdf", "Mean ± StdErr"]);
        // The key moved with its run's values (run 2 is 0.2 for Alanine).
        assert_eq!(series[0].2[0].0, 0.2);
    }

    /// Short line segments in the last frame with the whiskers' colour; the grid lines are as tall as the plot.
    fn whisker_lines(h: &egui_kittest::Harness<'_, Done>) -> usize {
        fn count(s: &egui::Shape, color: egui::Color32) -> usize {
            match s {
                egui::Shape::Vec(v) => v.iter().map(|s| count(s, color)).sum(),
                egui::Shape::LineSegment { stroke, points } => {
                    (stroke.color == color && (points[0] - points[1]).length() < 10.0) as usize
                }
                _ => 0,
            }
        }
        let color = egui::Visuals::dark().text_color();
        (h.output().shapes.iter())
            .map(|c| count(&c.shape, color))
            .sum()
    }

    #[test]
    fn hiding_a_run_hides_its_bars_and_whiskers_only() {
        use egui_kittest::kittest::{NodeT, Queryable};
        let mut d = done(Method::Formula, &["s", "s"]);
        d.paths = ["x/s.1.cdf", "x/s.2.cdf"].map(String::from).to_vec();
        let mut h = egui_kittest::Harness::new_ui_state(
            |ui, d: &mut Done| {
                let f = &d.mids.compounds[0].fragments[0];
                fragment_plot(ui, f, (&d.columns, &d.paths), &[], (0, 0));
            },
            d,
        );
        h.run();
        let toggled = |h: &egui_kittest::Harness<'_, Done>| {
            let boxes = h.get_all_by_role_and_label(egui::accesskit::Role::CheckBox, "s");
            boxes
                .map(|n| n.accesskit_node().toggled())
                .collect::<Vec<_>>()
        };
        // Two legend entries although both runs are called `s`: a legend by name would merge them.
        let on = Some(egui::accesskit::Toggled::True);
        let off = Some(egui::accesskit::Toggled::False);
        assert_eq!(toggled(&h), [on, on]);
        let both = whisker_lines(&h);
        let click = |h: &mut egui_kittest::Harness<'_, Done>, i: usize| {
            let node = (h.get_all_by_role_and_label(egui::accesskit::Role::CheckBox, "s")).nth(i);
            node.unwrap().click();
            h.run();
        };
        click(&mut h, 1);
        assert_eq!(toggled(&h), [on, off]);
        let one = whisker_lines(&h);
        click(&mut h, 0);
        assert_eq!(toggled(&h), [off, off]);
        // Each run has the same number of whisker lines.
        let none = whisker_lines(&h);
        assert!(one > 0, "no whiskers drawn");
        assert_eq!((both, none), (2 * one, 0));
    }

    #[test]
    fn group_column_follows_the_language_but_not_key_or_export() {
        let d = done(Method::Uncorrected, &[GROUP_COLUMN]);
        let f = &d.mids.compounds[0].fragments[0];
        let series = |l| {
            crate::lang::set(l);
            let s = plot_series(f, &d.columns, &d.paths, &[]);
            (s[0].0, s[0].1, shown_name(&d.columns, &d.paths, 0))
        };
        assert_eq!(
            series(crate::lang::Lang::En),
            ("labelled", "labelled", "labelled")
        );
        assert_eq!(
            series(crate::lang::Lang::De),
            ("labelled", "markiert", "markiert")
        );
        // Rebuilt while German is on: the export must stay English.
        let german = d.mids_text();
        crate::lang::set(crate::lang::Lang::En);
        assert!(!german.contains("markiert"));
        assert_eq!(german, d.text);
    }

    #[test]
    fn group_export_follows_the_moved_substances() {
        let mut d = done(Method::Uncorrected, &[GROUP_COLUMN]);
        d.apply(Move::Compound(1, 0));
        assert!(d.text.starts_with("Glycine_73 (M0)\t0.4\t"));
        assert!(d.stats.is_empty() && d.stats_text.is_none());
    }
}
