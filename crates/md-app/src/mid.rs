//! MID window (F16, old `MIDWizard`): scans labelled (and unlabelled) runs on a worker thread
//! against the library's substances with fragments and shows their mass isotopomer distributions
//! (`md_core::mid`).

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::{Duration, Instant};

use crate::layout::ToolWindow;
use eframe::egui;
use egui_plot::{Bar, BarChart, BoxElem, BoxPlot, BoxSpread, Legend, Plot};
use md_core::batch::{QuantSettings, ShiftReference, complete_scan, scan_runs};
use md_core::isotope::IsotopeTable;
use md_core::library::LibraryCompound;
use md_core::mid::{FragmentMid, Method, Mids, determine, group_stat, mid_library};
use md_core::ri::RiCalibration;
use md_core::search::PreparedLibrary;
use md_core::settings::GcmsSettings;

use crate::identify::Library;
use crate::spawn;

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
    mids: Mids,
    /// Per-run methods: mean and standard error of the labelled runs per substance, fragment and
    /// isotopomer (old group `GL`); empty for the group methods.
    stats: Vec<Vec<Vec<(f64, f64)>>>,
    /// The old program's text exports: the MIDs, and for the per-run methods the statistics.
    text: String,
    stats_text: Option<String>,
    errors: Vec<String>,
    took: Duration,
}

/// A move in the result view: run columns (old "Change column order…") or substances with their
/// diagrams (old "Change order of diagrams…"), from one position to another.
#[derive(Clone, Copy)]
enum Move {
    Run(usize, usize),
    Compound(usize, usize),
}

impl Done {
    /// Applies a move to the table, the diagrams and both exports; the values stay as determined.
    fn apply(&mut self, m: Move) {
        match m {
            Move::Run(a, b) => {
                self.columns.swap(a, b);
                if !self.paths.is_empty() {
                    self.paths.swap(a, b);
                }
                (self.mids.compounds.iter_mut())
                    .flat_map(|c| &mut c.fragments)
                    .for_each(|f| f.fits.swap(a, b));
            }
            Move::Compound(a, b) => self.mids.compounds.swap(a, b),
        }
        self.refresh();
    }

    /// Statistics and text exports of the MIDs in their current order.
    fn refresh(&mut self) {
        let per_run = self.method.per_run();
        self.stats = (self.mids.compounds.iter())
            .filter(|_| per_run)
            .map(|c| {
                (c.fragments.iter())
                    .map(|f| (0..f.isotopomers).map(|k| group_stat(f, k)).collect())
                    .collect()
            })
            .collect();
        self.text = match per_run {
            true => (self.mids).run_text(&self.lib, self.method, &self.columns),
            false => self.mids.group_text(&self.lib),
        };
        self.stats_text = per_run.then(|| self.mids.stats_text(&self.lib));
    }
}

struct Job {
    total: usize,
    done: Arc<AtomicUsize>,
    cancel: Arc<AtomicBool>,
    started: Instant,
    rx: mpsc::Receiver<Result<Done, String>>,
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
    job: Option<Job>,
    result: Option<Done>,
    note: String,
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
            note: String::new(),
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
        egui::Window::new("Mass isotopomer distributions")
            .tool_window(ctx, 3, [800.0, 600.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui, lib, cal, reference));
        self.open = open;
    }

    fn poll(&mut self, ctx: &egui::Context) {
        let Some(job) = &self.job else { return };
        match job.rx.try_recv() {
            Ok(Ok(done)) => (self.result, self.note) = (Some(done), String::new()),
            Ok(Err(e)) => self.note = e,
            Err(mpsc::TryRecvError::Empty) => {
                ctx.request_repaint_after(Duration::from_millis(200));
                return;
            }
            Err(mpsc::TryRecvError::Disconnected) => {
                self.note = "MID determination crashed.".into()
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
            false => vec!["labelled".into()],
        };
        let column_paths: Vec<String> = match method.per_run() {
            true => (self.labelled.iter())
                .map(|p| p.display().to_string())
                .collect(),
            false => Vec::new(),
        };
        let (nl, nu) = (self.labelled.len(), unlabelled.len());
        let rescan = self.rescan;
        let (done, cancel) = (
            Arc::new(AtomicUsize::new(0)),
            Arc::new(AtomicBool::new(false)),
        );
        let (l, d, c, tracer, table) = (
            lib.clone(),
            done.clone(),
            cancel.clone(),
            self.tracer.trim().to_owned(),
            self.table.clone(),
        );
        let rx = spawn(ctx, move || {
            let started = Instant::now();
            let mids_lib = mid_library(&l.compounds);
            if mids_lib.is_empty() {
                return Err("No substance of the library has a fragment with sum formula (MID fragments in the library editor).".into());
            }
            // The shift reference is an index into the scanned library, which keeps only
            // substances with fragments.
            s.shift = (shift_ref.as_ref())
                .and_then(|(n, cutoff)| ShiftReference::find(&mids_lib, n, *cutoff));
            let prepared = match l.prepared.is_rt_only() {
                true => PreparedLibrary::new_rt_only(&mids_lib, &l.settings),
                false => PreparedLibrary::new(&mids_lib, &l.settings),
            };
            let cancelled = || "Cancelled.".to_owned();
            let (mut runs, mut errors) =
                scan_runs(&paths, &prepared, &s, 4, &d, &c).ok_or_else(cancelled)?;
            // The quant ions are the fragment ions, so the rescan fills those that were not deconvoluted.
            let more = complete_scan(
                &mut runs,
                &mids_lib,
                &s,
                &QuantSettings::default(),
                rescan,
                4,
                &c,
            )
            .ok_or_else(cancelled)?;
            errors.extend(more.errors);
            if let Some((n, _)) = shift_ref.filter(|_| s.shift.is_none()) {
                errors.push(format!(
                    "RT shift reference {n} is not a substance with MID fragments; scanned without RT shift."
                ));
            }
            // shortcut: a run that fails to load is dropped by scan_runs, which shifts the group
            // indices; the runs are matched back by path.
            let index = |p: &PathBuf| runs.iter().position(|r| &r.path == p);
            let labelled: Vec<usize> = paths[..nl].iter().filter_map(index).collect();
            let unlabelled: Vec<usize> = paths[nl..nl + nu].iter().filter_map(index).collect();
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
                    f.fits = (paths[..nl].iter())
                        .map(|p| index(p).and_then(|_| fits.next().flatten()))
                        .collect();
                }
            }
            let mut done = Done {
                method,
                lib: mids_lib,
                columns,
                paths: column_paths,
                mids,
                stats: Vec::new(),
                text: String::new(),
                stats_text: None,
                errors,
                took: started.elapsed(),
            };
            done.refresh();
            Ok(done)
        });
        self.job = Some(Job {
            total: nl + nu,
            done,
            cancel,
            started: Instant::now(),
            rx,
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
        ui.horizontal(|ui| match &self.job {
            Some(job) => {
                let d = job.done.load(Relaxed);
                ui.add(
                    egui::ProgressBar::new(d as f32 / job.total.max(1) as f32)
                        .text(format!(
                            "{d} / {} runs, {:.0?}",
                            job.total,
                            job.started.elapsed()
                        ))
                        .desired_width(300.0),
                );
                if ui.button("Cancel").clicked() {
                    job.cancel.store(true, Relaxed);
                }
            }
            None => {
                let problem = self.problem();
                // Not while the library is rebuilt for new search settings, as in the batch window.
                let ready =
                    lib.filter(|l| problem.is_none() && l.settings == self.applied.search());
                if ui
                    .add_enabled(ready.is_some(), egui::Button::new("Start"))
                    .clicked()
                    && let Some(l) = ready
                {
                    self.start(ui.ctx(), l, cal, reference);
                }
                match (lib, problem) {
                    (None, _) => ui.label("Open a library first (File → Open library)."),
                    (_, Some(p)) => ui.label(p),
                    (Some(l), None) if l.settings != self.applied.search() => {
                        ui.label("Library is being prepared for the new search settings.")
                    }
                    (Some(l), None) => ui.label(format!(
                        "Library {}, {}",
                        l.name,
                        if cal.is_some() {
                            "RI calibrated"
                        } else {
                            "no RI calibration"
                        }
                    )),
                };
            }
        });
        if !self.note.is_empty() {
            ui.label(&self.note);
        }
        let Some(r) = &mut self.result else { return };
        ui.separator();
        result_view(ui, r, &mut self.note);
    }

    /// Why the MIDs cannot be determined yet.
    fn problem(&self) -> Option<&'static str> {
        let tracer = self.tracer.trim();
        if self.labelled.is_empty() {
            Some("Add labelled runs.")
        } else if self.kind == Kind::Reference && self.unlabelled.is_empty() {
            Some("Add unlabelled runs for the reference correction.")
        } else if !tracer.is_empty() && self.table.isotopes(tracer).is_none() {
            Some("The tracer is not an element of the isotope table.")
        } else if tracer.is_empty() && matches!(self.kind, Kind::Library | Kind::Formula) {
            Some("This correction needs a tracer element.")
        } else {
            None
        }
    }

    fn settings(&mut self, ui: &mut egui::Ui) {
        egui::Grid::new("mid settings").show(ui, |ui| {
            ui.label("Correction");
            ui.vertical(|ui| {
                // A grid column starts as narrow as its first frame; the short titles never wrap.
                ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
                for (kind, title, hint) in [
                    (Kind::Uncorrected, "Uncorrected MIDs", "Mean and confidence interval of the labelled runs' normalised fragment spectra."),
                    (Kind::Reference, "Corrected using reference", "Labelled runs against unlabelled runs of the same substances."),
                    (Kind::Library, "Corrected using library", "Every labelled run against the library spectrum, corrected for the tracer's natural M+1."),
                    (Kind::Formula, "Corrected using sum formulas", "Every labelled run against the natural isotope pattern of the fragment's sum formula."),
                ] {
                    ui.radio_value(&mut self.kind, kind, title).on_hover_text(hint);
                }
            });
            ui.end_row();
            let label = ui.label("Tracer element");
            ui.add(egui::TextEdit::singleline(&mut self.tracer).desired_width(40.0))
                .labelled_by(label.id);
            ui.end_row();
            if self.kind == Kind::Reference {
                // Pre-fill the tracer's natural M+1 abundance whenever the tracer changes.
                if self.abundance_for != self.tracer {
                    self.abundance_for = self.tracer.clone();
                    if let Some(i) = self.table.isotopes(self.tracer.trim()).and_then(|i| i.get(1)) {
                        self.abundance = i.abundance * 100.0;
                    }
                }
                ui.label("M1 correction");
                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.m1, "natural abundance");
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
            let hint = "Search the ion traces again for fragment ions the deconvolution missed";
            let rescan = ui.label("Rescan missing ions").on_hover_text(hint);
            ui.add(egui::Checkbox::without_text(&mut self.rescan))
                .labelled_by(rescan.id)
                .on_hover_text(hint);
            ui.end_row();
        });
        runs(ui, "Labelled runs", &mut self.labelled);
        if self.kind == Kind::Reference {
            runs(ui, "Unlabelled runs", &mut self.unlabelled);
        }
    }
}

fn runs(ui: &mut egui::Ui, title: &str, paths: &mut Vec<PathBuf>) {
    ui.horizontal(|ui| {
        ui.label(format!("{title}: {}", paths.len()));
        if ui.button("Add…").clicked()
            && let Some(p) = rfd::FileDialog::new()
                .add_filter("NetCDF", &["cdf", "CDF", "nc"])
                .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
                .pick_files()
        {
            paths.extend(p);
            paths.sort();
            paths.dedup();
        }
        if ui.button("Clear").clicked() {
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
fn result_view(ui: &mut egui::Ui, r: &mut Done, note: &mut String) {
    ui.horizontal_wrapped(|ui| {
        ui.label(format!(
            "{} with MIDs in {}",
            match r.mids.compounds.len() {
                1 => "1 substance".into(),
                n => format!("{n} substances"),
            },
            crate::view::took(r.took)
        ));
        export_buttons(ui, "MIDs", "mids.txt", &r.text, note);
        if let Some(text) = &r.stats_text {
            export_buttons(ui, "statistics", "mid-statistics.txt", text, note);
        }
    });
    let error = ui.visuals().error_fg_color;
    if r.mids.compounds.is_empty() {
        ui.colored_label(
            error,
            "No compounds found during MID determination. Check library and RI window.",
        );
    }
    for e in r.errors.iter().chain(&r.mids.errors) {
        ui.colored_label(error, e);
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
                if up.on_hover_text("Move substance up").clicked() {
                    mv = Some(Move::Compound(i, i - 1));
                }
                let down = ui.add_enabled(i + 1 < n, egui::Button::new("▼").small());
                if down.on_hover_text("Move substance down").clicked() {
                    mv = Some(Move::Compound(i, i + 1));
                }
            });
            state.show_body_indented(&head.response, ui, |ui| {
                for (j, f) in c.fragments.iter().enumerate() {
                    let stats = r.stats.get(i).map_or(&[][..], |s| &s[j]);
                    ui.strong(format!("Fragment {}–{}", f.lo, f.hi));
                    let id = (c.compound, j);
                    mv = fragment_grid(ui, f, (&r.columns, &r.paths), stats, id).or(mv);
                    fragment_plot(ui, f, &r.columns, stats, id);
                }
            });
        }
    });
    if let Some(m) = mv {
        r.apply(m);
    }
}

/// "Export …" and "Copy …" of one of the old program's text exports.
fn export_buttons(ui: &mut egui::Ui, what: &str, file: &str, text: &str, note: &mut String) {
    if ui.button(format!("Export {what}…")).clicked()
        && let Some(p) = rfd::FileDialog::new()
            .add_filter("Text", &["txt", "tsv"])
            .set_file_name(file)
            .save_file()
    {
        *note = match std::fs::write(&p, text) {
            Ok(()) => format!("Saved {}", p.display()),
            Err(e) => format!("Export failed: {e}"),
        };
    }
    if ui
        .button(format!("Copy {what}"))
        .on_hover_text("Tab-separated, as the export file")
        .clicked()
    {
        ui.ctx().copy_text(text.to_owned());
        *note = format!("Copied the {what}.");
    }
}

/// Cell text of a value, a dash for none (NaN in the old program).
fn cell(v: Option<f64>) -> String {
    v.filter(|v| v.is_finite())
        .map_or("–".into(), |v| format!("{v:.4}"))
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
    let n = columns.len();
    egui::Grid::new(id).striped(true).show(ui, |ui| {
        ui.label("");
        for (i, c) in columns.iter().enumerate() {
            ui.horizontal(|ui| {
                if i > 0 && ui.small_button("‹").on_hover_text("Move run left").clicked() {
                    mv = Some(Move::Run(i, i - 1));
                }
                let head = ui.strong(c);
                if let Some(path) = paths.get(i) {
                    head.on_hover_text(path);
                }
                if i + 1 < n && ui.small_button("›").on_hover_text("Move run right").clicked() {
                    mv = Some(Move::Run(i, i + 1));
                }
            });
        }
        if !stats.is_empty() {
            ui.strong("Mean").on_hover_text(
                "Mean of the labelled runs; shares of 0 and runs without MID do not count",
            );
            ui.strong("StdErr").on_hover_text(
                "Standard error of the mean (SD / √n), none for fewer than two runs",
            );
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
                    shown.on_hover_text(
                        "Below 0 after the correction: noise around 0, not a real share",
                    );
                }
            }
            if let Some((mean, se)) = stats.get(k) {
                ui.label(cell(Some(*mean)));
                ui.label(cell(Some(*se)));
            }
            ui.end_row();
        }
        ui.label("R²").on_hover_text(
            "How well the isotopomers explain the measured spectra; 1 is a perfect fit. It is always 1 \
             when there are no more measured values than isotopomers, as with one library spectrum. \
             – means no fit: uncorrected MIDs are not fitted.",
        );
        for fit in &f.fits {
            ui.label(cell(fit.as_ref().and_then(|x| x.r2)));
        }
        ui.end_row();
    });
    mv
}

/// Bars per isotopomer: one per run or group, and for the per-run methods the mean; whiskers show
/// the CI of a group or the standard error of the mean.
fn fragment_plot(
    ui: &mut egui::Ui,
    f: &FragmentMid,
    columns: &[String],
    stats: &[(f64, f64)],
    id: (usize, usize),
) {
    let n = f.isotopomers;
    let fits = columns.iter().zip(&f.fits);
    // Per series: name and per isotopomer the value and its error.
    let mut series: Vec<(&str, Vec<_>)> = fits
        .map(|(c, fit)| {
            let v = (0..n).map(|k| match fit {
                Some(x) => (
                    x.mids.get(k).copied().unwrap_or(f64::NAN),
                    x.cis.get(k).copied().flatten(),
                ),
                None => (f64::NAN, None),
            });
            (c.as_str(), v.collect())
        })
        .collect();
    if !stats.is_empty() {
        series.push((
            "Mean ± StdErr",
            stats.iter().map(|s| (s.0, Some(s.1))).collect(),
        ));
    }
    let width = 0.8 / series.len().max(1) as f64;
    let whisker = ui.visuals().text_color();
    Plot::new(("mid plot", id))
        .height(160.0)
        .legend(Legend::default())
        .allow_drag(false)
        .allow_zoom(false)
        .allow_scroll(false)
        .y_axis_label("Share")
        .x_axis_formatter(move |m, _| {
            let k = m.value.round();
            match (m.value - k).abs() < 1e-6 && k >= 0.0 && (k as usize) < n {
                true => format!("M{k}"),
                false => String::new(),
            }
        })
        .show(ui, |p| {
            let mut errors = Vec::new();
            for (j, (name, values)) in series.iter().enumerate() {
                let x = |k: usize| k as f64 - 0.4 + width * (j as f64 + 0.5);
                let bars = (values.iter().enumerate())
                    .filter(|(_, v)| v.0.is_finite())
                    .map(|(k, v)| Bar::new(x(k), v.0).width(width * 0.9).name(format!("M{k}")));
                p.bar_chart(BarChart::new(*name, bars.collect()).color(crate::color(j)));
                for (k, (v, e)) in values.iter().enumerate() {
                    if let Some(e) = e.filter(|e| e.is_finite() && v.is_finite()) {
                        errors.push(
                            BoxElem::new(x(k), BoxSpread::new(v - e, *v, *v, *v, v + e))
                                .box_width(0.0)
                                .whisker_width(width * 0.5)
                                .stroke(egui::Stroke::new(1.0, whisker)),
                        );
                    }
                }
            }
            p.box_plot(BoxPlot::new("", errors).allow_hover(false));
        });
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
            mids: Mids {
                compounds: vec![compound(0), compound(1)],
                errors: Vec::new(),
            },
            stats: Vec::new(),
            text: String::new(),
            stats_text: None,
            errors: Vec::new(),
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
    fn group_export_follows_the_moved_substances() {
        let mut d = done(Method::Uncorrected, &["labelled"]);
        d.apply(Move::Compound(1, 0));
        assert!(d.text.starts_with("Glycine_73 (M0)\t0.4\t"));
        assert!(d.stats.is_empty() && d.stats_text.is_none());
    }
}
