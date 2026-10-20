//! MetaboliteDetector desktop app (egui). Analysis lives in `md-core`; this crate only draws.

use std::path::PathBuf;
use std::sync::{Arc, mpsc};
use std::time::Instant;

use eframe::egui;
use md_core::Chromatogram;
use md_core::deconv::Component;

mod batch;
mod calibration;
mod central;
mod export;
mod identify;
mod isotopes;
mod jobs;
mod library_editor;
mod loading;
mod menu;
mod results;
mod run_list;
mod settings;
mod statistics;
mod status;
mod style;
#[cfg(test)]
mod ui_tests;
mod view;

use jobs::{poll_job, spawn};

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("MetaboliteDetector")
            .with_inner_size([1200.0, 800.0]),
        ..Default::default()
    };
    eframe::run_native(
        "MetaboliteDetector",
        options,
        Box::new(|cc| {
            style::apply(&cc.egui_ctx);
            let mut app = App::default();
            for path in std::env::args_os().skip(1).map(PathBuf::from) {
                if identify::is_library(&path) {
                    app.open_library(&cc.egui_ctx, path);
                } else {
                    app.open(&cc.egui_ctx, path);
                }
            }
            Ok(Box::new(app))
        }),
    )
}

#[derive(Clone, Copy, PartialEq)]
enum Trace {
    Tic,
    Sic,
    Both,
}

struct App {
    runs: Vec<Run>,
    /// Index into `runs` whose scan cursor and spectrum are shown.
    active: usize,
    loading: Vec<(PathBuf, mpsc::Receiver<Result<Loaded, String>>)>,
    load_started: Option<Instant>,
    load_note: String,
    error: Option<String>,
    trace: Trace,
    mz: f32,
    tolerance: f32,
    ri: calibration::Calibration,
    settings: settings::Settings,
    isotopes: isotopes::Isotopes,
    library: Option<Arc<identify::Library>>,
    library_job: Option<identify::LibraryJob>,
    /// Search settings whose library preparation crashed; not retried until they change.
    library_failed: Option<md_core::search::SearchSettings>,
    editor: library_editor::LibraryEditor,
    batch: batch::Batch,
    statistics: statistics::Statistics,
    results: results::Results,
    export: export::Export,
    view: view::View,
}

impl Default for App {
    fn default() -> Self {
        Self {
            runs: Vec::new(),
            active: 0,
            loading: Vec::new(),
            load_started: None,
            load_note: String::new(),
            error: None,
            trace: Trace::Tic,
            // TMS fragment; many runs start scanning above m/z 73.
            mz: 147.0,
            tolerance: 0.5,
            ri: calibration::Calibration::default(),
            settings: settings::Settings::default(),
            isotopes: isotopes::Isotopes::default(),
            library: None,
            library_job: None,
            library_failed: None,
            editor: library_editor::LibraryEditor::default(),
            batch: batch::Batch::default(),
            statistics: statistics::Statistics::default(),
            results: results::Results::default(),
            export: export::Export::default(),
            view: view::View::default(),
        }
    }
}

/// (m/z, tolerance) of a SIC.
type SicKey = (f32, f32);

/// A loaded run with warnings about its side files (`.cmp`, `.calib`).
type Loaded = (Run, Vec<String>);

struct Run {
    path: PathBuf,
    name: String,
    data: Arc<Chromatogram>,
    /// Scan times in minutes and TIC, cached once per file.
    minutes: Vec<f64>,
    tic: Vec<f64>,
    /// SIC for the (m/z, tolerance) in `sic_for`; recomputed on a worker thread when that changes.
    sic: Vec<f64>,
    sic_for: Option<SicKey>,
    sic_job: Option<(SicKey, mpsc::Receiver<Vec<f64>>)>,
    /// Deconvoluted substances with the applied settings, ordered by retention time.
    components: Vec<Component>,
    /// Ion mass of each peak the components were fitted from (`FittedPeak::peak`).
    peak_masses: Vec<u32>,
    /// Calibration saved with the run by the old program, handed to the calibration window on load.
    saved_calibration: Option<md_core::ri::RiCalibration>,
    analysis: settings::Analysis,
    scan: usize,
    /// Selected component; its deconvoluted spectrum replaces the raw one.
    component: Option<usize>,
    visible: bool,
    /// Library hits of `components` once a library is open.
    ident: Option<identify::Identified>,
    ident_job: Option<identify::IdentJob>,
    /// Counts every replacement of `components` and `ident`, so caches can tell a new result from the old.
    rev: u64,
}

/// Same colour per run in every plot and in the run list (egui_plot's golden-ratio hues).
fn color(i: usize) -> egui::Color32 {
    egui::ecolor::Hsva::new(i as f32 * 0.618_034, 0.85, 0.5, 1.0).into()
}

impl eframe::App for App {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.poll_jobs(ui.ctx());
        self.menu_bar(ui);
        self.windows(ui.ctx());
        self.status_bar(ui);
        self.run_list(ui);
        if let Some(run) = self.runs.get(self.active) {
            identify::hits_panel(ui, run);
        }
        self.central_panel(ui);
    }
}

impl App {
    /// Collects finished background work before drawing.
    fn poll_jobs(&mut self, ctx: &egui::Context) {
        self.poll_loading();
        if let Some(e) = self.settings.update_runs(ctx, &mut self.runs) {
            self.error = Some(e);
        }
        self.update_identification(ctx);
        self.poll_export();
    }

    /// The floating windows opened from the Analysis menu.
    fn windows(&mut self, ctx: &egui::Context) {
        self.ri
            .window(ctx, &self.runs, self.active, self.library.as_deref());
        self.settings.window(ctx);
        self.batch.analysis = self.settings.key();
        // Not while the library is rebuilt after a switch of RT-only: its anchors would not match.
        let cal = (self.ri.cal.as_deref()).filter(|_| {
            (self.library.as_ref()).is_some_and(|l| l.prepared.is_rt_only() == self.ri.rt_only)
        });
        self.batch.window(ctx, self.library.as_ref(), cal);
        self.statistics
            .window(ctx, self.batch.shown_table(), self.batch.z_scored());
        self.isotopes.window(ctx);
        self.results.window(ctx, &mut self.runs, self.active);
        self.info_window(ctx);
        self.library_editor(ctx);
    }
}
