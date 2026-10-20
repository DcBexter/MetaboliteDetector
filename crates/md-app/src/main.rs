//! MetaboliteDetector desktop app (egui). Analysis lives in `md-core`; this crate only draws.
// Release builds on Windows open no console window next to the app (shortcuts, double-click).
#![cfg_attr(not(debug_assertions), windows_subsystem = "windows")]

use std::path::PathBuf;
use std::sync::{Arc, mpsc};
use std::time::Instant;

use eframe::egui;
use md_core::Chromatogram;
use md_core::deconv::Component;
use md_core::peaks::Peak;

mod batch;
mod calibration;
mod central;
mod detail;
mod detected_peaks;
mod export;
mod export_all;
mod identify;
mod isotopes;
mod jobs;
mod layout;
mod library_editor;
mod loading;
mod menu;
mod mid;
mod persist;
mod plot_export;
mod results;
mod run_list;
mod settings;
mod shortcuts;
mod statistics;
mod status;
mod style;
mod tof_calibration;
#[cfg(test)]
mod ui_tests;
mod update;
mod vector;
mod view;

use jobs::{poll_job, spawn};

fn main() -> eframe::Result {
    // Without our config dir eframe would pick its own; name it, so `restore` checks that file too.
    let file = (persist::config_dir().or_else(|| eframe::storage_dir("MetaboliteDetector")))
        .map(|d| d.join("app.ron"));
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("MetaboliteDetector")
            .with_inner_size([1200.0, 800.0]),
        persistence_path: file.clone(),
        ..Default::default()
    };
    let result = eframe::run_native(
        "MetaboliteDetector",
        options,
        Box::new(move |cc| {
            style::apply(&cc.egui_ctx);
            let mut app = App::default();
            persist::restore(&mut app, &cc.egui_ctx, cc.storage, file.as_deref());
            app.updater.start(&cc.egui_ctx);
            for path in std::env::args_os().skip(1).map(PathBuf::from) {
                if identify::is_library(&path) {
                    app.open_library(&cc.egui_ctx, path);
                } else {
                    app.open(&cc.egui_ctx, path);
                }
            }
            Ok(Box::new(app))
        }),
    );
    // Without a console (release on Windows) a failed start (e.g. no OpenGL) would otherwise say nothing.
    if let Err(e) = &result {
        rfd::MessageDialog::new()
            .set_level(rfd::MessageLevel::Error)
            .set_title("MetaboliteDetector")
            .set_description(format!("The window could not be opened: {e}"))
            .show();
    }
    result
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
    tof: tof_calibration::TofWindow,
    library: Option<Arc<identify::Library>>,
    library_job: Option<identify::LibraryJob>,
    /// Search settings whose library preparation crashed; not retried until they change.
    library_failed: Option<md_core::search::SearchSettings>,
    nist: identify::Nist,
    editor: library_editor::LibraryEditor,
    batch: batch::Batch,
    statistics: statistics::Statistics,
    mid: mid::MidWindow,
    results: results::Results,
    export: export::Export,
    export_all: export_all::ExportAll,
    view: view::View,
    spectrum_mode: detail::Mode,
    updater: update::Updater,
    panels: layout::Panels,
    saved: persist::Saved,
    /// Detected Peaks also draws the modelled profile of each ion.
    model_profile: bool,
    about_open: bool,
    /// Run whose TIC zoom is shown; another run starts unzoomed.
    plotted: Option<std::sync::Weak<md_core::Chromatogram>>,
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
            tof: tof_calibration::TofWindow::default(),
            library: None,
            library_job: None,
            library_failed: None,
            nist: identify::Nist::default(),
            editor: library_editor::LibraryEditor::default(),
            batch: batch::Batch::default(),
            statistics: statistics::Statistics::default(),
            mid: mid::MidWindow::default(),
            results: results::Results::default(),
            export: export::Export::default(),
            export_all: export_all::ExportAll::default(),
            view: view::View::default(),
            spectrum_mode: detail::Mode::default(),
            updater: update::Updater::default(),
            panels: layout::Panels::default(),
            saved: persist::Saved::default(),
            model_profile: false,
            about_open: false,
            plotted: None,
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
    /// The peaks themselves (S/N, quality, sharpness), indexed like `peak_masses`.
    peaks: Vec<Peak>,
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
        self.shortcuts(ui.ctx());
        self.windows(ui.ctx());
        self.status_bar(ui);
        self.left_panel(ui);
        if self.panels.info
            && let Some(run) = self.runs.get(self.active)
            && detail::panel(
                ui,
                run,
                &self.nist,
                (&mut self.mz, &mut self.trace),
                self.ri.for_run(run),
            )
        {
            self.nist_search();
        }
        self.central_panel(ui);
        // Last, so the frame holds every plot that asked for an export.
        self.poll_export(ui.ctx());
    }

    fn save(&mut self, storage: &mut dyn eframe::Storage) {
        storage.set_string(persist::KEY, persist::text(self));
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
        self.update_nist(ctx);
    }

    /// The floating windows opened from the Analysis and Help menus.
    fn windows(&mut self, ctx: &egui::Context) {
        self.ri
            .window(ctx, &self.runs, self.active, self.library.as_deref());
        let libraries = [
            (self.library.as_ref()).map(|l| l.path.as_path()),
            (self.nist.lib.as_ref()).map(|l| l.path.as_path()),
        ];
        if let Some(which) = self.settings.window(ctx, libraries) {
            self.pick_library(ctx, which);
        }
        self.batch
            .use_settings(&self.settings.applied, self.ri.rt_only);
        // Not while the library is rebuilt after a switch of RT-only: its anchors would not match.
        let cal = (self.ri.cal.as_deref()).filter(|_| {
            (self.library.as_ref()).is_some_and(|l| l.prepared.is_rt_only() == self.ri.rt_only)
        });
        let reference = (self.ri.reference.as_deref()).map(|n| (n, self.ri.ref_cutoff));
        self.batch
            .window(ctx, self.library.as_ref(), cal, reference);
        self.mid.applied.clone_from(&self.settings.applied);
        self.mid.window(ctx, self.library.as_ref(), cal, reference);
        self.statistics
            .window(ctx, self.batch.shown_table(), self.batch.z_scored());
        self.isotopes.window(ctx);
        self.tof.window(ctx);
        self.info_window(ctx);
        self.library_editor(ctx);
        self.updater.window(ctx);
        self.about_window(ctx);
        self.export_all_window(ctx);
    }
}
