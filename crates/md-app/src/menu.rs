//! Top menu bar: File, Analysis, View and Help.

use eframe::egui;

use crate::App;
use crate::shortcuts::{Action, REPO};

impl App {
    pub fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| self.file_menu(ui));
                ui.menu_button("Analysis", |ui| self.analysis_menu(ui));
                ui.menu_button("View", |ui| {
                    for (on, name, hint) in self.panels.toggles() {
                        ui.checkbox(on, name).on_hover_text(hint);
                    }
                });
                ui.menu_button("Help", |ui| self.help_menu(ui));
            });
        });
    }

    pub fn open_chromatograms_button(&mut self, ui: &mut egui::Ui) {
        self.action_item(
            ui,
            "Open chromatograms…",
            "Load one or more GC-MS runs; each is analysed in the background",
            Action::Open,
            true,
        );
    }

    pub fn pick_chromatograms(&mut self, ctx: &egui::Context) {
        if let Some(paths) = rfd::FileDialog::new()
            .add_filter("NetCDF (ANDI-MS)", &["cdf", "CDF", "nc"])
            .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
            .add_filter("FastFlight raw data (MsData.FFT)", &["FFT", "fft"])
            .pick_files()
        {
            for path in paths {
                self.open(ctx, path);
            }
        }
    }

    /// A menu entry that shows the shortcut of `action` and runs it on a click.
    fn action_item(
        &mut self,
        ui: &mut egui::Ui,
        label: &str,
        hint: &str,
        action: Action,
        enabled: bool,
    ) {
        let shortcut = ui.ctx().format_shortcut(&action.shortcut());
        let button = egui::Button::new(label).shortcut_text(shortcut);
        if ui
            .add_enabled(enabled, button)
            .on_hover_text(hint)
            .on_disabled_hover_text(if self.runs.is_empty() {
                "Open a run first"
            } else {
                "Select a component first"
            })
            .clicked()
        {
            // A click closes the menu by itself; outside a menu there is nothing to close.
            self.run_action(&ui.ctx().clone(), action);
        }
    }

    fn help_menu(&mut self, ui: &mut egui::Ui) {
        self.action_item(
            ui,
            "Manual",
            "The README with screenshots, in the browser",
            Action::Manual,
            true,
        );
        if ui
            .button("About MetaboliteDetector…")
            .on_hover_text("Version, licence, authors and how to cite")
            .clicked()
        {
            self.about_open = true;
        }
        ui.separator();
        if ui
            .button("Open log folder")
            .on_hover_text("The folder where the app keeps its files; it writes no log yet")
            .clicked()
            && let Some(dir) = crate::update::config_path()
                .and_then(|p| p.parent().map(std::path::Path::to_path_buf))
        {
            // Opening fails quietly if the folder does not exist yet; make it first.
            let _ = std::fs::create_dir_all(&dir);
            // open_url would start the web browser on a file: address on Windows.
            if cfg!(windows) {
                let _ = std::process::Command::new("explorer").arg(&dir).spawn();
            } else {
                ui.ctx()
                    .open_url(egui::OpenUrl::new_tab(dir.to_string_lossy()));
            }
        }
        self.updater.menu(ui);
        ui.separator();
        crate::guided::switch(ui, &mut self.guided.mode);
    }

    pub fn about_window(&mut self, ctx: &egui::Context) {
        egui::Window::new("About MetaboliteDetector")
            .open(&mut self.about_open)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.heading("MetaboliteDetector");
                ui.label(format!("Version {}", env!("CARGO_PKG_VERSION")));
                ui.label("MIT licence");
                ui.separator();
                ui.label("Rewrite (2026): by Milan Burgmann.Original program, version 3.3 (2009): by Milan Burgmann, and the authors named in the citation below, Karsten Hiller et al., TU Braunschweig.");
                ui.label("Please cite: Hiller K, Hangebrauk J, Jäger C, Spura J, Schreiber K, Schomburg D. MetaboliteDetector: Comprehensive Analysis Tool for Targeted and Nontargeted GC/MS Based Metabolome Analysis. Anal. Chem. 2009, 81(9), 3429–3439.");
                ui.hyperlink_to("Source code", REPO);
            });
    }

    fn file_menu(&mut self, ui: &mut egui::Ui) {
        self.open_chromatograms_button(ui);
        self.recent_menu(ui);
        self.open_library_button(ui);
        self.open_nist_button(ui);
        let has_run = !self.runs.is_empty();
        let has_component = (self.runs.get(self.active)).is_some_and(|r| r.component.is_some());
        self.action_item(
            ui,
            "Save components…",
            "Writes the components of the active run as a .cmp file (old format)",
            Action::SaveComponents,
            has_run,
        );
        self.action_item(
            ui,
            "Export results as CSV…",
            "Save the component list of the active run as CSV",
            Action::ExportCsv,
            has_run,
        );
        self.action_item(
            ui,
            "Add component to library…",
            "The selected component of the active run becomes a library entry",
            Action::AddToLibrary,
            has_component,
        );
        if ui
            .add_enabled(!self.runs.is_empty(), egui::Button::new("Close all"))
            .on_hover_text("Close every open run")
            .on_disabled_hover_text("No runs are open")
            .clicked()
        {
            (self.runs, self.active) = (Vec::new(), 0);
        }
        // Export needs a run; without one the section is empty and gets no second separator.
        if !self.runs.is_empty() {
            ui.separator();
        }
        if !self.runs.is_empty() {
            ui.menu_button("Export as NetCDF", |ui| self.netcdf_menu(ui));
        }
        if ui
            .add_enabled(
                !self.runs.is_empty(),
                egui::Button::new("Export all runs as NetCDF…"),
            )
            .on_hover_text("Write every open run into one folder")
            .clicked()
        {
            self.export_all.open = true;
        }
        self.export_menu(ui);
        ui.separator();
        if ui
            .button("Library editor…")
            .on_hover_text("Create or change a library: compounds, RI, quantification ions")
            .clicked()
        {
            self.editor.open = true;
        }
        ui.separator();
        self.action_item(ui, "Quit", "Close the program", Action::Quit, true);
    }

    /// The three ways the old program wrote a run: all points, one point per peak, or per nominal mass.
    fn netcdf_menu(&mut self, ui: &mut egui::Ui) {
        use md_core::centroid::{CentroidSettings, ExportKind, reduce};
        for kind in ExportKind::ALL {
            let hint = match kind {
                ExportKind::Profile => {
                    "Save every point of the active run as an ANDI-MS NetCDF file"
                }
                ExportKind::Centroid => "Save one point per detected peak (centroided spectra)",
                ExportKind::Nominal => "Save the highest intensity per nominal mass (bins of 1)",
            };
            let Some(r) = self.runs.get(self.active) else {
                return;
            };
            if ui
                .button(format!("{}…", kind.label()))
                .on_hover_text(hint)
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("NetCDF (ANDI-MS)", &["cdf"])
                    .set_file_name(format!("{}{}.cdf", r.name, kind.suffix()))
                    .save_file()
            {
                // shortcut: writes on the UI thread (~0.2 s per 10M points); spawn if runs grow.
                let result = if kind == ExportKind::Profile {
                    md_core::andi::write(&path, &r.data)
                } else {
                    let reduced = reduce(
                        &r.data,
                        kind,
                        &CentroidSettings::from(&self.settings.applied),
                    );
                    md_core::andi::write_with_tic(&path, &reduced, &r.data.tic())
                };
                self.error = result.err().map(|e| format!("{}: {e}", path.display()));
                ui.close();
            }
        }
    }

    /// Setup first, then the active run, then all runs, then the stand-alone calculator.
    fn analysis_menu(&mut self, ui: &mut egui::Ui) {
        type Toggle = fn(&mut App) -> &mut bool;
        // The flag marks entries that need an open run.
        let groups: [&[(&str, &str, Toggle, bool)]; 4] = [
            &[
                (
                    "Settings…",
                    "Peak detection, deconvolution and library search parameters",
                    |a| &mut a.settings.open,
                    false,
                ),
                (
                    "RI calibration…",
                    "Retention index calibration from alkanes, a calibration mix or a .calib file",
                    |a| &mut a.ri.open,
                    false,
                ),
            ],
            &[(
                "Run info…",
                "File, instrument and scan details of the active run",
                |a| &mut a.view.info_open,
                true,
            )],
            &[
                (
                    "Batch quantification…",
                    "Quantify the library's substances across many runs",
                    |a| &mut a.batch.open,
                    false,
                ),
                (
                    "Mass isotopomer distributions…",
                    "Labelling experiments: isotopomer distributions of the library's fragments across runs",
                    |a| &mut a.mid.open,
                    false,
                ),
                (
                    "Biomarker statistics…",
                    "Group comparison, box plots and PCA of the batch table",
                    |a| &mut a.statistics.open,
                    false,
                ),
            ],
            &[
                (
                    "Natural isotopes…",
                    "Isotope distribution (M+0, M+1, …) of a molecular formula",
                    |a| &mut a.isotopes.open,
                    false,
                ),
                (
                    "TOF calibration…",
                    "m/z calibration of FastFlight raw data from reference ions; saved as .mzcalib",
                    |a| &mut a.tof.open,
                    false,
                ),
            ],
        ];
        self.flows_menu(ui);
        for group in groups {
            ui.separator();
            for (name, hint, open, needs_run) in group {
                // Run info shows nothing without a run; the others explain their empty state.
                if ui
                    .add_enabled(
                        !needs_run || !self.runs.is_empty(),
                        egui::Button::new(*name),
                    )
                    .on_hover_text(*hint)
                    .on_disabled_hover_text("Open a run first")
                    .clicked()
                {
                    *open(self) = true;
                }
            }
        }
    }
}
