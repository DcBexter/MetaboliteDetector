//! Top menu bar: File and Analysis.

use eframe::egui;

use crate::App;

impl App {
    pub fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button("File", |ui| self.file_menu(ui));
                ui.menu_button("Analysis", |ui| self.analysis_menu(ui));
            });
        });
    }

    pub fn open_chromatograms_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .button("Open chromatograms…")
            .on_hover_text("Load one or more GC-MS runs; each is analysed in the background")
            .clicked()
            && let Some(paths) = rfd::FileDialog::new()
                .add_filter("NetCDF (ANDI-MS)", &["cdf", "CDF", "nc"])
                .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
                .add_filter("FastFlight raw data (MsData.FFT)", &["FFT", "fft"])
                .pick_files()
        {
            for path in paths {
                self.open(ui.ctx(), path);
            }
        }
    }

    fn file_menu(&mut self, ui: &mut egui::Ui) {
        self.open_chromatograms_button(ui);
        self.open_library_button(ui);
        self.open_nist_button(ui);
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
        if let Some(r) = self.runs.get(self.active)
            && ui
                .button("Export as NetCDF…")
                .on_hover_text("Save the active run as an ANDI-MS NetCDF file")
                .clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("NetCDF (ANDI-MS)", &["cdf"])
                .set_file_name(format!("{}.cdf", r.name))
                .save_file()
        {
            // shortcut: writes on the UI thread (~50 ms per example run); spawn if runs grow.
            self.error = md_core::andi::write(&path, &r.data)
                .err()
                .map(|e| format!("{}: {e}", path.display()));
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
            &[
                (
                    "Run info…",
                    "File, instrument and scan details of the active run",
                    |a| &mut a.view.info_open,
                    true,
                ),
                (
                    "Results of the active run…",
                    "Table of every component with its library hit; sortable, exportable as CSV",
                    |a| &mut a.results.open,
                    false,
                ),
            ],
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
            &[(
                "Natural isotopes…",
                "Isotope distribution (M+0, M+1, …) of a molecular formula",
                |a| &mut a.isotopes.open,
                false,
            )],
        ];
        for (i, group) in groups.into_iter().enumerate() {
            if i > 0 {
                ui.separator();
            }
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
