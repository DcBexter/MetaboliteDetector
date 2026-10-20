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

    fn file_menu(&mut self, ui: &mut egui::Ui) {
        if ui.button("Open chromatograms…").clicked()
            && let Some(paths) = rfd::FileDialog::new()
                .add_filter("NetCDF (ANDI-MS)", &["cdf", "CDF", "nc"])
                .add_filter("Old MetaboliteDetector (.bin + .idx)", &["bin"])
                .pick_files()
        {
            for path in paths {
                self.open(ui.ctx(), path);
            }
        }
        if let Some(r) = self.runs.get(self.active)
            && ui.button("Export as NetCDF…").clicked()
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
        self.library_menu(ui);
        if ui.button("Close all").clicked() {
            (self.runs, self.active) = (Vec::new(), 0);
        }
    }

    fn analysis_menu(&mut self, ui: &mut egui::Ui) {
        if ui.button("RI calibration…").clicked() {
            self.ri.open = true;
        }
        if ui.button("Settings…").clicked() {
            self.settings.open = true;
        }
        if ui.button("Batch quantification…").clicked() {
            self.batch.open = true;
        }
        if ui.button("Run info…").clicked() {
            self.view.info_open = true;
        }
        if ui.button("Natural isotopes…").clicked() {
            self.isotopes.open = true;
        }
        if ui.button("Results of the active run…").clicked() {
            self.results.open = true;
        }
    }
}
