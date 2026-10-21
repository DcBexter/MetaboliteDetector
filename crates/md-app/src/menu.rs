//! Top menu bar: File, Analysis, View and Help.

use eframe::egui;

use crate::App;
use crate::lang::{t, tf};
use crate::shortcuts::{Action, REPO};

impl App {
    pub fn menu_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::top("menu").show(ui, |ui| {
            egui::MenuBar::new().ui(ui, |ui| {
                ui.menu_button(t("File", "Datei"), |ui| self.file_menu(ui));
                ui.menu_button(t("Analysis", "Analyse"), |ui| self.analysis_menu(ui));
                ui.menu_button(t("View", "Ansicht"), |ui| {
                    for (on, name, hint) in self.panels.toggles() {
                        ui.checkbox(on, name).on_hover_text(hint);
                    }
                });
                ui.menu_button(t("Help", "Hilfe"), |ui| self.help_menu(ui));
            });
        });
    }

    pub fn open_chromatograms_button(&mut self, ui: &mut egui::Ui) {
        self.action_item(
            ui,
            t("Open chromatograms…", "Chromatogramme öffnen…"),
            t(
                "Load one or more GC-MS runs; each is analysed in the background",
                "Einen oder mehrere GC-MS-Läufe laden; jeder wird im Hintergrund analysiert",
            ),
            Action::Open,
            true,
        );
    }

    pub fn pick_chromatograms(&mut self, ctx: &egui::Context) {
        if let Some(paths) = rfd::FileDialog::new()
            .add_filter("NetCDF (ANDI-MS)", &["cdf", "CDF", "nc"])
            .add_filter(
                t(
                    "Old MetaboliteDetector (.bin + .idx)",
                    "Alter MetaboliteDetector (.bin + .idx)",
                ),
                &["bin"],
            )
            .add_filter(
                t(
                    "FastFlight raw data (MsData.FFT)",
                    "FastFlight-Rohdaten (MsData.FFT)",
                ),
                &["FFT", "fft"],
            )
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
                t("Open a run first", "Erst einen Lauf öffnen")
            } else {
                t("Select a component first", "Erst eine Komponente wählen")
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
            t("Manual", "Handbuch"),
            t(
                "The README with screenshots, in the browser",
                "Das README mit Bildschirmfotos, im Browser",
            ),
            Action::Manual,
            true,
        );
        if ui
            .button(t("About MetaboliteDetector…", "Über MetaboliteDetector…"))
            .on_hover_text(t(
                "Version, licence, authors and how to cite",
                "Version, Lizenz, Autoren und Zitierweise",
            ))
            .clicked()
        {
            self.about_open = true;
        }
        ui.separator();
        if ui
            .button(t("Open log folder", "Log-Ordner öffnen"))
            .on_hover_text(t(
                "The folder where the app keeps its files; it writes no log yet",
                "Der Ordner, in dem das Programm seine Dateien ablegt; es schreibt noch kein Log",
            ))
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
        crate::lang::switch(ui);
    }

    pub fn about_window(&mut self, ctx: &egui::Context) {
        egui::Window::new(t("About MetaboliteDetector", "Über MetaboliteDetector"))
            .id(egui::Id::unique("About MetaboliteDetector"))
            .open(&mut self.about_open)
            .collapsible(false)
            .resizable(false)
            .show(ctx, |ui| {
                ui.heading(t("MetaboliteDetector", "MetaboliteDetector"));
                ui.label(tf!("Version {}", "Version {}", env!("CARGO_PKG_VERSION")));
                ui.label(t("MIT licence", "MIT-Lizenz"));
                ui.separator();
                ui.label(t(
                    "Rewrite (2026): by Milan Burgmann. Original program, version 3.3 (2009): by Milan Burgmann, and the authors named in the citation below, Karsten Hiller et al., TU Braunschweig.",
                    "Neuentwicklung (2026): von Milan Burgmann. Originalprogramm, Version 3.3 (2009): von Milan Burgmann und den unten im Zitat genannten Autoren, Karsten Hiller et al., TU Braunschweig.",
                ));
                ui.label(t("Please cite: Hiller K, Hangebrauk J, Jäger C, Spura J, Schreiber K, Schomburg D. MetaboliteDetector: Comprehensive Analysis Tool for Targeted and Nontargeted GC/MS Based Metabolome Analysis. Anal. Chem. 2009, 81(9), 3429–3439.", "Bitte zitieren: Hiller K, Hangebrauk J, Jäger C, Spura J, Schreiber K, Schomburg D. MetaboliteDetector: Comprehensive Analysis Tool for Targeted and Nontargeted GC/MS Based Metabolome Analysis. Anal. Chem. 2009, 81(9), 3429–3439."));
                ui.hyperlink_to(t("Source code", "Quellcode"), REPO);
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
            t("Save components…", "Komponenten speichern…"),
            t(
                "Writes the components of the active run as a .cmp file (old format)",
                "Schreibt die Komponenten des aktiven Laufs als .cmp-Datei (altes Format)",
            ),
            Action::SaveComponents,
            has_run,
        );
        self.action_item(
            ui,
            t("Export results as CSV…", "Ergebnisse als CSV exportieren…"),
            t(
                "Save the component list of the active run as CSV",
                "Die Komponentenliste des aktiven Laufs als CSV speichern",
            ),
            Action::ExportCsv,
            has_run,
        );
        self.action_item(
            ui,
            t(
                "Add component to library…",
                "Komponente zur Bibliothek hinzufügen…",
            ),
            t(
                "The selected component of the active run becomes a library entry",
                "Die gewählte Komponente des aktiven Laufs wird ein Bibliothekseintrag",
            ),
            Action::AddToLibrary,
            has_component,
        );
        if ui
            .add_enabled(
                !self.runs.is_empty(),
                egui::Button::new(t("Close all", "Alle schließen")),
            )
            .on_hover_text(t("Close every open run", "Jeden offenen Lauf schließen"))
            .on_disabled_hover_text(t("No runs are open", "Keine Läufe geöffnet"))
            .clicked()
        {
            (self.runs, self.active) = (Vec::new(), 0);
        }
        // Export needs a run; without one the section is empty and gets no second separator.
        if !self.runs.is_empty() {
            ui.separator();
        }
        if !self.runs.is_empty() {
            ui.menu_button(t("Export as NetCDF", "Als NetCDF exportieren"), |ui| {
                self.netcdf_menu(ui)
            });
        }
        if ui
            .add_enabled(
                !self.runs.is_empty(),
                egui::Button::new(t(
                    "Export all runs as NetCDF…",
                    "Alle Läufe als NetCDF exportieren…",
                )),
            )
            .on_hover_text(t(
                "Write every open run into one folder",
                "Jeden offenen Lauf in einen Ordner schreiben",
            ))
            .clicked()
        {
            self.export_all.open = true;
        }
        self.export_menu(ui);
        ui.separator();
        if ui
            .button(t("Library editor…", "Bibliothekseditor…"))
            .on_hover_text(t(
                "Create or change a library: compounds, RI, quantification ions",
                "Eine Bibliothek anlegen oder ändern: Substanzen, RI, Quantifizierungsionen",
            ))
            .clicked()
        {
            self.editor.open = true;
        }
        ui.separator();
        self.action_item(
            ui,
            t("Quit", "Beenden"),
            t("Close the program", "Das Programm schließen"),
            Action::Quit,
            true,
        );
    }

    /// The three ways the old program wrote a run: all points, one point per peak, or per nominal mass.
    fn netcdf_menu(&mut self, ui: &mut egui::Ui) {
        use md_core::centroid::{CentroidSettings, ExportKind, reduce};
        for kind in ExportKind::ALL {
            let (label, hint) = match kind {
                ExportKind::Profile => (
                    t("Profile…", "Profil…"),
                    t(
                        "Save every point of the active run as an ANDI-MS NetCDF file",
                        "Jeden Punkt des aktiven Laufs als ANDI-MS-NetCDF-Datei speichern",
                    ),
                ),
                ExportKind::Centroid => (
                    t("Centroid…", "Zentroid…"),
                    t(
                        "Save one point per detected peak (centroided spectra)",
                        "Einen Punkt pro erkanntem Peak speichern (zentroidierte Spektren)",
                    ),
                ),
                ExportKind::Nominal => (
                    t("Nominal mass…", "Nominalmasse…"),
                    t(
                        "Save the highest intensity per nominal mass (bins of 1)",
                        "Die höchste Intensität pro Nominalmasse speichern (Klassen der Breite 1)",
                    ),
                ),
            };
            let Some(r) = self.runs.get(self.active) else {
                return;
            };
            if ui.button(label).on_hover_text(hint).clicked()
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
                    t("Settings…", "Einstellungen…"),
                    t(
                        "Peak detection, deconvolution and library search parameters",
                        "Parameter für Peakerkennung, Dekonvolution und Bibliothekssuche",
                    ),
                    |a| &mut a.settings.open,
                    false,
                ),
                (
                    t("RI calibration…", "RI-Kalibrierung…"),
                    t(
                        "Retention index calibration from alkanes, a calibration mix or a .calib file",
                        "Retentionsindex-Kalibrierung aus Alkanen, einem Kalibriermix oder einer .calib-Datei",
                    ),
                    |a| &mut a.ri.open,
                    false,
                ),
            ],
            &[(
                t("Run info…", "Laufinfo…"),
                t(
                    "File, instrument and scan details of the active run",
                    "Datei-, Geräte- und Scan-Details des aktiven Laufs",
                ),
                |a| &mut a.view.info_open,
                true,
            )],
            &[
                (
                    t("Batch quantification…", "Batch-Quantifizierung…"),
                    t(
                        "Quantify the library's substances across many runs",
                        "Die Substanzen der Bibliothek über viele Läufe quantifizieren",
                    ),
                    |a| &mut a.batch.open,
                    false,
                ),
                (
                    t("Mass isotopomer distributions…", "Isotopomere (MID)…"),
                    t(
                        "Labelling experiments: isotopomer distributions of the library's fragments across runs",
                        "Markierungsexperimente: Isotopomerverteilungen der Fragmente der Bibliothek über Läufe",
                    ),
                    |a| &mut a.mid.open,
                    false,
                ),
                (
                    t("Biomarker statistics…", "Biomarker-Statistik…"),
                    t(
                        "Group comparison, box plots and PCA of the batch table",
                        "Gruppenvergleich, Boxplots und PCA der Batch-Tabelle",
                    ),
                    |a| &mut a.statistics.open,
                    false,
                ),
            ],
            &[
                (
                    t("Natural isotopes…", "Natürliche Isotope…"),
                    t(
                        "Isotope distribution (M+0, M+1, …) of a molecular formula",
                        "Isotopenverteilung (M+0, M+1, …) einer Summenformel",
                    ),
                    |a| &mut a.isotopes.open,
                    false,
                ),
                (
                    t("TOF calibration…", "TOF-Kalibrierung…"),
                    t(
                        "m/z calibration of FastFlight raw data from reference ions; saved as .mzcalib",
                        "m/z-Kalibrierung von FastFlight-Rohdaten aus Referenzionen; gespeichert als .mzcalib",
                    ),
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
                    .on_disabled_hover_text(t("Open a run first", "Erst einen Lauf öffnen"))
                    .clicked()
                {
                    *open(self) = true;
                }
            }
        }
    }
}
