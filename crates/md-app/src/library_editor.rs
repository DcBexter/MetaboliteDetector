//! Library editor (F18): edits a copy of the open library, or a new one, and saves it as `.mdlib`
//! (docs/architecture.md, Libraries). It never writes `.lbr` or `.msl` files.

use std::path::{Path, PathBuf};

use crate::layout::ToolWindow;
use eframe::egui;
use egui_plot::{Bar, BarChart, Plot};
use md_core::library::{LibraryCompound, MslType, save_mdlib, write_msl_as};
use md_core::library_add::check;

mod add;
mod extras;

use crate::App;
use crate::errors;
use crate::identify::Library;
use crate::lang::{t, tf};
use add::{AddDialog, RtDialog, Source};
use extras::NewFragment;

#[derive(Default)]
pub struct LibraryEditor {
    pub open: bool,
    pub(crate) compounds: Vec<LibraryCompound>,
    /// `.mdlib` file the compounds came from or were last saved to; Save writes only there.
    path: Option<PathBuf>,
    /// Library name in the title; empty for a new library.
    name: String,
    selected: Option<usize>,
    filter: String,
    /// Quant ions as typed; they replace the compound's ions whenever they parse.
    ions: String,
    /// Key of the feature to add.
    new_key: String,
    /// Inputs of the "Add fragment" row.
    frag: NewFragment,
    dirty: bool,
    /// Replacement of the edited library that waits for "Discard" of unsaved changes.
    pending: Option<Replace>,
    error: Option<String>,
    /// Entries named by the last save or export because they are incomplete.
    warning: Option<String>,
    /// "Add component" and "Copy retention times" while they are open.
    add: Option<AddDialog>,
    rt: Option<RtDialog>,
}

#[derive(Clone, Copy)]
enum Replace {
    New,
    Open,
}

fn is_mdlib(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("mdlib"))
}

fn ions_text(ions: &[f32]) -> String {
    let ions: Vec<_> = ions.iter().map(f32::to_string).collect();
    ions.join(" ")
}

/// Quant ions are a sorted set like the old `std::set<float>`.
/// A "known" checkbox and, when ticked, a number field. Unknown is stored as -1 (the old program's
/// value); `unit` converts the stored value to what the field shows. True when the value changed.
fn optional_value(ui: &mut egui::Ui, value: &mut f64, unit: f64, speed: f64, hover: &str) -> bool {
    ui.horizontal(|ui| {
        let mut changed = false;
        let mut known = *value >= 0.0;
        if ui
            .checkbox(&mut known, t("known", "bekannt"))
            .on_hover_text(hover)
            .changed()
        {
            *value = if known { 0.0 } else { -1.0 };
            changed = true;
        }
        if known {
            let mut shown = *value / unit;
            let drag = egui::DragValue::new(&mut shown)
                .speed(speed)
                .range(0.0..=f64::MAX)
                .clamp_existing_to_range(false);
            if ui.add(drag).changed() {
                *value = shown * unit;
                changed = true;
            }
        }
        changed
    })
    .inner
}

fn parse_ions(text: &str) -> Option<Vec<f32>> {
    let mut ions = (text.split([' ', ',', ';']).filter(|s| !s.is_empty()))
        .map(|s| s.parse::<f32>().ok().filter(|m| m.is_finite() && *m > 0.0))
        .collect::<Option<Vec<_>>>()?;
    ions.sort_by(f32::total_cmp);
    ions.dedup();
    Some(ions)
}

impl LibraryEditor {
    pub(crate) fn select(&mut self, i: Option<usize>) {
        self.selected = i;
        self.ions = i.map_or(String::new(), |i| ions_text(&self.compounds[i].quant_ions));
    }

    /// Replaces the edited library, asking first when that would drop unsaved changes.
    fn ask(&mut self, what: Replace, library: Option<&Library>) {
        if self.dirty {
            self.pending = Some(what);
        } else {
            self.replace(what, library);
        }
    }

    fn replace(&mut self, what: Replace, library: Option<&Library>) {
        (self.compounds, self.path, self.name) = match (what, library) {
            (Replace::Open, Some(lib)) => (
                lib.compounds.clone(),
                Some(lib.path.clone()).filter(|p| is_mdlib(p)),
                lib.name.clone(),
            ),
            _ => (Vec::new(), None, String::new()),
        };
        (self.dirty, self.pending, self.error, self.warning) = (false, None, None, None);
        self.select(None);
    }

    /// Names the entries with problems (empty when all are fine); saving and exporting still write them.
    fn incomplete(&self) -> Option<String> {
        let names: Vec<_> = (self.compounds.iter())
            .filter(|c| !check(c).is_empty())
            .map(|c| c.name.as_str())
            .collect();
        const SHOWN: usize = 3;
        let more = names.len().saturating_sub(SHOWN);
        (!names.is_empty()).then(|| {
            let (n, shown) = (names.len(), names[..names.len().min(SHOWN)].join(", "));
            if more > 0 {
                tf!(
                    "{n} incomplete entries were written: {shown} and {more} more (filter the list to find them).",
                    "{n} unvollständige Einträge wurden geschrieben: {shown} und {more} weitere (Liste filtern, um sie zu finden).",
                )
            } else {
                tf!(
                    "{n} incomplete entries were written: {shown} (filter the list to find them).",
                    "{n} unvollständige Einträge wurden geschrieben: {shown} (Liste filtern, um sie zu finden).",
                )
            }
        })
    }

    fn save(&mut self, path: PathBuf) -> Option<PathBuf> {
        if let Err(e) = save_mdlib(&path, &self.compounds) {
            self.error = Some(format!("{}: {e}", path.display()));
            return None;
        }
        self.name = path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into();
        (self.dirty, self.error, self.path) = (false, None, Some(path.clone()));
        self.warning = self.incomplete();
        Some(path)
    }

    /// "Export MSL" with the three types of the old export dialog; only `.msl` files are written.
    /// `fwhm` is `SF_FWHM` of the centroiding.
    fn export_menu(&mut self, ui: &mut egui::Ui, stem: &str, fwhm: f64) {
        ui.menu_button(t("Export MSL", "MSL exportieren"), |ui| {
            for (kind, label) in [
                (MslType::Nominal, "MSL Nominal"),
                (MslType::Profile, "MSL Profile"),
                (MslType::Fwhm, "MSL FWHM"),
            ] {
                if ui.button(label).clicked() {
                    ui.close();
                    let Some(path) = rfd::FileDialog::new()
                        .add_filter(t("MSL library files", "MSL-Bibliotheksdateien"), &["msl"])
                        .set_file_name(format!("{stem}.msl"))
                        .save_file()
                    else {
                        return;
                    };
                    let text = write_msl_as(&self.compounds, kind, fwhm);
                    self.error = if path
                        .extension()
                        .is_some_and(|e| e.eq_ignore_ascii_case("msl"))
                    {
                        std::fs::write(&path, text)
                            .err()
                            .map(|e| format!("{}: {e}", path.display()))
                    } else {
                        let p = path.display();
                        Some(tf!("{p}: export as .msl.", "{p}: als .msl exportieren."))
                    };
                    self.warning = if self.error.is_none() {
                        self.incomplete()
                    } else {
                        None
                    };
                }
            }
        });
    }

    /// Draws the editor; returns the file it saved, if any.
    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        library: Option<&Library>,
        source: Option<Source>,
        fwhm: f64,
    ) -> Option<PathBuf> {
        let mut saved = None;
        // Wrapped: the German labels do not fit one line in a narrow window.
        ui.horizontal_wrapped(|ui| {
            if ui.button(t("New library", "Neue Bibliothek")).clicked() {
                self.ask(Replace::New, library);
            }
            if ui
                .add_enabled(
                    library.is_some(),
                    egui::Button::new(t("Edit open library", "Offene Bibliothek bearbeiten")),
                )
                .on_disabled_hover_text(t(
                    "Open a library first (File menu).",
                    "Erst eine Bibliothek öffnen (Menü Datei).",
                ))
                .clicked()
            {
                self.ask(Replace::Open, library);
            }
            if let Some(path) = self.path.clone()
                && ui
                    .button(t("Save", "Speichern"))
                    .on_hover_text(path.display().to_string())
                    .clicked()
            {
                saved = self.save(path);
            }
            let stem = Path::new(&self.name)
                .file_stem()
                .unwrap_or_default()
                .to_string_lossy()
                .into_owned();
            self.export_menu(ui, &stem, fwhm);
            if ui
                .button(t("Save as…", "Speichern unter…"))
                .on_hover_text(t(
                    "Write a .mdlib file. For an .msl file use Export MSL.",
                    "Schreibt eine .mdlib-Datei. Für eine .msl-Datei MSL exportieren verwenden.",
                ))
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .set_title(t(
                        "Save as .mdlib (for .msl use Export MSL in the editor)",
                        "Als .mdlib speichern (für .msl im Editor MSL exportieren verwenden)",
                    ))
                    .add_filter(t("MetaboliteDetector library", "MetaboliteDetector-Bibliothek"), &["mdlib"])
                    .set_file_name(format!("{stem}.mdlib"))
                    .save_file()
            {
                // The native dialog asked before replacing a file; only a .mdlib target is written.
                if is_mdlib(&path) {
                    saved = self.save(path);
                } else {
                    let p = path.display();
                    self.error = Some(tf!(
                        "{p}: save as .mdlib (use Export MSL for .msl); .lbr and .msl files are never overwritten.",
                        "{p}: als .mdlib speichern (für .msl MSL exportieren verwenden); .lbr- und .msl-Dateien werden nie überschrieben.",
                    ));
                }
            }
        });
        if let Some(what) = self.pending {
            ui.horizontal(|ui| {
                ui.colored_label(
                    ui.visuals().warn_fg_color,
                    t(
                        "Unsaved changes will be lost.",
                        "Ungespeicherte Änderungen gehen verloren.",
                    ),
                );
                if ui.button(t("Discard", "Verwerfen")).clicked() {
                    self.replace(what, library);
                }
                if ui.button(t("Cancel", "Abbrechen")).clicked() {
                    self.pending = None;
                }
            });
        }
        if let Some(e) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, errors::text(e));
        }
        if let Some(w) = &self.warning {
            ui.colored_label(ui.visuals().warn_fg_color, w);
        }
        ui.separator();
        self.add_windows(ui.ctx());

        crate::layout::list_and_details(
            ui,
            "library compounds",
            self,
            |ed, ui| ed.list(ui, source.as_ref()),
            |ed, ui| match ed.selected {
                Some(i) => {
                    egui::ScrollArea::vertical().show(ui, |ui| ed.details(ui, i));
                }
                None => {
                    ui.label(t("Select a compound.", "Eine Substanz auswählen."));
                }
            },
        );
        saved
    }

    fn list(&mut self, ui: &mut egui::Ui, source: Option<&Source>) {
        ui.horizontal_wrapped(|ui| {
            if ui.button(t("New compound", "Neue Substanz")).clicked() {
                self.compounds
                    // Data in the library file, so English like the file formats.
                    .push(LibraryCompound::new("New compound"));
                self.dirty = true;
                self.select(Some(self.compounds.len() - 1));
            }
            let src = source.filter(|s| s.component().is_some());
            if ui
                .add_enabled(
                    src.is_some(),
                    egui::Button::new(t("Add component…", "Komponente hinzufügen…")),
                )
                .on_hover_text(t(
                    "The selected component of the active run: spectrum, RT, RI and quant ions.",
                    "Die ausgewählte Komponente des aktiven Laufs: Spektrum, RT, RI und Quant-Ionen.",
                ))
                .clicked()
                && let Some(src) = src
            {
                self.start_add(src);
            }
            if ui
                .add_enabled(source.is_some(), egui::Button::new(t("Copy RTs…", "RTs kopieren…")))
                .on_hover_text(t(
                    "Copy RT and RI of the assigned components of the active run into this library.",
                    "Kopiert RT und RI der zugeordneten Komponenten des aktiven Laufs in diese Bibliothek.",
                ))
                .clicked()
                && let Some(src) = source
            {
                self.start_rt_copy(ui.ctx(), src);
            }
            if let Some(i) = self.selected
                && ui.button(t("Delete", "Löschen")).clicked()
            {
                self.compounds.remove(i);
                self.dirty = true;
                self.select(None);
            }
        });
        ui.add(
            egui::TextEdit::singleline(&mut self.filter)
                .hint_text(t("Filter by name", "Nach Name filtern"))
                .desired_width(f32::INFINITY),
        );
        let filter = self.filter.to_lowercase();
        let shown: Vec<usize> = (0..self.compounds.len())
            .filter(|&i| self.compounds[i].name.to_lowercase().contains(&filter))
            .collect();
        let (s, n) = (shown.len(), self.compounds.len());
        ui.label(tf!("{s} of {n} compounds", "{s} von {n} Substanzen"));
        let row = ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y;
        egui::ScrollArea::vertical().show_rows(ui, row, shown.len(), |ui, rows| {
            for &i in &shown[rows] {
                let c = &self.compounds[i];
                let mut name = egui::RichText::new(&c.name);
                let issues = check(c);
                if !issues.is_empty() {
                    name = name.color(ui.visuals().warn_fg_color);
                }
                let r = ui.selectable_label(self.selected == Some(i), name);
                let r = if issues.is_empty() {
                    r
                } else {
                    r.on_hover_text(
                        issues
                            .iter()
                            .map(ToString::to_string)
                            .collect::<Vec<_>>()
                            .join(
                                "
",
                            ),
                    )
                };
                if r.clicked() {
                    self.select(Some(i));
                }
            }
        });
    }

    fn details(&mut self, ui: &mut egui::Ui, i: usize) {
        let c = &mut self.compounds[i];
        let mut changed = false;
        egui::Grid::new("compound fields")
            .num_columns(2)
            .show(ui, |ui| {
                let name = ui.label(t("Name", "Name")).id;
                changed |= (ui.text_edit_singleline(&mut c.name))
                    .labelled_by(name)
                    .changed();
                ui.end_row();

                // The old program stores -1 for an unknown RT or RI; 0 is a real value.
                ui.label(t("RT (min)", "RT (min)"));
                changed |= optional_value(
                    ui,
                    &mut c.spectrum.rt_ms,
                    60_000.0,
                    0.01,
                    t(
                        "Untick when the retention time is unknown. 0 is a real value.",
                        "Abwählen, wenn die Retentionszeit unbekannt ist. 0 ist ein echter Wert.",
                    ),
                );
                ui.end_row();

                // The search treats a negative RI as unknown.
                ui.label("RI");
                changed |= optional_value(
                    ui,
                    &mut c.spectrum.ri,
                    1.0,
                    1.0,
                    t(
                        "Untick when the retention index is unknown. 0 is a real value.",
                        "Abwählen, wenn der Retentionsindex unbekannt ist. 0 ist ein echter Wert.",
                    ),
                );
                ui.end_row();

                ui.label(t("Resolution", "Auflösung"));
                ui.label(c.spectrum.bin_width.to_string()).on_hover_text(t(
                    "Bin width of the spectrum in m/z; the library file fixes it.",
                    "Bin-Breite des Spektrums in m/z; die Bibliotheksdatei legt sie fest.",
                ));
                ui.end_row();

                ui.label(t("Quant ions", "Quant-Ionen"));
                ui.horizontal(|ui| {
                    if ui
                        .add(
                            egui::TextEdit::singleline(&mut self.ions)
                                .hint_text(t("e.g. 174 248", "z. B. 174 248")),
                        )
                        .changed()
                        && let Some(ions) = parse_ions(&self.ions)
                    {
                        c.quant_ions = ions;
                        changed = true;
                    }
                    if parse_ions(&self.ions).is_none() {
                        ui.colored_label(
                            ui.visuals().error_fg_color,
                            t("m/z values > 0", "m/z-Werte > 0"),
                        );
                    }
                });
                ui.end_row();
                ui.label("");
                ui.vertical(|ui| {
                    ui.weak(t(
                        "m/z values separated by spaces, commas or semicolons.",
                        "m/z-Werte, getrennt durch Leerzeichen, Kommas oder Semikolons.",
                    ));
                    for issue in check(c) {
                        ui.colored_label(ui.visuals().warn_fg_color, issue.to_string());
                    }
                });
                ui.end_row();

                let fields = [
                    (
                        t("MID fragments", "MID-Fragmente"),
                        "FRAG",
                        false,
                        t(
                            "(m/z/carbons)…, e.g. (174/2)(248/2)",
                            "(m/z/Kohlenstoffe)…, z. B. (174/2)(248/2)",
                        ),
                    ),
                    (t("Comment", "Kommentar"), "COMMENT", true, ""),
                ];
                for (label, key, multiline, hint) in fields {
                    let name = ui.label(label).id;
                    if feature_field(ui, c, key, multiline, hint, name) {
                        changed = true;
                        // Clearing the field drops the feature instead of saving an empty one.
                        c.features.retain(|k, v| k != key || !v.is_empty());
                    }
                    ui.end_row();
                }
            });

        egui::CollapsingHeader::new(t("Other features", "Weitere Merkmale"))
            .id_salt("Other features")
            .show(ui, |ui| {
                changed |= other_features(ui, c, &mut self.new_key);
            });
        egui::CollapsingHeader::new(t("Automatic quant ions", "Automatische Quant-Ionen"))
            .id_salt("automatic quant ions")
            .show(ui, |ui| {
                changed |= extras::auto_ions(ui, c, &mut self.ions);
            })
            .header_response
            .on_hover_text(t(
                "Ions the program chose for this substance when it was added.",
                "Ionen, die das Programm beim Hinzufügen für diese Substanz gewählt hat.",
            ));
        egui::CollapsingHeader::new(t("Ion frequency", "Ionenhäufigkeit"))
            .id_salt("ion frequency")
            .show(ui, |ui| extras::ion_frequency(ui, c))
            .header_response
            .on_hover_text(t(
                "How often each ion occurred in the spectra this substance was built from.",
                "Wie oft jedes Ion in den Spektren vorkam, aus denen diese Substanz gebildet wurde.",
            ));
        egui::CollapsingHeader::new(t("MID fragments (sum formulas)", "MID-Fragmente (Summenformeln)"))
            .id_salt("mid fragments")
            .show(ui, |ui| {
                changed |= extras::mid_fragments(ui, c, &mut self.frag);
            })
            .header_response
            .on_hover_text(t(
                "Sum formula and MID size per fragment, read by the MID wizard. The FRAG field above is a separate list, used by the FRAG export.",
                "Summenformel und MID-Größe je Fragment, gelesen vom MID-Assistenten. Das Feld FRAG oben ist eine eigene Liste für den FRAG-Export.",
            ));

        let s = &c.spectrum;
        let n = s.masses.len();
        ui.label(tf!("Spectrum: {n} ions", "Spektrum: {n} Ionen"));
        let bars = (s.masses.iter().zip(&s.values))
            .map(|(m, v)| Bar::new(*m as f64, *v as f64).width(0.3))
            .collect();
        Plot::new("library spectrum")
            .height(220.0)
            .x_axis_label("m/z")
            .y_axis_label(t("Intensity", "Intensität"))
            .show(ui, |p| p.bar_chart(BarChart::new("library spectrum", bars)));
        self.dirty |= changed;
    }
}

/// Features other than FRAG and COMMENT: edit, remove, add by key.
fn other_features(ui: &mut egui::Ui, c: &mut LibraryCompound, new_key: &mut String) -> bool {
    let mut changed = false;
    let other: Vec<String> = (c.features.keys())
        .filter(|k| !matches!(k.as_str(), "FRAG" | "COMMENT"))
        .cloned()
        .collect();
    egui::Grid::new("features").striped(true).show(ui, |ui| {
        for k in other {
            let name = ui.label(&k).id;
            changed |= feature_field(ui, c, &k, false, "", name);
            if ui
                .small_button("×")
                .on_hover_text(t("Remove", "Entfernen"))
                .clicked()
            {
                c.features.remove(&k);
                changed = true;
            }
            ui.end_row();
        }
    });
    ui.horizontal(|ui| {
        ui.add(egui::TextEdit::singleline(new_key).hint_text(t("Key", "Schlüssel")));
        let key = new_key.trim().to_uppercase();
        if ui
            .add_enabled(
                // FRAG and COMMENT have their own fields above.
                !matches!(key.as_str(), "" | "FRAG" | "COMMENT") && !c.features.contains_key(&key),
                egui::Button::new(t("Add feature", "Merkmal hinzufügen")),
            )
            .clicked()
        {
            c.features.insert(key, String::new());
            new_key.clear();
            changed = true;
        }
    });
    changed
}

/// Edits one feature's value, adding the feature on the first edit; `name` is the row's label.
fn feature_field(
    ui: &mut egui::Ui,
    c: &mut LibraryCompound,
    key: &str,
    multiline: bool,
    hint: &str,
    name: egui::Id,
) -> bool {
    let mut value = c.features.get(key).cloned().unwrap_or_default();
    let edit = if multiline {
        egui::TextEdit::multiline(&mut value).desired_rows(3)
    } else {
        egui::TextEdit::singleline(&mut value)
    };
    let changed = ui.add(edit.hint_text(hint)).labelled_by(name).changed();
    if changed {
        c.features.insert(key.into(), value);
    }
    changed
}

impl App {
    /// Opens "Add component" and "Copy RTs" for the active run, as their buttons do; for tests.
    #[cfg(test)]
    pub fn open_editor_dialogs(&mut self, ctx: &egui::Context) {
        let run = &self.runs[self.active];
        let src = Source {
            run,
            ri: self.ri.for_run(run),
            quant: self.settings.applied.quant(),
            search: self.settings.applied.search(),
        };
        self.editor.start_add(&src);
        self.editor.start_rt_copy(ctx, &src);
    }

    /// Ctrl+L: opens the editor and "Add component" for the selected component of the active run.
    pub fn add_component_to_library(&mut self) {
        let Some(run) = self.runs.get(self.active) else {
            return;
        };
        let source = Source {
            run,
            ri: self.ri.for_run(run),
            quant: self.settings.applied.quant(),
            search: self.settings.applied.search(),
        };
        if source.component().is_some() {
            self.editor.open = true;
            self.editor.start_add(&source);
        }
    }

    pub fn library_editor(&mut self, ctx: &egui::Context) {
        if !self.editor.open {
            return;
        }
        let run = self.runs.get(self.active);
        let source = run.map(|run| Source {
            run,
            ri: self.ri.for_run(run),
            quant: self.settings.applied.quant(),
            search: self.settings.applied.search(),
        });
        let fwhm = self.settings.applied.sf_fwhm;
        let ed = &mut self.editor;
        let name = if ed.name.is_empty() {
            t("new library", "neue Bibliothek")
        } else {
            &ed.name
        };
        let title = if ed.dirty {
            tf!(
                "Library editor — {name} (unsaved)",
                "Bibliothekseditor — {name} (ungespeichert)"
            )
        } else {
            tf!("Library editor — {name}", "Bibliothekseditor — {name}")
        };
        let (mut open, mut saved) = (true, None);
        egui::Window::new(title)
            .tool_window(ctx, 8, [760.0, 520.0])
            .id(egui::Id::unique("library editor"))
            .open(&mut open)
            .show(ctx, |ui| {
                saved = ed.contents(ui, self.library.as_deref(), source, fwhm);
            });
        ed.open = open;
        // The app searches with the saved library from now on.
        if let Some(path) = saved {
            self.open_library(ctx, path);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_core::library::read_mdlib;
    use md_core::search::{PreparedLibrary, SearchSettings};

    fn library(path: &str) -> Library {
        let mut c = LibraryCompound::new("Alanine 2TMS");
        (c.spectrum.masses, c.spectrum.values) = (vec![73.0, 116.0], vec![1000, 800]);
        let s = SearchSettings::default();
        Library {
            path: path.into(),
            name: path.into(),
            prepared: PreparedLibrary::new(std::slice::from_ref(&c), &s),
            settings: s,
            compounds: vec![c],
        }
    }

    #[test]
    fn quant_ions_parse_as_a_sorted_set() {
        assert_eq!(
            parse_ions("248, 174 174;73.5"),
            Some(vec![73.5, 174.0, 248.0])
        );
        assert_eq!(parse_ions(""), Some(vec![]));
        assert_eq!(parse_ions("174 x"), None);
        assert_eq!(parse_ions("-1"), None);
        assert_eq!(parse_ions("NaN"), None);
    }

    #[test]
    fn save_targets_only_mdlib_files() {
        let mut ed = LibraryEditor::default();
        ed.replace(Replace::Open, Some(&library("old.lbr")));
        assert_eq!((ed.path.as_ref(), ed.compounds.len()), (None, 1));
        ed.replace(Replace::Open, Some(&library("mine.MDLIB")));
        assert_eq!(ed.path, Some(PathBuf::from("mine.MDLIB")));

        // Unsaved edits wait for Discard before another library replaces them.
        ed.dirty = true;
        ed.ask(Replace::New, None);
        assert_eq!(ed.compounds.len(), 1);
        ed.replace(ed.pending.unwrap(), None);
        assert!(ed.compounds.is_empty() && ed.path.is_none() && !ed.dirty);
    }

    #[test]
    fn edited_library_saves_and_reopens() {
        let mut ed = LibraryEditor::default();
        ed.replace(Replace::Open, Some(&library("old.lbr")));
        let c = &mut ed.compounds[0];
        c.spectrum.ri = 1105.3;
        c.quant_ions = parse_ions("116 73").unwrap();
        c.features.insert("FRAG".into(), "(116/2)".into());
        c.features.insert("COMMENT".into(), "checked".into());
        ed.dirty = true;
        let path = std::env::temp_dir().join(format!("md-app-{}.mdlib", std::process::id()));
        ed.compounds.push(LibraryCompound::new("Empty"));
        let saved = ed.save(path.clone());
        assert!(ed.warning.as_ref().is_some_and(|w| w.contains("Empty")));
        let again = read_mdlib(&path);
        std::fs::remove_file(&path).unwrap();
        assert_eq!(saved.as_ref(), Some(&path));
        assert_eq!((ed.path.as_ref(), ed.dirty), (Some(&path), false));
        assert_eq!(again.unwrap(), ed.compounds);
    }

    #[test]
    fn saving_names_incomplete_entries() {
        let mut ed = LibraryEditor::default();
        ed.compounds.push(LibraryCompound::new("Test"));
        assert!(ed.incomplete().is_some_and(|w| w.contains("Test")));
        for i in 0..5 {
            ed.compounds.push(LibraryCompound::new(&format!("N{i}")));
        }
        let w = ed.incomplete().unwrap();
        assert!(w.starts_with("6 incomplete") && w.contains("and 3 more") && !w.contains("N4"));
        ed.compounds.clear();
        assert_eq!(ed.incomplete(), None);
    }

    #[test]
    fn viewing_a_compound_changes_nothing() {
        let mut app = App::default();
        app.editor.replace(Replace::Open, Some(&library("old.lbr")));
        // Values the fields could clamp or drop just by being drawn.
        let c = &mut app.editor.compounds[0];
        c.spectrum.rt_ms = -1.0;
        c.features.insert("COMMENT".into(), String::new());
        let before = app.editor.compounds.clone();
        app.editor.open = true;
        app.editor.select(Some(0));
        let ctx = egui::Context::default();
        for _ in 0..2 {
            let mut out = ctx.run_ui(egui::RawInput::default(), |ui| app.library_editor(ui.ctx()));
            out.textures_delta.clear();
        }
        assert!(app.editor.open && !app.editor.dirty);
        assert_eq!(app.editor.compounds, before);
    }
}
