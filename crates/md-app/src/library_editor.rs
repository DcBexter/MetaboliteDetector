//! Library editor (F18): edits a copy of the open library, or a new one, and saves it as `.mdlib`
//! (docs/architecture.md, Libraries). It never writes `.lbr` or `.msl` files.

use std::path::{Path, PathBuf};

use eframe::egui;
use egui_plot::{Bar, BarChart, Plot};
use md_core::library::{LibraryCompound, save_mdlib};

use crate::App;
use crate::identify::Library;

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
    dirty: bool,
    /// Replacement of the edited library that waits for "Discard" of unsaved changes.
    pending: Option<Replace>,
    error: Option<String>,
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
        (self.dirty, self.pending, self.error) = (false, None, None);
        self.select(None);
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
        Some(path)
    }

    /// Draws the editor; returns the file it saved, if any.
    fn contents(
        &mut self,
        ui: &mut egui::Ui,
        library: Option<&Library>,
        component: Option<LibraryCompound>,
    ) -> Option<PathBuf> {
        let mut saved = None;
        ui.horizontal(|ui| {
            if ui.button("New library").clicked() {
                self.ask(Replace::New, library);
            }
            if ui
                .add_enabled(library.is_some(), egui::Button::new("Edit open library"))
                .on_disabled_hover_text("Open a library first (File menu).")
                .clicked()
            {
                self.ask(Replace::Open, library);
            }
            if let Some(path) = self.path.clone()
                && ui
                    .button("Save")
                    .on_hover_text(path.display().to_string())
                    .clicked()
            {
                saved = self.save(path);
            }
            let stem = Path::new(&self.name).file_stem().unwrap_or_default();
            if ui.button("Save as…").clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("MetaboliteDetector library", &["mdlib"])
                    .set_file_name(format!("{}.mdlib", stem.to_string_lossy()))
                    .save_file()
            {
                // The native dialog asked before replacing a file; only a .mdlib target is written.
                if is_mdlib(&path) {
                    saved = self.save(path);
                } else {
                    self.error = Some(format!(
                        "{}: save as .mdlib; .lbr and .msl files are never overwritten.",
                        path.display()
                    ));
                }
            }
        });
        if let Some(what) = self.pending {
            ui.horizontal(|ui| {
                ui.colored_label(ui.visuals().warn_fg_color, "Unsaved changes will be lost.");
                if ui.button("Discard").clicked() {
                    self.replace(what, library);
                }
                if ui.button("Cancel").clicked() {
                    self.pending = None;
                }
            });
        }
        if let Some(e) = &self.error {
            ui.colored_label(ui.visuals().error_fg_color, e);
        }
        ui.separator();

        crate::layout::list_and_details(
            ui,
            "library compounds",
            self,
            |ed, ui| ed.list(ui, component),
            |ed, ui| match ed.selected {
                Some(i) => {
                    egui::ScrollArea::vertical().show(ui, |ui| ed.details(ui, i));
                }
                None => {
                    ui.label("Select a compound.");
                }
            },
        );
        saved
    }

    fn list(&mut self, ui: &mut egui::Ui, component: Option<LibraryCompound>) {
        ui.horizontal_wrapped(|ui| {
            if ui.button("New compound").clicked() {
                self.compounds.push(LibraryCompound::new("New compound"));
                self.dirty = true;
                self.select(Some(self.compounds.len() - 1));
            }
            if ui
                .add_enabled(component.is_some(), egui::Button::new("Add component"))
                .on_hover_text("The selected component of the active run: spectrum, RT and RI.")
                .clicked()
                && let Some(c) = component
            {
                self.compounds.push(c);
                self.dirty = true;
                self.select(Some(self.compounds.len() - 1));
            }
            if let Some(i) = self.selected
                && ui.button("Delete").clicked()
            {
                self.compounds.remove(i);
                self.dirty = true;
                self.select(None);
            }
        });
        ui.add(
            egui::TextEdit::singleline(&mut self.filter)
                .hint_text("Filter by name")
                .desired_width(f32::INFINITY),
        );
        let filter = self.filter.to_lowercase();
        let shown: Vec<usize> = (0..self.compounds.len())
            .filter(|&i| self.compounds[i].name.to_lowercase().contains(&filter))
            .collect();
        ui.label(format!(
            "{} of {} compounds",
            shown.len(),
            self.compounds.len()
        ));
        let row = ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y;
        egui::ScrollArea::vertical().show_rows(ui, row, shown.len(), |ui, rows| {
            for &i in &shown[rows] {
                let name = &self.compounds[i].name;
                if ui
                    .selectable_label(self.selected == Some(i), name)
                    .clicked()
                {
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
                ui.label("Name");
                changed |= ui.text_edit_singleline(&mut c.name).changed();
                ui.end_row();

                ui.label("RT (min)");
                let mut minutes = c.spectrum.rt_ms / 60_000.0;
                if ui
                    .add(
                        egui::DragValue::new(&mut minutes)
                            .speed(0.01)
                            .range(0.0..=f64::MAX)
                            .clamp_existing_to_range(false),
                    )
                    .changed()
                {
                    c.spectrum.rt_ms = minutes * 60_000.0;
                    changed = true;
                }
                ui.end_row();

                // The search treats a negative RI as unknown.
                ui.label("RI");
                ui.horizontal(|ui| {
                    let mut known = c.spectrum.ri >= 0.0;
                    if ui.checkbox(&mut known, "known").changed() {
                        c.spectrum.ri = if known { 0.0 } else { -1.0 };
                        changed = true;
                    }
                    if known {
                        changed |= ui
                            .add(
                                egui::DragValue::new(&mut c.spectrum.ri)
                                    .range(0.0..=f64::MAX)
                                    .clamp_existing_to_range(false),
                            )
                            .changed();
                    }
                });
                ui.end_row();

                ui.label("Quant ions");
                ui.horizontal(|ui| {
                    if ui
                        .add(egui::TextEdit::singleline(&mut self.ions).hint_text("e.g. 174 248"))
                        .changed()
                        && let Some(ions) = parse_ions(&self.ions)
                    {
                        c.quant_ions = ions;
                        changed = true;
                    }
                    if parse_ions(&self.ions).is_none() {
                        ui.colored_label(ui.visuals().error_fg_color, "m/z values > 0");
                    }
                });
                ui.end_row();

                let fields = [
                    (
                        "MID fragments",
                        "FRAG",
                        false,
                        "(m/z/carbons)…, e.g. (174/2)(248/2)",
                    ),
                    ("Comment", "COMMENT", true, ""),
                ];
                for (label, key, multiline, hint) in fields {
                    ui.label(label);
                    if feature_field(ui, c, key, multiline, hint) {
                        changed = true;
                        // Clearing the field drops the feature instead of saving an empty one.
                        c.features.retain(|k, v| k != key || !v.is_empty());
                    }
                    ui.end_row();
                }
            });

        ui.collapsing("Other features", |ui| {
            changed |= other_features(ui, c, &mut self.new_key);
        });

        let s = &c.spectrum;
        ui.label(format!("Spectrum: {} ions", s.masses.len()));
        let bars = (s.masses.iter().zip(&s.values))
            .map(|(m, v)| Bar::new(*m as f64, *v as f64).width(0.3))
            .collect();
        Plot::new("library spectrum")
            .height(220.0)
            .x_axis_label("m/z")
            .y_axis_label("Intensity")
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
            ui.label(&k);
            changed |= feature_field(ui, c, &k, false, "");
            if ui.small_button("✕").on_hover_text("Remove").clicked() {
                c.features.remove(&k);
                changed = true;
            }
            ui.end_row();
        }
    });
    ui.horizontal(|ui| {
        ui.add(egui::TextEdit::singleline(new_key).hint_text("Key"));
        let key = new_key.trim().to_uppercase();
        if ui
            .add_enabled(
                // FRAG and COMMENT have their own fields above.
                !matches!(key.as_str(), "" | "FRAG" | "COMMENT") && !c.features.contains_key(&key),
                egui::Button::new("Add feature"),
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

/// Edits one feature's value, adding the feature on the first edit.
fn feature_field(
    ui: &mut egui::Ui,
    c: &mut LibraryCompound,
    key: &str,
    multiline: bool,
    hint: &str,
) -> bool {
    let mut value = c.features.get(key).cloned().unwrap_or_default();
    let edit = if multiline {
        egui::TextEdit::multiline(&mut value).desired_rows(3)
    } else {
        egui::TextEdit::singleline(&mut value)
    };
    let changed = ui.add(edit.hint_text(hint)).changed();
    if changed {
        c.features.insert(key.into(), value);
    }
    changed
}

impl App {
    pub fn library_editor(&mut self, ctx: &egui::Context) {
        if !self.editor.open {
            return;
        }
        let component = self.selected_component();
        let ed = &mut self.editor;
        let name = if ed.name.is_empty() {
            "new library"
        } else {
            &ed.name
        };
        let title = format!(
            "Library editor — {name}{}",
            if ed.dirty { " (unsaved)" } else { "" }
        );
        let (mut open, mut saved) = (true, None);
        egui::Window::new(title)
            .id(egui::Id::unique("library editor"))
            .open(&mut open)
            .default_size([760.0, 520.0])
            .show(ctx, |ui| {
                saved = ed.contents(ui, self.library.as_deref(), component);
            });
        ed.open = open;
        // The app searches with the saved library from now on.
        if let Some(path) = saved {
            self.open_library(ctx, path);
        }
    }

    /// The selected component of the active run as a new library compound.
    fn selected_component(&self) -> Option<LibraryCompound> {
        let run = self.runs.get(self.active)?;
        let k = &run.components[run.component?];
        let rt_s = run.data.scan_times_s[k.apex_scan as usize];
        let mut c = LibraryCompound::new(&format!("{} {:.3} min", run.name, rt_s / 60.0));
        let s = &mut c.spectrum;
        (s.scan, s.rt_ms) = (k.apex_scan, rt_s * 1000.0);
        s.ri = self.ri.for_run(run).map_or(-1.0, |cal| cal.ri(rt_s));
        (s.masses, s.values) = (k.spectrum.iter())
            .map(|&(m, h)| (m as f32, i32::try_from(h).unwrap_or(i32::MAX)))
            .unzip();
        // Like an MSL import: one source spectrum, so it is also the combined spectrum.
        c.comb_spec = (s.masses.iter().copied())
            .zip(s.values.iter().copied())
            .collect();
        c.ion_freq = s.masses.iter().map(|&m| (m, 1)).collect();
        c.sources.push(format!("Added from {}", run.name));
        Some(c)
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
        let saved = ed.save(path.clone());
        let again = read_mdlib(&path);
        std::fs::remove_file(&path).unwrap();
        assert_eq!(saved.as_ref(), Some(&path));
        assert_eq!((ed.path.as_ref(), ed.dirty), (Some(&path), false));
        assert_eq!(again.unwrap(), ed.compounds);
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
