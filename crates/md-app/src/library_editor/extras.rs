//! Compound fields that came with the `.lbr` versions 3 to 5: automatic quant ions, ion frequency
//! and the MID fragments (sum formulas and MID sizes) the MID wizard reads.

use eframe::egui;
use md_core::isotope::parse_formula;
use md_core::library::LibraryCompound;

use super::ions_text;

/// Inputs of the row that adds a MID fragment.
#[derive(Default)]
pub struct NewFragment {
    begin: u32,
    end: u32,
    formula: String,
    size: u32,
}

const FORMULA_WIDTH: f32 = 120.0;

/// What is wrong with a sum formula, if anything.
fn formula_error(text: &str) -> Option<String> {
    if text.trim().is_empty() {
        return Some("Enter a sum formula, e.g. C7H20NSi2".into());
    }
    parse_formula(text.trim()).err()
}

/// Automatic quant ions: show them, take them over as fixed quant ions or clear them.
/// `ions` is the quant ion text of the editor, kept in step. True when the compound changed.
pub fn auto_ions(ui: &mut egui::Ui, c: &mut LibraryCompound, ions: &mut String) -> bool {
    if c.auto_quant_ions.is_empty() {
        ui.label("None");
        return false;
    }
    let mut changed = false;
    ui.label(ions_text(&c.auto_quant_ions));
    ui.horizontal(|ui| {
        if ui
            .button("Use as quant ions")
            .on_hover_text("Add them to the quant ions of this substance.")
            .clicked()
        {
            c.quant_ions.extend_from_slice(&c.auto_quant_ions);
            c.quant_ions.sort_by(f32::total_cmp);
            c.quant_ions.dedup();
            *ions = ions_text(&c.quant_ions);
            changed = true;
        }
        if ui
            .button("Clear automatic ions")
            .on_hover_text("Forget the ions the program chose when it added the substance.")
            .clicked()
        {
            c.auto_quant_ions.clear();
            changed = true;
        }
    });
    changed
}

/// How often each ion occurred in the spectra this substance was built from (read only).
pub fn ion_frequency(ui: &mut egui::Ui, c: &LibraryCompound) {
    ui.label(format!("Built from {} spectra.", c.spec_count));
    if c.ion_freq.is_empty() {
        return;
    }
    let text: Vec<_> = (c.ion_freq.iter())
        .map(|(m, n)| format!("{m} ×{n}"))
        .collect();
    ui.label(text.join("   "));
}

/// The MID size stored for a fragment, 0 (one share per ion) when none is.
fn mid_size(c: &LibraryCompound, begin: u32, end: u32) -> u32 {
    (c.mid_sizes.iter())
        .find(|m| (m.0, m.1) == (begin, end))
        .map_or(0, |m| m.2)
}

fn set_mid_size(c: &mut LibraryCompound, begin: u32, end: u32, size: u32) {
    match c.mid_sizes.iter_mut().find(|m| (m.0, m.1) == (begin, end)) {
        Some(m) => m.2 = size,
        None => c.mid_sizes.push((begin, end, size)),
    }
}

fn has_fragment(c: &LibraryCompound, begin: u32, end: u32) -> bool {
    c.formulas.iter().any(|f| (f.0, f.1) == (begin, end))
}

fn ion_field(ui: &mut egui::Ui, value: &mut u32, name: egui::Id) -> bool {
    let drag = egui::DragValue::new(value).range(0..=u32::MAX);
    ui.add(drag).labelled_by(name).changed()
}

/// One fragment: changing start or end moves its MID size with it, and a change that would repeat
/// another fragment's start and end is dropped. `names` are the column header labels.
fn fragment_row(
    ui: &mut egui::Ui,
    c: &mut LibraryCompound,
    i: usize,
    names: [egui::Id; 4],
) -> bool {
    let (b0, e0) = (c.formulas[i].0, c.formulas[i].1);
    let (mut b, mut e, mut size) = (b0, e0, mid_size(c, b0, e0));
    let mut changed = ion_field(ui, &mut b, names[0]) | ion_field(ui, &mut e, names[1]);
    if (b, e) != (b0, e0) {
        if has_fragment(c, b, e) {
            (b, e) = (b0, e0);
            changed = false;
        } else {
            let size = mid_size(c, b0, e0);
            c.formulas[i].0 = b;
            c.formulas[i].1 = e;
            c.mid_sizes.retain(|m| (m.0, m.1) != (b0, e0));
            set_mid_size(c, b, e, size);
        }
    }
    let error = formula_error(&c.formulas[i].2);
    let mut edit = egui::TextEdit::singleline(&mut c.formulas[i].2);
    if error.is_some() {
        edit = edit.text_color(ui.visuals().error_fg_color);
    }
    // A grid cell is only as wide as its header, so the width is set here.
    let field = [FORMULA_WIDTH, ui.spacing().interact_size.y];
    let r = ui.add_sized(field, edit).labelled_by(names[2]);
    changed |= r.changed();
    if let Some(error) = error {
        r.on_hover_text(error);
    }
    let drag = egui::DragValue::new(&mut size).range(0..=u32::MAX);
    if ui.add(drag).labelled_by(names[3]).changed() {
        set_mid_size(c, b, e, size);
        changed = true;
    }
    changed
}

/// The MID fragments of a substance: edit, remove or add `start / end / formula / MID size`.
/// True when the compound changed.
pub fn mid_fragments(ui: &mut egui::Ui, c: &mut LibraryCompound, new: &mut NewFragment) -> bool {
    let mut changed = false;
    let mut remove = None;
    if !c.formulas.is_empty() {
        egui::Grid::new("mid fragments")
            .striped(true)
            .show(ui, |ui| {
                let hints = [
                    "First ion of the fragment",
                    "Last ion of the fragment",
                    "Sum formula of the fragment, e.g. C7H20NSi2",
                    "Number of isotopomers in the MID; 0 means one share per ion",
                ];
                let heads = ["Start", "End", "Formula", "MID size"];
                let names = std::array::from_fn(|k| ui.strong(heads[k]).on_hover_text(hints[k]).id);
                ui.end_row();
                for i in 0..c.formulas.len() {
                    changed |= fragment_row(ui, c, i, names);
                    if ui.small_button("×").on_hover_text("Remove").clicked() {
                        remove = Some(i);
                    }
                    ui.end_row();
                }
            });
    }
    if let Some(i) = remove {
        let (b, e, _) = c.formulas.remove(i);
        c.mid_sizes.retain(|m| (m.0, m.1) != (b, e));
        changed = true;
    }
    ui.horizontal(|ui| {
        ui.add(egui::DragValue::new(&mut new.begin).prefix("start "));
        ui.add(egui::DragValue::new(&mut new.end).prefix("end "));
        ui.add(
            egui::TextEdit::singleline(&mut new.formula)
                .desired_width(FORMULA_WIDTH)
                .hint_text("C7H20NSi2"),
        );
        ui.add(egui::DragValue::new(&mut new.size).prefix("size "));
        let problem = if new.begin == 0 || new.end == 0 {
            Some("Start and end are ions above 0.".to_string())
        } else if has_fragment(c, new.begin, new.end) {
            Some("A fragment with this start and end exists.".to_string())
        } else {
            formula_error(&new.formula)
        };
        let add = ui.add_enabled(problem.is_none(), egui::Button::new("Add fragment"));
        if add.clicked() {
            c.formulas
                .push((new.begin, new.end, new.formula.trim().to_string()));
            c.mid_sizes.push((new.begin, new.end, new.size));
            c.formulas.sort_by_key(|f| (f.0, f.1));
            c.mid_sizes.sort_by_key(|m| (m.0, m.1));
            *new = NewFragment::default();
            changed = true;
        } else if let Some(problem) = problem {
            add.on_disabled_hover_text(problem);
        }
    });
    changed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formulas_are_checked_with_the_isotope_parser() {
        assert!(formula_error("C7H20NSi2").is_none());
        assert!(formula_error("").is_some());
        assert!(formula_error("c7").is_some());
    }

    /// Adds a fragment through the dialog, saves, reloads, then removes it again.
    #[test]
    fn fragments_are_added_saved_reloaded_and_removed() {
        use egui_kittest::Harness;
        use egui_kittest::kittest::Queryable;
        let mut c = LibraryCompound::new("x");
        c.formulas.push((200, 210, "C6H12O2".into()));
        c.mid_sizes.push((200, 210, 6));
        let new = NewFragment {
            begin: 73,
            end: 147,
            formula: " C7H20NSi2 ".into(),
            size: 4,
        };
        let mut h = Harness::new_ui_state(
            |ui, s: &mut (LibraryCompound, NewFragment)| {
                mid_fragments(ui, &mut s.0, &mut s.1);
            },
            (c, new),
        );
        h.run();
        h.get_by_label("Add fragment").click();
        h.run();
        let c = h.state().0.clone();
        // Sorted by start, formula trimmed, MID size kept, inputs cleared.
        assert_eq!(
            c.formulas,
            [(73, 147, "C7H20NSi2".into()), (200, 210, "C6H12O2".into())]
        );
        assert_eq!(c.mid_sizes, [(73, 147, 4), (200, 210, 6)]);
        assert_eq!(h.state().1.formula, "");

        let dir = std::env::temp_dir().join(format!("md-frag-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("lib.mdlib");
        md_core::library::save_mdlib(&path, std::slice::from_ref(&c)).unwrap();
        let back = md_core::library::read_mdlib(&path).unwrap();
        std::fs::remove_dir_all(&dir).unwrap();
        assert_eq!(back[0].formulas, c.formulas);
        assert_eq!(back[0].mid_sizes, c.mid_sizes);

        // The first "×" removes the first fragment and its MID size.
        h.get_all_by_label("×").next().unwrap().click();
        h.run();
        let c = &h.state().0;
        assert_eq!(c.formulas, [(200, 210, "C6H12O2".into())]);
        assert_eq!(c.mid_sizes, [(200, 210, 6)]);
    }

    #[test]
    fn a_repeated_start_and_end_cannot_be_added() {
        use egui_kittest::Harness;
        use egui_kittest::kittest::Queryable;
        let mut c = LibraryCompound::new("x");
        c.formulas.push((73, 147, "C3H9Si".into()));
        let new = NewFragment {
            begin: 73,
            end: 147,
            formula: "C2".into(),
            size: 0,
        };
        let mut h = Harness::new_ui_state(
            |ui, s: &mut (LibraryCompound, NewFragment)| {
                mid_fragments(ui, &mut s.0, &mut s.1);
            },
            (c, new),
        );
        h.run();
        h.get_by_label("Add fragment").click();
        h.run();
        assert_eq!(h.state().0.formulas.len(), 1);
    }

    #[test]
    fn automatic_ions_become_quant_ions_or_are_cleared() {
        use egui_kittest::Harness;
        use egui_kittest::kittest::Queryable;
        let mut c = LibraryCompound::new("x");
        c.quant_ions = vec![73.0, 147.0];
        c.auto_quant_ions = vec![116.0, 147.0];
        let state = (c, String::new());
        let mut h = Harness::new_ui_state(
            |ui, s: &mut (LibraryCompound, String)| {
                auto_ions(ui, &mut s.0, &mut s.1);
            },
            state,
        );
        h.run();
        h.get_by_label("Use as quant ions").click();
        h.run();
        assert_eq!(h.state().0.quant_ions, [73.0, 116.0, 147.0]);
        assert_eq!(h.state().1, "73 116 147");
        h.get_by_label("Clear automatic ions").click();
        h.run();
        assert!(h.state().0.auto_quant_ions.is_empty());
    }

    #[test]
    fn mid_sizes_follow_their_fragment() {
        let mut c = LibraryCompound::new("x");
        c.formulas.push((73, 147, "C3H9Si".into()));
        assert_eq!(mid_size(&c, 73, 147), 0);
        set_mid_size(&mut c, 73, 147, 2);
        set_mid_size(&mut c, 73, 147, 3);
        assert_eq!(c.mid_sizes, [(73, 147, 3)]);
        assert!(has_fragment(&c, 73, 147) && !has_fragment(&c, 73, 148));
    }
}
