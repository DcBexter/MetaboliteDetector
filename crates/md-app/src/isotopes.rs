//! Natural isotope calculator window: formula in, M+0, M+1, … out.

use crate::errors;
use crate::lang::{t, tf};
use crate::layout::ToolWindow;
use eframe::egui;
use md_core::isotope::{IsotopeTable, parse_formula};

pub struct Isotopes {
    pub open: bool,
    formula: String,
    /// Number of ions shown, `M+0` … `M+len-1`.
    len: usize,
    table: IsotopeTable,
}

impl Default for Isotopes {
    fn default() -> Self {
        Self {
            open: false,
            formula: "C6H12O6".into(),
            len: 6,
            table: IsotopeTable::standard(),
        }
    }
}

impl Isotopes {
    pub fn window(&mut self, ctx: &egui::Context) {
        let mut open = self.open;
        egui::Window::new(t("Natural isotopes", "Natürliche Isotope"))
            .id(egui::Id::unique("Natural isotopes"))
            .tool_window(ctx, 5, [340.0, 420.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui));
        self.open = open;
    }

    fn contents(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            let label = ui.label(t("Formula", "Summenformel"));
            ui.text_edit_singleline(&mut self.formula)
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Molecular formula, e.g. C6H12O6 or C3H9NSi (element symbols with counts)",
                    "Summenformel, z. B. C6H12O6 oder C3H9NSi (Elementsymbole mit Anzahl)",
                ));
            let label = ui.label(t("Ions", "Ionen"));
            ui.add(egui::DragValue::new(&mut self.len).range(1..=50))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "How many ions are shown: M+0 up to M+(n−1)",
                    "Wie viele Ionen angezeigt werden: M+0 bis M+(n−1)",
                ));
        });
        match parse_formula(self.formula.trim()).and_then(|f| self.table.pattern(&f, self.len)) {
            Err(e) => {
                ui.colored_label(ui.visuals().error_fg_color, errors::text(&e));
            }
            Ok(p) => {
                ui.label(tf!(
                    "Nominal mass {}, exact mass {:.4}",
                    "Nominale Masse {}, exakte Masse {:.4}",
                    p.nominal_mass,
                    p.monoisotopic_mass
                ));
                egui::Grid::new("isotope pattern")
                    .striped(true)
                    .show(ui, |ui| {
                        for (i, h) in [
                            t("Ion", "Ion"),
                            "m/z",
                            t("Abundance (%)", "Häufigkeit (%)"),
                            t("Relative to largest (%)", "Relativ zum größten (%)"),
                        ]
                        .into_iter()
                        .enumerate()
                        {
                            // Numbers are right-aligned, so their headers are too.
                            match i {
                                0 => {
                                    ui.strong(h);
                                }
                                _ => right(ui, |ui| ui.strong(h)),
                            }
                        }
                        ui.end_row();
                        let max = p.abundance.iter().copied().fold(0.0, f64::max);
                        for (i, a) in p.abundance.iter().enumerate() {
                            ui.label(format!("M+{i}"));
                            right(ui, |ui| ui.label(format!("{}", p.nominal_mass + i as u64)));
                            right(ui, |ui| ui.label(format!("{:.4}", 100.0 * a)));
                            let rel = if max > 0.0 { 100.0 * a / max } else { 0.0 };
                            right(ui, |ui| ui.label(format!("{rel:.2}")));
                            ui.end_row();
                        }
                    });
            }
        }
    }
}

/// A table cell with its content flush right, so that the decimal points line up.
fn right(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> egui::Response) {
    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
        add(ui);
    });
}
