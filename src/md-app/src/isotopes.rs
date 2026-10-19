//! Natural isotope calculator window: formula in, M+0, M+1, … out.

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
        egui::Window::new("Natural isotopes")
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui));
        self.open = open;
    }

    fn contents(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.label("Formula");
            ui.text_edit_singleline(&mut self.formula);
            ui.label("Ions");
            ui.add(egui::DragValue::new(&mut self.len).range(1..=50));
        });
        match parse_formula(self.formula.trim()).and_then(|f| self.table.pattern(&f, self.len)) {
            Err(e) => {
                ui.colored_label(ui.visuals().error_fg_color, e);
            }
            Ok(p) => {
                ui.label(format!(
                    "Nominal mass {}, exact mass {:.4}",
                    p.nominal_mass, p.monoisotopic_mass
                ));
                egui::Grid::new("isotope pattern")
                    .striped(true)
                    .show(ui, |ui| {
                        for h in ["Ion", "m/z", "Abundance (%)", "Relative to largest (%)"] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        let max = p.abundance.iter().copied().fold(0.0, f64::max);
                        for (i, a) in p.abundance.iter().enumerate() {
                            ui.label(format!("M+{i}"));
                            ui.label(format!("{}", p.nominal_mass + i as u64));
                            ui.label(format!("{:.4}", 100.0 * a));
                            ui.label(format!(
                                "{:.2}",
                                if max > 0.0 { 100.0 * a / max } else { 0.0 }
                            ));
                            ui.end_row();
                        }
                    });
            }
        }
    }
}
