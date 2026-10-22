//! Left panel: the open runs with visibility, colour, selection and close buttons, and the component
//! list of the active run below them.

use eframe::egui;

use crate::lang::t;
use crate::{App, color, layout};

impl App {
    pub fn left_panel(&mut self, ui: &mut egui::Ui) {
        let p = self.panels;
        if self.runs.is_empty() || !(p.runs || p.components) {
            return;
        }
        layout::dock(egui::Panel::left("left"), ui, 300.0).show(ui, |ui| {
            if p.runs && p.components {
                let h = ui.available_height();
                egui::Panel::top("runs")
                    .resizable(true)
                    .frame(egui::Frame::NONE.inner_margin(egui::Margin {
                        bottom: 6,
                        ..Default::default()
                    }))
                    .default_size(h * 0.3)
                    .size_range(60.0_f32.min(h * 0.3)..=h * 0.7)
                    .show(ui, |ui| self.run_list(ui));
                // The splitter's grab zone reaches 3 px into the list below.
                ui.add_space(6.0);
            } else if p.runs {
                self.run_list(ui);
            }
            if p.components {
                let cal = (self.runs.get(self.active)).and_then(|r| self.ri.for_run(r).cloned());
                let cutoff = self.settings.applied().ls_cutoff;
                (self.results).panel(ui, self.runs.get_mut(self.active), cal.as_ref(), cutoff);
            }
        });
    }

    fn run_list(&mut self, ui: &mut egui::Ui) {
        ui.strong(t("Runs", "Läufe"));
        egui::ScrollArea::both()
            .id_salt("run list")
            .auto_shrink([false, true])
            .show(ui, |ui| self.run_rows(ui));
    }

    fn run_rows(&mut self, ui: &mut egui::Ui) {
        let mut close = None;
        for (i, run) in self.runs.iter_mut().enumerate() {
            ui.horizontal(|ui| {
                let visible = ui
                    .checkbox(&mut run.visible, "")
                    .on_hover_text(t("Show in plots", "In den Plots zeigen"));
                ui.colored_label(color(i), "■");
                let name = ui
                    .selectable_label(i == self.active, &run.name)
                    .on_hover_text(run.path.display().to_string());
                if name.clicked() {
                    self.active = i;
                }
                visible.labelled_by(name.id);
                // At the row end, away from the name, so selecting a run never closes one.
                if ui
                    .small_button("×")
                    .on_hover_text(t("Close", "Schließen"))
                    .clicked()
                {
                    close = Some(i);
                }
            });
        }
        if let Some(i) = close {
            self.runs.remove(i);
            if self.active > i || self.active >= self.runs.len() {
                self.active = self.active.saturating_sub(1);
            }
        }
    }
}
