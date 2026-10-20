//! Left panel: the open runs with visibility, colour, selection and close buttons.

use eframe::egui;

use crate::{App, color};

impl App {
    pub fn run_list(&mut self, ui: &mut egui::Ui) {
        if self.runs.is_empty() {
            return;
        }
        crate::layout::side_panel(ui, egui::Panel::left, "runs", |ui| {
            let mut close = None;
            for (i, run) in self.runs.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    let visible = ui
                        .checkbox(&mut run.visible, "")
                        .on_hover_text("Show in plots");
                    ui.colored_label(color(i), "■");
                    let name = ui
                        .selectable_label(i == self.active, &run.name)
                        .on_hover_text(run.path.display().to_string());
                    if name.clicked() {
                        self.active = i;
                    }
                    visible.labelled_by(name.id);
                    // At the row end, away from the name, so selecting a run never closes one.
                    if ui.small_button("×").on_hover_text("Close").clicked() {
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
        });
    }
}
