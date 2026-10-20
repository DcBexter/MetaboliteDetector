//! Bottom status bar: progress, last error or the state of the active run.

use eframe::egui;

use crate::{App, Run};

impl App {
    pub fn status_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("status").show(ui, |ui| {
            if !self.loading.is_empty() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Loading {} file(s)…", self.loading.len()));
                });
            } else if self.runs.iter().any(|r| r.analysis.busy()) {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Re-analysing with the new settings…");
                });
            } else if let Some(e) = self.error.clone() {
                ui.horizontal(|ui| {
                    if ui.small_button("×").on_hover_text("Dismiss").clicked() {
                        self.error = None;
                    }
                    let text = egui::RichText::new(&e).color(ui.visuals().error_fg_color);
                    let label = egui::Label::new(text)
                        .truncate()
                        .show_tooltip_when_elided(false);
                    ui.add(label).on_hover_text(e);
                });
            } else if let Some(r) = self.runs.get(self.active) {
                // One line: a wrapped status changes the bar's height, which egui applies a frame
                // late, so the second line would hang below the window for that frame (#200).
                let status = self.run_status(r);
                let label = egui::Label::new(&status)
                    .truncate()
                    .show_tooltip_when_elided(false);
                ui.add(label)
                    .on_hover_text(format!("{status}\n{}", r.path.display()));
            } else {
                ui.label("No runs open.");
            }
        });
    }

    fn run_status(&self, r: &Run) -> String {
        let shown = match r.component {
            Some(k) => format!(
                "component {} of {} ({} ions)",
                k + 1,
                r.components.len(),
                r.components[k].spectrum.len()
            ),
            None => format!("{} components (click a dot)", r.components.len()),
        };
        format!(
            "{} — {} scans — scan {} at {:.3} {}{} (←/→ to step) — {shown} — {}",
            r.name,
            r.data.scan_count(),
            r.scan + 1,
            self.view.time(r.minutes[r.scan]),
            self.view.unit.label(),
            self.ri.ri_note(r, r.scan),
            self.load_note
        )
    }
}
