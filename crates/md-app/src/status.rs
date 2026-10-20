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
            } else if let Some(n) = self.hidden_list_note().map(str::to_owned) {
                ui.horizontal(|ui| {
                    if ui.small_button("×").on_hover_text("Dismiss").clicked() {
                        self.results.set_note(String::new());
                    }
                    ui.add(
                        egui::Label::new(&n)
                            .truncate()
                            .show_tooltip_when_elided(false),
                    )
                    .on_hover_text(n);
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

    /// Result of Ctrl+S / Ctrl+E (or a button) while the component list, which shows it, is hidden.
    fn hidden_list_note(&self) -> Option<&str> {
        Some(self.results.note()).filter(|n| !self.panels.components && !n.is_empty())
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
            "{} — {} scans — scan {} at {:.3} {}{} (←/→ scan, ↑/↓ component) — {shown} — {}",
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn result_note_reaches_the_status_bar_only_while_the_list_is_hidden() {
        let mut app = App::default();
        app.results.set_note("Save failed: disk full".into());
        assert_eq!(app.hidden_list_note(), None);
        app.panels.components = false;
        assert_eq!(app.hidden_list_note(), Some("Save failed: disk full"));
    }
}
