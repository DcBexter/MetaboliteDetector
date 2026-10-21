//! Bottom status bar: progress, last error or the state of the active run.

use eframe::egui;

use crate::errors;
use crate::lang::{t, tf};
use crate::{App, Run};

impl App {
    pub fn status_bar(&mut self, ui: &mut egui::Ui) {
        egui::Panel::bottom("status").show(ui, |ui| {
            if !self.loading.is_empty() {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(tf!(
                        "Loading {} file(s)…",
                        "Lade {} Datei(en)…",
                        self.loading.len()
                    ));
                });
            } else if self.runs.iter().any(|r| r.analysis.busy()) {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(t(
                        "Re-analysing with the new settings…",
                        "Neu analysieren mit den neuen Einstellungen…",
                    ));
                });
            } else if let Some(e) = self.error.as_deref().map(errors::text) {
                ui.horizontal(|ui| {
                    if ui
                        .small_button("×")
                        .on_hover_text(t("Dismiss", "Ausblenden"))
                        .clicked()
                    {
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
                let mut status = self.run_status(r);
                let note = self.hidden_list_note().map(errors::text);
                if let Some(n) = &note {
                    status = format!("{n} — {status}");
                }
                ui.horizontal(|ui| {
                    if note.is_some()
                        && ui
                            .small_button("×")
                            .on_hover_text(t("Dismiss", "Ausblenden"))
                            .clicked()
                    {
                        self.results.set_note(String::new());
                    }
                    let label = egui::Label::new(&status)
                        .truncate()
                        .show_tooltip_when_elided(false);
                    ui.add(label)
                        .on_hover_text(format!("{status}\n{}", r.path.display()));
                });
            } else {
                ui.label(t("No runs open.", "Keine Läufe geöffnet."));
            }
        });
    }

    /// Result of Ctrl+S / Ctrl+E (or a button) while the component list, which shows it, is hidden.
    fn hidden_list_note(&self) -> Option<&str> {
        Some(self.results.note()).filter(|n| !self.panels.components && !n.is_empty())
    }

    fn run_status(&self, r: &Run) -> String {
        let shown = match r.component {
            Some(k) => tf!(
                "component {} of {} ({} ions)",
                "Komponente {} von {} ({} Ionen)",
                k + 1,
                r.components.len(),
                r.components[k].spectrum.len()
            ),
            None => tf!(
                "{} components (click a dot)",
                "{} Komponenten (einen Punkt anklicken)",
                r.components.len()
            ),
        };
        tf!(
            "{} — {} scans — scan {} at {:.3} {}{} (←/→ scan, ↑/↓ component) — {shown} — {}",
            "{} — {} Scans — Scan {} bei {:.3} {}{} (←/→ Scan, ↑/↓ Komponente) — {shown} — {}",
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

    /// The note must stand in the drawn bar, not only in the chooser (#326).
    #[test]
    fn the_drawn_status_bar_shows_the_note_while_the_list_is_hidden() {
        use egui_kittest::Harness;
        use egui_kittest::kittest::Queryable;
        let mut app = App::default();
        app.runs.push(crate::results::tests::run());
        app.results.set_note("Save failed: disk full".into());
        app.panels.components = false;
        let mut h = Harness::builder()
            .with_context_setup(crate::style::apply)
            .build_ui_state(|ui, app: &mut App| app.status_bar(ui), app);
        h.run();
        h.get_by_label_contains("Save failed: disk full");
        h.state_mut().panels.components = true;
        h.run();
        assert!(
            h.query_by_label_contains("Save failed: disk full")
                .is_none()
        );
    }
}
