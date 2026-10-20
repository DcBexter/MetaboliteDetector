//! Bottom status bar: progress, last error or the state of the active run.

use eframe::egui;

use crate::{App, Run};

impl App {
    pub fn status_bar(&self, ui: &mut egui::Ui) {
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
            } else if let Some(e) = &self.error {
                ui.colored_label(ui.visuals().error_fg_color, e);
            } else if let Some(r) = self.runs.get(self.active) {
                ui.label(self.run_status(r));
            } else {
                ui.label("Open one or more chromatograms (NetCDF or old .bin) via File → Open.");
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
            r.path.display(),
            r.data.scan_count(),
            r.scan + 1,
            self.view.time(r.minutes[r.scan]),
            self.view.unit.label(),
            self.ri.ri_note(r, r.scan),
            self.load_note
        )
    }
}
