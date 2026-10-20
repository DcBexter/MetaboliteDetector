//! Chromatogram view options (derivative, time unit, auto Y, score filter) and the run info window.

use eframe::egui;
use md_core::TimeUnit;

use crate::{App, Run};

pub struct View {
    pub unit: TimeUnit,
    pub derivative: bool,
    pub auto_y: bool,
    /// Show names with score ≤ 0.7, > 0.7 and > 0.9 (old label switches); a name shows if any matching switch is on.
    pub score_filter: [bool; 3],
    pub info_open: bool,
}

impl Default for View {
    fn default() -> Self {
        Self {
            unit: TimeUnit::default(),
            derivative: false,
            auto_y: true,
            score_filter: [true; 3],
            info_open: false,
        }
    }
}

impl View {
    /// Whether the old label switches let a name with this score through.
    pub fn shows(&self, score: f64) -> bool {
        let [low, mid, high] = self.score_filter;
        (score <= 0.7 && low) || (score > 0.7 && mid) || (score > 0.9 && high)
    }

    /// Seconds of the run's minute-based times in the chosen unit.
    pub fn time(&self, minutes: f64) -> f64 {
        self.unit.from_seconds(minutes * 60.0)
    }

    pub fn minutes(&self, t: f64) -> f64 {
        self.unit.to_seconds(t) / 60.0
    }

    pub fn controls(&mut self, ui: &mut egui::Ui) {
        ui.checkbox(&mut self.derivative, "1st derivative")
            .on_hover_text(
                "Plot the slope per scan instead of the intensity; hides component marks",
            );
        ui.checkbox(&mut self.auto_y, "Auto Y")
            .on_hover_text("Fit the Y axis to the visible time range");
        egui::ComboBox::from_id_salt("time unit")
            .selected_text(format!("Time: {}", self.unit.label()))
            .show_ui(ui, |ui| {
                for u in TimeUnit::ALL {
                    ui.selectable_value(&mut self.unit, u, u.label());
                }
            });
        ui.label("Show names with score");
        for (on, text) in self
            .score_filter
            .iter_mut()
            .zip(["≤ 0.7", "> 0.7", "> 0.9"])
        {
            ui.checkbox(on, text)
                .on_hover_text("A name is shown if any switch matching its library score is on");
        }
    }
}

impl App {
    pub fn info_window(&mut self, ctx: &egui::Context) {
        let Some(run) = self.runs.get(self.active) else {
            return;
        };
        let mut open = self.view.info_open;
        egui::Window::new("Run info")
            .open(&mut open)
            .default_width(420.0)
            .show(ctx, |ui| self.run_info(ui, run));
        self.view.info_open = open;
    }

    fn run_info(&self, ui: &mut egui::Ui, r: &Run) {
        let d = &r.data;
        let times = &d.scan_times_s;
        // Scans are m/z-ascending, so the extremes sit at the scan ends.
        let first = (0..d.scan_count()).filter_map(|s| d.spectrum(s).0.first().copied());
        let last = (0..d.scan_count()).filter_map(|s| d.spectrum(s).0.last().copied());
        let lowest = first.fold(f32::INFINITY, f32::min);
        let highest = last.fold(f32::NEG_INFINITY, f32::max);
        let peak = self.settings.key().0;
        let yes = |b: bool| if b { "yes" } else { "no" };
        let rows = [
            ("Source file", r.path.display().to_string()),
            ("Scans", d.scan_count().to_string()),
            (
                "First / last RT",
                format!(
                    "{:.3} / {:.3} min",
                    times.first().unwrap_or(&0.0) / 60.0,
                    times.last().unwrap_or(&0.0) / 60.0
                ),
            ),
            (
                "Measurement time",
                format!(
                    "{:.1} s",
                    times.last().unwrap_or(&0.0) - times.first().unwrap_or(&0.0)
                ),
            ),
            ("Lowest mass", format!("{lowest:.2}")),
            (
                "Mass traces (nominal)",
                format!("{}", (highest.round() - lowest.round() + 1.0).max(0.0)),
            ),
            ("Noise factor", format!("{}", md_core::peaks::NOISE_FACTOR)),
            ("RI calibration", yes(self.ri.cal.is_some()).into()),
            (
                "RT correction",
                yes(r.ident.as_ref().is_some_and(|i| i.2.is_some())).into(),
            ),
            (
                "Peak detection",
                format!(
                    "slope begin {} / end {}, min height {} noise units, {} bins",
                    peak.threshold_begin,
                    peak.threshold_end,
                    peak.min_peak_height,
                    peak.number_bins
                ),
            ),
        ];
        egui::Grid::new("run info").striped(true).show(ui, |ui| {
            for (k, v) in rows {
                ui.strong(k);
                ui.add(egui::Label::new(v).wrap());
                ui.end_row();
            }
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn score_switches_pick_their_bands() {
        let only = |score_filter| View {
            score_filter,
            ..View::default()
        };
        let high = only([false, false, true]);
        assert!(high.shows(0.95) && !high.shows(0.8) && !high.shows(0.5));
        let low = only([true, false, false]);
        assert!(low.shows(0.7) && !low.shows(0.71));
    }
}
