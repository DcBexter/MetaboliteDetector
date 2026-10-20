//! Chromatogram view options (derivative, time unit, auto Y, score filter) and the run info window.

use crate::layout::ToolWindow;
use eframe::egui;
use md_core::TimeUnit;

use crate::{App, Run};

pub struct View {
    pub unit: TimeUnit,
    /// Unit the zoomed plot bounds are still in for one frame after a change of `unit`.
    pub bounds_in: Option<TimeUnit>,
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
            bounds_in: None,
            derivative: false,
            auto_y: true,
            score_filter: [true; 3],
            info_open: false,
        }
    }
}

/// Large intensities compactly, the same in every chromatogram plot: `12 M`, `1.2 k`, `-3`.
/// `step` is the tick spacing, so neighbouring ticks never share a label.
pub fn compact(v: f64, step: f64) -> String {
    let units = [(1e9, " G"), (1e6, " M"), (1e3, " k")];
    // Values below 10 000 stay plain; a value that rounds up to 1000 of its unit moves to the next one.
    let mut unit = units
        .iter()
        .find(|(d, _)| v.abs() >= d * if *d == 1e3 { 10.0 } else { 1.0 });
    if let Some((d, _)) = unit
        && (v / d).abs() >= 999.5
        && *d < 1e9
    {
        unit = units.iter().rev().find(|(u, _)| *u == d * 1000.0);
    }
    let (div, suffix) = unit.copied().unwrap_or((1.0, ""));
    let decimals = (-(step.abs() / div).log10()).ceil().clamp(0.0, 6.0) as usize;
    format!("{:.*}{suffix}", decimals, v / div)
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
        let before = self.unit;
        let time = ui.label("Time:");
        egui::ComboBox::from_id_salt("time unit")
            .selected_text(self.unit.label())
            .show_ui(ui, |ui| {
                for u in TimeUnit::ALL {
                    ui.selectable_value(&mut self.unit, u, u.label());
                }
            })
            .response
            .labelled_by(time.id)
            .on_hover_text("Unit of the time axis");
        if self.unit != before {
            self.bounds_in = Some(before);
        }
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
            .tool_window(ctx, 7, [420.0, 420.0])
            .open(&mut open)
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
            (
                "Lowest mass",
                if lowest.is_finite() {
                    format!("{lowest:.2}")
                } else {
                    "–".into()
                },
            ),
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
    fn large_values_are_compact() {
        assert_eq!(compact(12_000_000.0, 2e6), "12 M");
        assert_eq!(compact(1_200_000.0, 2e5), "1.2 M");
        assert_eq!(compact(25_000.0, 5e3), "25 k");
        assert_eq!(compact(9_000.0, 1e3), "9000");
        assert_eq!(compact(-0.5, 0.1), "-0.5");
        // Fine ticks keep distinct labels; no "1000 k".
        assert_eq!(compact(1_220_000.0, 2e4), "1.22 M");
        assert_ne!(compact(1_200_000.0, 2e4), compact(1_220_000.0, 2e4));
        assert_eq!(compact(999_999.0, 1e5), "1.0 M");
    }

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
