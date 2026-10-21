//! Chromatogram view options (derivative, time unit, auto Y, score filter) and the run info window.

use crate::lang::{t, tf};
use crate::layout::ToolWindow;
use eframe::egui;
use md_core::TimeUnit;

use crate::identify::score_color;
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
    /// Draw a crosshair at the pointer in the chromatogram plots.
    pub crosshair: bool,
    /// Toolbar press for the chromatogram plots, applied in the frame it was made.
    pub nav: Option<Nav>,
    /// The same for the spectrum plot.
    pub spectrum_nav: Option<Nav>,
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
            crosshair: true,
            nav: None,
            spectrum_nav: None,
        }
    }
}

/// A toolbar button that moves or scales the visible x range of a plot.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Nav {
    In,
    Out,
    /// Back to the whole data; the plot does this itself through `Plot::reset`.
    All,
    Left,
    Right,
}

impl Nav {
    /// Zooms around the middle of the visible range or moves it by a quarter of its width.
    pub fn apply(self, p: &mut egui_plot::PlotUi) {
        let b = p.plot_bounds();
        let (lo, hi) = (b.min(), b.max());
        let mid = egui_plot::PlotPoint::new((lo[0] + hi[0]) / 2.0, (lo[1] + hi[1]) / 2.0);
        let width = (hi[0] - lo[0]) as f32;
        match self {
            Nav::In => p.zoom_bounds(egui::vec2(2.0, 1.0), mid),
            Nav::Out => p.zoom_bounds(egui::vec2(0.5, 1.0), mid),
            Nav::Left => p.translate_bounds(egui::vec2(-width / 4.0, 0.0)),
            Nav::Right => p.translate_bounds(egui::vec2(width / 4.0, 0.0)),
            Nav::All => {}
        }
    }
}

/// Zoom and move buttons and the "?" with the controls of a plot; `help` is one line per control.
pub fn nav_buttons(ui: &mut egui::Ui, help: &[&str]) -> Option<Nav> {
    let mut nav = None;
    for (n, text, hover) in [
        (
            Nav::In,
            "+",
            t(
                "Zoom in on the middle of the plot",
                "In die Mitte des Plots zoomen",
            ),
        ),
        (Nav::Out, "−", t("Zoom out", "Herauszoomen")),
        (
            Nav::All,
            t("All", "Alles"),
            t(
                "Show all data again (double-click on the plot does the same)",
                "Wieder alle Daten zeigen (ein Doppelklick auf den Plot tut dasselbe)",
            ),
        ),
        (
            Nav::Left,
            "◀",
            t(
                "Move the view to the left",
                "Die Ansicht nach links verschieben",
            ),
        ),
        (
            Nav::Right,
            "▶",
            t(
                "Move the view to the right",
                "Die Ansicht nach rechts verschieben",
            ),
        ),
    ] {
        if ui.button(text).on_hover_text(hover).clicked() {
            nav = Some(n);
        }
    }
    let r = ui
        .button("?")
        .on_hover_text(t("How to use this plot", "So bedienen Sie diesen Plot"));
    egui::Popup::from_toggle_button_response(&r)
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
        .show(|ui| {
            ui.set_max_width(340.0);
            for line in help {
                ui.label(*line);
            }
        });
    if nav.is_some() {
        // The decimation of the next frame needs the new range.
        ui.ctx().request_repaint();
    }
    nav
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

/// How long work took, the same in every status text: `80 ms`, `12.3 s`.
pub fn took(d: std::time::Duration) -> String {
    match d.as_millis() {
        ms @ 0..1000 => format!("{ms} ms"),
        _ => format!("{:.1} s", d.as_secs_f64()),
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
        ui.checkbox(&mut self.derivative, t("1st derivative", "1. Ableitung"))
            .on_hover_text(t(
                "Plot the slope per scan instead of the intensity; hides component marks",
                "Die Steigung pro Scan statt der Intensität zeichnen; blendet die Komponentenmarken aus",
            ));
        ui.checkbox(&mut self.auto_y, t("Auto Y", "Auto-Y"))
            .on_hover_text(t(
                "Fit the Y axis to the visible time range",
                "Die Y-Achse an den sichtbaren Zeitbereich anpassen",
            ));
        let before = self.unit;
        let time = ui.label(t("Time:", "Zeit:"));
        egui::ComboBox::from_id_salt("time unit")
            .selected_text(self.unit.label())
            .show_ui(ui, |ui| {
                for u in TimeUnit::ALL {
                    ui.selectable_value(&mut self.unit, u, u.label());
                }
            })
            .response
            .labelled_by(time.id)
            .on_hover_text(t("Unit of the time axis", "Einheit der Zeitachse"));
        if self.unit != before {
            self.bounds_in = Some(before);
        }
        ui.toggle_value(&mut self.crosshair, t("Crosshair", "Fadenkreuz"))
            .on_hover_text(t(
                "Draw cross lines at the pointer in the chromatograms",
                "Fadenkreuz am Mauszeiger in den Chromatogrammen zeichnen",
            ));
        ui.add(egui::Label::new(t("Names with score", "Namen mit Score")).extend());
        // The triangles wear the colours of the names they switch.
        let weak = ui.visuals().weak_text_color();
        for ((on, text), shade) in self
            .score_filter
            .iter_mut()
            .zip(["≤ 0.7", "> 0.7", "> 0.9"])
            .zip([0.5, 0.8, 0.95])
        {
            let color = if *on { score_color(shade) } else { weak };
            let label = egui::RichText::new(format!("▲ {text}")).color(color);
            if ui
                .button(label)
                .on_hover_text(t(
                    "A name is shown if any switch matching its library score is on",
                    "Ein Name wird gezeigt, wenn ein zu seinem Bibliotheks-Score passender Schalter an ist",
                ))
                .clicked()
            {
                *on = !*on;
            }
        }
    }
}

impl App {
    pub fn info_window(&mut self, ctx: &egui::Context) {
        let Some(run) = self.runs.get(self.active) else {
            return;
        };
        let mut open = self.view.info_open;
        egui::Window::new(t("Run info", "Laufinfo"))
            .id(egui::Id::unique("Run info"))
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
        let yes = |b: bool| if b { t("yes", "ja") } else { t("no", "nein") };
        let rows = [
            (t("Source file", "Quelldatei"), r.path.display().to_string()),
            (t("Scans", "Scans"), d.scan_count().to_string()),
            (
                t("First / last RT", "Erste / letzte RT"),
                format!(
                    "{:.3} / {:.3} min",
                    times.first().unwrap_or(&0.0) / 60.0,
                    times.last().unwrap_or(&0.0) / 60.0
                ),
            ),
            (
                t("Measurement time", "Messdauer"),
                format!(
                    "{:.1} s",
                    times.last().unwrap_or(&0.0) - times.first().unwrap_or(&0.0)
                ),
            ),
            (
                t("Lowest mass", "Kleinste Masse"),
                if lowest.is_finite() {
                    format!("{lowest:.2}")
                } else {
                    "–".into()
                },
            ),
            (
                t("Mass traces (nominal)", "Massenspuren (nominal)"),
                format!("{}", (highest.round() - lowest.round() + 1.0).max(0.0)),
            ),
            (
                t("Noise factor", "Rauschfaktor"),
                format!("{}", md_core::peaks::NOISE_FACTOR),
            ),
            (
                t("RI calibration", "RI-Kalibrierung"),
                yes(self.ri.cal.is_some()).into(),
            ),
            (
                t("RT correction", "RT-Korrektur"),
                yes(r.ident.as_ref().is_some_and(|i| i.2.is_some())).into(),
            ),
            (
                t("Peak detection", "Peakerkennung"),
                // The names of the Settings dialog.
                tf!(
                    "start slope {}, end slope {}, minimum S/N {}, {} apex steps per scan",
                    "Steigung am Peakbeginn {}, Steigung am Peakende {}, Minimum-S/N {}, {} Apex-Schritte pro Scan",
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
                // A grid column starts as narrow as its first frame; without a floor a long path
                // wraps into a one-word-wide column.
                ui.scope(|ui| {
                    ui.set_min_width(240.0);
                    ui.add(egui::Label::new(v).wrap());
                });
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
    fn durations_read_as_ms_or_s() {
        use std::time::Duration;
        assert_eq!(took(Duration::from_micros(80_460)), "80 ms");
        assert_eq!(took(Duration::from_millis(999)), "999 ms");
        assert_eq!(took(Duration::from_millis(12_340)), "12.3 s");
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
