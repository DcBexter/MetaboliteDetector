//! Central panel: trace selection, TIC/SIC plots, scan cursor and the spectrum below them.

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, Points, Text, VLine};
use md_core::decimate;

use crate::{App, Run, Trace, color, export, identify};

/// One chromatogram plot with a line per run and `marks` as dots; returns the clicked retention
/// time (min).
#[allow(clippy::too_many_arguments)]
fn trace_plot(
    ui: &mut egui::Ui,
    id: &str,
    label: &str,
    height: f32,
    runs: &[Run],
    cursor_min: f64,
    trace: fn(&Run) -> &[f64],
    marks: &[[f64; 2]],
    labels: &[identify::Label],
) -> Option<f64> {
    let link = egui::Id::unique("chromatograms");
    let plot = Plot::new(id)
        .height(height)
        // Same width in both plots, so their linked x axes line up.
        .y_axis_min_width(90.0)
        .link_axis(link, [true, false])
        .link_cursor(link, [true, false])
        .x_axis_label("Retention time (min)")
        .y_axis_label(label)
        .show(ui, |p| {
            // Draw at most ~2 points per pixel of the visible range, whatever the file size.
            let (lo, hi) = if p.auto_bounds().x {
                (f64::NEG_INFINITY, f64::INFINITY)
            } else {
                (p.plot_bounds().min()[0], p.plot_bounds().max()[0])
            };
            let width = p.response().rect.width().max(100.0) as usize;
            for (i, run) in runs.iter().enumerate() {
                let ys = trace(run);
                // A SIC is empty until its first job finishes.
                if run.visible && ys.len() == run.minutes.len() {
                    let points = decimate(&run.minutes, ys, lo, hi, width);
                    p.line(Line::new(run.name.as_str(), points).color(color(i)));
                }
            }
            if !marks.is_empty() {
                p.points(Points::new("components", marks.to_vec()).radius(3.0));
            }
            for (at, name, color) in labels {
                let text = egui::RichText::new(name).color(*color);
                p.text(Text::new("hits", (*at).into(), text).anchor(egui::Align2::LEFT_BOTTOM));
            }
            p.vline(VLine::new("", cursor_min));
            p.pointer_coordinate()
        });
    plot.response.clicked().then_some(plot.inner?.x)
}

impl App {
    pub fn central_panel(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default().show(ui, |ui| {
            self.trace_controls(ui);
            if self.trace != Trace::Tic {
                self.update_sics(ui.ctx());
            }
            if self.runs.get(self.active).is_none() {
                return;
            }
            self.step_scan(ui);
            let clicked = self.chromatogram_plots(ui);
            self.select_clicked(clicked);
            self.spectrum_plot(ui);
        });
    }

    fn trace_controls(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.trace, Trace::Tic, "TIC");
            ui.selectable_value(&mut self.trace, Trace::Sic, "SIC");
            ui.selectable_value(&mut self.trace, Trace::Both, "TIC + SIC");
            ui.add_enabled_ui(self.trace != Trace::Tic, |ui| {
                ui.label("m/z");
                ui.add(
                    egui::DragValue::new(&mut self.mz)
                        .range(1.0..=2000.0)
                        .max_decimals(2),
                );
                ui.label("±");
                ui.add(
                    egui::DragValue::new(&mut self.tolerance)
                        .range(0.0..=10.0)
                        .speed(0.05)
                        .max_decimals(2),
                );
            });
        });
    }

    /// Arrow keys move the scan cursor of the active run.
    fn step_scan(&mut self, ui: &egui::Ui) {
        let run = &mut self.runs[self.active];
        let last = run.data.scan_count() - 1;
        ui.input(|i| {
            if i.key_pressed(egui::Key::ArrowRight) {
                (run.scan, run.component) = ((run.scan + 1).min(last), None);
            }
            if i.key_pressed(egui::Key::ArrowLeft) {
                (run.scan, run.component) = (run.scan.saturating_sub(1), None);
            }
        });
    }

    /// The TIC and/or SIC plots; returns the clicked retention time (min).
    fn chromatogram_plots(&mut self, ui: &mut egui::Ui) -> Option<f64> {
        let run = &self.runs[self.active];
        let cursor = run.minutes[run.scan];
        // Both traces: a third each, so their axes keep enough room for tick labels.
        let plots = if self.trace == Trace::Both { 3.0 } else { 2.0 };
        let height = ui.available_height() / plots;
        let mut clicked = None;
        let top = ui.cursor().min;
        if self.trace != Trace::Sic {
            // Components of the active run as dots on its TIC.
            let marks: Vec<[f64; 2]> = run
                .components
                .iter()
                .map(|k| {
                    [
                        run.minutes[k.apex_scan as usize],
                        run.tic[k.apex_scan as usize],
                    ]
                })
                .collect();
            let labels = identify::labels(run);
            clicked = trace_plot(
                ui,
                "tic",
                "TIC",
                height,
                &self.runs,
                cursor,
                |r| &r.tic,
                &marks,
                &labels,
            );
        }
        if self.trace != Trace::Tic {
            let label = format!("SIC m/z {} ± {}", self.mz, self.tolerance);
            let sic = trace_plot(
                ui,
                "sic",
                &label,
                height,
                &self.runs,
                cursor,
                |r| &r.sic,
                &[],
                &[],
            );
            clicked = clicked.or(sic);
        }
        self.export.chromatograms =
            egui::Rect::from_min_max(top, egui::pos2(ui.max_rect().right(), ui.cursor().min.y));
        clicked
    }

    /// A click selects the scan at that time, snapping to a component apex within 3 scans.
    fn select_clicked(&mut self, clicked: Option<f64>) {
        let run = &mut self.runs[self.active];
        if let Some(s) = clicked.and_then(|x| run.data.scan_at(x * 60.0)) {
            run.component = (0..run.components.len())
                .min_by_key(|k| run.components[*k].apex_scan.abs_diff(s as u32))
                .filter(|k| run.components[*k].apex_scan.abs_diff(s as u32) <= 3);
            run.scan = run
                .component
                .map_or(s, |k| run.components[k].apex_scan as usize);
        }
    }

    fn spectrum_plot(&mut self, ui: &mut egui::Ui) {
        let run = &self.runs[self.active];
        export::norm_picker(ui, &mut self.export.norm);
        let top = ui.cursor().min;
        let (name, mz, y) = export::shown_spectrum(run, self.export.norm);
        let bars = mz.iter().zip(&y).map(|(m, i)| Bar::new(*m, *i).width(0.3));
        let label = match self.export.norm {
            md_core::spectrum::Normalization::None => "Intensity",
            _ => "Normalised intensity",
        };
        Plot::new("spectrum")
            .x_axis_label("m/z")
            .y_axis_label(format!("{label} — {}", run.name))
            .show(ui, |p| p.bar_chart(BarChart::new(name, bars.collect())));
        // Not the plot response: that is only the inner area, the PNG needs axes and titles too.
        self.export.spectrum =
            egui::Rect::from_min_max(top, egui::pos2(ui.max_rect().right(), ui.cursor().min.y));
    }
}
