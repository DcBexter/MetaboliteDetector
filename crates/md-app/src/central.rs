//! Central panel: trace selection, TIC/SIC plots, scan cursor and the spectrum below them.

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, Points, Text, VLine};
use md_core::decimate;

use crate::{App, Run, Trace, color, export, identify, view::View};

/// One chromatogram plot with a line per run and `marks` as dots; returns the clicked retention
/// time in the view's unit.
#[allow(clippy::too_many_arguments)]
fn trace_plot(
    ui: &mut egui::Ui,
    id: &str,
    label: &str,
    height: f32,
    runs: &[Run],
    view: &View,
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
        .x_axis_label(format!("Retention time ({})", view.unit.label()))
        .y_axis_label(label)
        .show(ui, |p| {
            // Draw at most ~2 points per pixel of the visible range, whatever the file size.
            let (lo, hi) = if p.auto_bounds().x {
                (f64::NEG_INFINITY, f64::INFINITY)
            } else {
                (
                    view.minutes(p.plot_bounds().min()[0]),
                    view.minutes(p.plot_bounds().max()[0]),
                )
            };
            let mut y_span: Option<(f64, f64)> = None;
            let width = p.response().rect.width().max(100.0) as usize;
            for (i, run) in runs.iter().enumerate() {
                let ys = trace(run);
                // A SIC is empty until its first job finishes.
                if run.visible && ys.len() == run.minutes.len() {
                    // shortcut: the derivative is recomputed every frame (O(scans)); cache per run if runs grow past ~1M scans.
                    let d;
                    let ys = if view.derivative {
                        d = md_core::first_derivative(ys);
                        &d
                    } else {
                        ys
                    };
                    // Decimate in minutes, then convert: the scan arrays stay untouched.
                    let mut points = decimate(&run.minutes, ys, lo, hi, width);
                    for pt in &mut points {
                        pt[0] = view.time(pt[0]);
                    }
                    for (a, b) in points.iter().map(|q| (q[1], q[1])) {
                        y_span = Some(y_span.map_or((a, b), |(l, h)| (l.min(a), h.max(b))));
                    }
                    p.line(Line::new(run.name.as_str(), points).color(color(i)));
                }
            }
            if view.auto_y
                && let Some((l, h)) = y_span
            {
                let pad = ((h - l) * 0.05).max(f64::EPSILON);
                p.set_plot_bounds_y(l - pad..=h + pad);
            }
            if !marks.is_empty() {
                let marks: Vec<[f64; 2]> = marks.iter().map(|m| [view.time(m[0]), m[1]]).collect();
                p.points(Points::new("components", marks).radius(3.0));
            }
            for (at, name, color) in labels {
                let text = egui::RichText::new(name).color(*color);
                let at = [view.time(at[0]), at[1]];
                p.text(Text::new("hits", at.into(), text).anchor(egui::Align2::LEFT_BOTTOM));
            }
            p.vline(VLine::new("", view.time(cursor_min)));
            p.pointer_coordinate()
        });
    plot.response
        .clicked()
        .then_some(view.minutes(plot.inner?.x))
}

impl App {
    pub fn central_panel(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default().show(ui, |ui| {
            self.trace_controls(ui);
            ui.horizontal_wrapped(|ui| {
                self.view.controls(ui);
                if ui.button("Run info").clicked() {
                    self.view.info_open = true;
                }
            });
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

    fn derivative_label(&self, label: String) -> String {
        if self.view.derivative {
            format!("d/dscan {label}")
        } else {
            label
        }
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
            // Components of the active run as dots on its TIC; not on the derivative, whose y differs.
            let marks: Vec<[f64; 2]> = run
                .components
                .iter()
                .filter(|_| !self.view.derivative)
                .map(|k| {
                    [
                        run.minutes[k.apex_scan as usize],
                        run.tic[k.apex_scan as usize],
                    ]
                })
                .collect();
            let labels = if self.view.derivative {
                Vec::new()
            } else {
                identify::labels(run, |s| self.view.shows(s))
            };
            clicked = trace_plot(
                ui,
                "tic",
                &self.derivative_label("TIC".into()),
                height,
                &self.runs,
                &self.view,
                cursor,
                |r| &r.tic,
                &marks,
                &labels,
            );
        }
        if self.trace != Trace::Tic {
            let label = format!("SIC m/z {} ± {}", self.mz, self.tolerance);
            let label = self.derivative_label(label);
            let sic = trace_plot(
                ui,
                "sic",
                &label,
                height,
                &self.runs,
                &self.view,
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
