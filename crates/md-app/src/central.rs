//! Central panel: trace selection, TIC/SIC plots, scan cursor and the spectrum below them.

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, Points, Text, VLine};

use crate::plot_export::Area;
use md_core::decimate;
use std::sync::Arc;

use crate::{App, Run, Trace, color, detail, export, identify, view::View};

/// Visible x range in minutes (unbounded when auto-fitted); after a unit change it also moves the
/// zoomed window into the new unit.
fn visible_minutes(p: &mut egui_plot::PlotUi, view: &View) -> (f64, f64) {
    if p.auto_bounds().x {
        return (f64::NEG_INFINITY, f64::INFINITY);
    }
    let b = p.plot_bounds();
    let old = view.bounds_in.unwrap_or(view.unit);
    let to_min = |x: f64| old.to_seconds(x) / 60.0;
    if old != view.unit {
        let conv = |x: f64| view.time(to_min(x));
        p.set_plot_bounds(egui_plot::PlotBounds::from_min_max(
            [conv(b.min()[0]), b.min()[1]],
            [conv(b.max()[0]), b.max()[1]],
        ));
    }
    (to_min(b.min()[0]), to_min(b.max()[0]))
}

/// Points of one run in the view's unit: at most ~2 per pixel, optionally as the derivative.
fn run_points(
    run: &Run,
    ys: &[f64],
    view: &View,
    (lo, hi): (f64, f64),
    width: usize,
) -> Vec<[f64; 2]> {
    // shortcut: the derivative is recomputed every frame (O(scans)); cache per run past ~1M scans.
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
    points
}

/// Rough width of one label character and the gap between labels, in points.
const LABEL_CHAR_PX: f32 = 7.0;
const LABEL_GAP_PX: f32 = 6.0;
/// A click this close to a component apex (along the time axis) selects the component.
const SNAP_PX: f64 = 10.0;

/// One chromatogram plot with a line per run and `marks` as dots, unzoomed when `reset`; returns
/// the clicked retention time and the snap radius `SNAP_PX`, both in minutes.
#[allow(clippy::too_many_arguments)]
fn trace_plot(
    ui: &mut egui::Ui,
    (id, reset): (&str, bool),
    label: &str,
    height: f32,
    runs: &[Run],
    view: &View,
    cursor_min: f64,
    trace: fn(&Run) -> &[f64],
    marks: &[[f64; 2]],
    labels: &[identify::Label],
) -> Option<(f64, f64)> {
    let link = egui::Id::unique("chromatograms");
    let mut plot = Plot::new(id);
    if reset {
        plot = plot.reset();
    }
    let plot = plot
        .height(height)
        // Same width in both plots, so their linked x axes line up.
        .y_axis_min_width(90.0)
        .link_axis(link, [true, false])
        .link_cursor(link, [true, false])
        .x_axis_label(format!("Retention time ({})", view.unit.label()))
        .y_axis_label(label)
        .y_axis_formatter(|m, _| crate::view::compact(m.value, m.step_size))
        .show(ui, |p| {
            // Draw at most ~2 points per pixel of the visible range, whatever the file size.
            let (lo, hi) = visible_minutes(p, view);
            let mut y_span: Option<(f64, f64)> = None;
            let width = p.response().rect.width().max(100.0) as usize;
            for (i, run) in runs.iter().enumerate() {
                let ys = trace(run);
                // A SIC is empty until its first job finishes.
                if run.visible && ys.len() == run.minutes.len() {
                    let points = run_points(run, ys, view, (lo, hi), width);
                    for q in &points {
                        let (l, h) = y_span.unwrap_or((q[1], q[1]));
                        y_span = Some((l.min(q[1]), h.max(q[1])));
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
            // Hide overlapping names, best score first (3.3 does the same).
            let b = p.plot_bounds();
            let px_per_t = p.response().rect.width() as f64 / (b.max()[0] - b.min()[0]).max(1e-9);
            let starts: Vec<f32> = labels
                .iter()
                .map(|l| ((view.time(l.0[0]) - b.min()[0]) * px_per_t) as f32)
                .collect();
            let widths: Vec<f32> = labels
                .iter()
                .map(|l| l.1.chars().count() as f32 * LABEL_CHAR_PX + LABEL_GAP_PX)
                .collect();
            let scores: Vec<f64> = labels.iter().map(|l| l.3).collect();
            let keep = identify::visible_labels(&starts, &widths, &scores);
            for ((at, name, color, _), _) in labels.iter().zip(&keep).filter(|(_, k)| **k) {
                let text = egui::RichText::new(name).color(*color);
                let at = [view.time(at[0]), at[1]];
                p.text(Text::new("hits", at.into(), text).anchor(egui::Align2::LEFT_BOTTOM));
            }
            p.vline(VLine::new("", view.time(cursor_min)));
            (p.pointer_coordinate(), view.minutes(SNAP_PX / px_per_t))
        });
    let (at, radius) = plot.inner;
    plot.response
        .clicked()
        .then_some((view.minutes(at?.x), radius))
}

impl App {
    pub fn central_panel(&mut self, ui: &mut egui::Ui) {
        egui::CentralPanel::default().show(ui, |ui| {
            if self.runs.is_empty() {
                return self.start_screen(ui);
            }
            if self.trace != Trace::Tic {
                self.update_sics(ui.ctx());
            }
            if self.runs.get(self.active).is_none() {
                return;
            }
            self.step_scan(ui);
            // The chromatogram on top, the spectrum and Detected Peaks side by side below; each
            // split can be dragged, and the top leaves the spectrum at least 160 px.
            let p = self.panels;
            // A hidden plot has nothing on screen to save as a picture.
            if !p.chromatogram {
                self.export.chromatograms = Area::default();
            }
            if !p.spectrum {
                self.export.spectrum = Area::default();
            }
            let below = p.spectrum || p.detected_peaks;
            let none = egui::Frame::NONE;
            if p.chromatogram && below {
                let h = ui.available_height();
                egui::Panel::top("chromatograms")
                    .resizable(true)
                    .frame(none.inner_margin(egui::Margin {
                        bottom: 6,
                        ..Default::default()
                    }))
                    .default_size(h * 0.55)
                    .size_range((h * 0.25)..=(h - 160.0).max(h * 0.5))
                    .show(ui, |ui| self.chromatogram_area(ui));
                // The splitter's grab zone reaches 3 px into the area below.
                ui.add_space(6.0);
            } else if p.chromatogram {
                self.chromatogram_area(ui);
            }
            if p.detected_peaks && p.spectrum {
                let w = ui.available_width();
                egui::Panel::right("detected peaks")
                    .resizable(true)
                    .frame(none.inner_margin(egui::Margin {
                        left: 6,
                        ..Default::default()
                    }))
                    .default_size(w * 0.4)
                    .size_range((w * 0.2)..=(w * 0.8))
                    .show(ui, |ui| self.detected_peaks(ui));
            } else if p.detected_peaks {
                self.detected_peaks(ui);
            }
            if p.spectrum {
                let right = if p.detected_peaks { 6 } else { 0 };
                (none.inner_margin(egui::Margin {
                    right,
                    ..Default::default()
                }))
                .show(ui, |ui| self.spectrum_plot(ui));
            }
        });
    }

    fn chromatogram_area(&mut self, ui: &mut egui::Ui) {
        self.trace_controls(ui);
        ui.horizontal_wrapped(|ui| {
            self.view.controls(ui);
            if ui.button("Run info").clicked() {
                self.view.info_open = true;
            }
        });
        let clicked = self.chromatogram_plots(ui);
        self.select_clicked(clicked);
    }

    /// Nothing open yet: the two ways in, instead of empty plots and controls. While the first runs
    /// load, a spinner; the status bar counts the files.
    fn start_screen(&mut self, ui: &mut egui::Ui) {
        ui.vertical_centered(|ui| {
            ui.add_space(ui.available_height() / 3.0);
            if !self.loading.is_empty() {
                ui.spinner();
                return;
            }
            ui.label(
                "Open GC-MS runs (NetCDF or old .bin) to see their chromatograms and spectra.",
            );
            ui.add_space(8.0);
            self.open_chromatograms_button(ui);
            self.open_library_button(ui);
        });
    }

    fn trace_controls(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut self.trace, Trace::Tic, "TIC")
                .on_hover_text("Total ion chromatogram: the sum of all masses per scan");
            ui.selectable_value(&mut self.trace, Trace::Sic, "SIC")
                .on_hover_text("Selected ion chromatogram: only the chosen m/z");
            ui.selectable_value(&mut self.trace, Trace::Both, "TIC + SIC")
                .on_hover_text("Both, one above the other with a shared time axis");
            // Each control on its own, not in a nested scope, so the row can wrap between them.
            let sic = self.trace != Trace::Tic;
            let label = ui.add_enabled(sic, egui::Label::new("m/z"));
            ui.add_enabled(
                sic,
                egui::DragValue::new(&mut self.mz)
                    .range(1.0..=2000.0)
                    .max_decimals(2),
            )
            .labelled_by(label.id)
            .on_hover_text("Ion whose intensity the SIC trace shows")
            .on_disabled_hover_text("Used by the SIC trace; pick SIC or TIC + SIC");
            let label = ui.add_enabled(sic, egui::Label::new("±"));
            ui.add_enabled(
                sic,
                egui::DragValue::new(&mut self.tolerance)
                    .range(0.0..=10.0)
                    .speed(0.05)
                    .max_decimals(2),
            )
            .labelled_by(label.id)
            .on_hover_text("Window around the m/z that the SIC trace adds up")
            .on_disabled_hover_text("Used by the SIC trace; pick SIC or TIC + SIC");
        });
    }

    pub fn derivative_label(&self, label: String) -> String {
        if self.view.derivative {
            format!("d/dscan {label}")
        } else {
            label
        }
    }

    /// Arrow keys move the scan cursor of the active run, unless a text field takes them.
    fn step_scan(&mut self, ui: &egui::Ui) {
        if ui.ctx().text_edit_focused() {
            return;
        }
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

    /// The TIC and/or SIC plots, filling the room left; returns the clicked retention time and the
    /// snap radius (min).
    fn chromatogram_plots(&mut self, ui: &mut egui::Ui) -> Option<(f64, f64)> {
        let run = &self.runs[self.active];
        let cursor = run.minutes[run.scan];
        // A zoom belongs to the run it was made on.
        let reset = !(self.plotted.as_ref()).is_some_and(|p| p.ptr_eq(&Arc::downgrade(&run.data)));
        if reset {
            self.plotted = Some(Arc::downgrade(&run.data));
        }
        let plots = if self.trace == Trace::Both { 2.0 } else { 1.0 };
        let gap = ui.spacing().item_spacing.y * (plots - 1.0);
        let height = ((ui.available_height() - gap) / plots).max(64.0);
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
                ("tic", reset),
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
                ("sic", reset),
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
        self.export.chromatograms = Area::since(ui, top);
        self.view.bounds_in = None;
        clicked
    }

    /// A click selects the scan at that time, snapping to the nearest component apex within
    /// `radius` (min).
    fn select_clicked(&mut self, clicked: Option<(f64, f64)>) {
        let run = &mut self.runs[self.active];
        let Some((x, radius)) = clicked else {
            return;
        };
        if let Some(s) = run.data.scan_at(x * 60.0) {
            let off = |k: &usize| (run.minutes[run.components[*k].apex_scan as usize] - x).abs();
            run.component = (0..run.components.len())
                .min_by(|a, b| off(a).total_cmp(&off(b)))
                .filter(|k| off(k) <= radius);
            run.scan = run
                .component
                .map_or(s, |k| run.components[k].apex_scan as usize);
        }
    }

    /// The scan and its retention time, typed in or dragged; the nearest scan to a time is taken.
    fn scan_entry(&mut self, ui: &mut egui::Ui) {
        let run = &mut self.runs[self.active];
        let mut scan = None;
        ui.horizontal_wrapped(|ui| {
            let label = ui.label("Scan");
            let mut n = run.scan + 1;
            if ui
                .add(egui::DragValue::new(&mut n).range(1..=run.data.scan_count()))
                .labelled_by(label.id)
                .on_hover_text("Scan whose spectrum is shown; ←/→ step through the scans")
                .changed()
            {
                scan = Some(n - 1);
            }
            let label = ui.label(format!("RT ({})", self.view.unit.label()));
            let mut t = self.view.time(run.minutes[run.scan]);
            if ui
                .add(egui::DragValue::new(&mut t).max_decimals(3))
                .labelled_by(label.id)
                .on_hover_text("Retention time of the shown scan; the nearest scan is taken")
                .changed()
            {
                scan = run.data.scan_at(self.view.minutes(t) * 60.0);
            }
        });
        if let Some(s) = scan {
            (run.scan, run.component) = (s, None);
        }
    }

    fn spectrum_plot(&mut self, ui: &mut egui::Ui) {
        // The controls scroll in a narrow window, so the plot keeps most of the height.
        let room = ui.available_height() * 0.4;
        egui::ScrollArea::both()
            .id_salt("spectrum controls")
            .max_height(room)
            .auto_shrink([false, true])
            .show(ui, |ui| {
                self.scan_entry(ui);
                export::norm_picker(ui, &mut self.export.norm);
                let run = &self.runs[self.active];
                let (title, ..) = export::shown_spectrum(run, self.export.norm);
                let cmp = detail::comparison(run, self.export.norm);
                if let Some(e) =
                    detail::controls(ui, &mut self.spectrum_mode, &title, cmp.as_ref(), &run.name)
                {
                    self.error = Some(e);
                }
            });
        let run = &self.runs[self.active];
        let (title, mz, y) = export::shown_spectrum(run, self.export.norm);
        let cmp = detail::comparison(run, self.export.norm);
        // Without a library hit only the found spectrum can be drawn.
        let mode = cmp
            .as_ref()
            .map_or(detail::Mode::Normal, |_| self.spectrum_mode);
        let top = ui.cursor().min;
        let label = match self.export.norm {
            md_core::spectrum::Normalization::None => "Intensity",
            _ => "Normalised intensity",
        };
        let (ions, charts) = match (&cmp, mode) {
            (Some((c, lib)), detail::Mode::Split | detail::Mode::Difference) => {
                (c.mz.clone(), detail::charts(mode, c, &lib.name))
            }
            _ => {
                let bars = mz.iter().zip(&y).map(|(m, i)| Bar::new(*m, *i).width(0.3));
                let chart = BarChart::new(title.clone(), bars.collect());
                (mz.clone(), vec![chart])
            }
        };
        let mut plot = Plot::new("spectrum")
            .x_axis_label("m/z")
            .y_axis_label(format!("{label} — {}", run.name));
        if mode != detail::Mode::Normal {
            plot = plot.legend(egui_plot::Legend::default());
        }
        let shown = plot.show(ui, |p| {
            for chart in charts {
                p.bar_chart(chart);
            }
            p.pointer_coordinate()
        });
        // A click on the plot sets the SIC to the nearest ion.
        if shown.response.clicked()
            && !shown.response.double_clicked()
            && let Some(at) = shown.inner
            && let Some(ion) = ions
                .iter()
                .copied()
                .min_by(|a, b| (a - at.x).abs().total_cmp(&(b - at.x).abs()))
        {
            detail::pick_ion(&mut self.mz, &mut self.trace, ion);
        }
        // Not the plot response: that is only the inner area, the PNG needs axes and titles too.
        self.export.spectrum = Area::since(ui, top);
    }
}

#[cfg(test)]
mod tests {
    use crate::App;

    #[test]
    fn a_click_snaps_to_the_nearest_component_within_the_radius() {
        let mut app = App::default();
        // Components at 1 and 2 min.
        app.runs.push(crate::results::tests::run());
        app.select_clicked(Some((1.7, 0.2)));
        assert_eq!((app.runs[0].component, app.runs[0].scan), (None, 2));
        app.select_clicked(Some((1.85, 0.2)));
        assert_eq!((app.runs[0].component, app.runs[0].scan), (Some(1), 2));
        app.select_clicked(Some((0.2, 0.25)));
        assert_eq!((app.runs[0].component, app.runs[0].scan), (None, 0));
    }
}
