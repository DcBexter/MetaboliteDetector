//! Central panel: trace selection, TIC/SIC plots, scan cursor and the spectrum below them.

use eframe::egui;
use egui_plot::{Bar, BarChart, Line, Plot, PlotPoint, Points, Text, VLine};

use crate::lang::{t, tf};
use crate::plot_export::Area;
use md_core::decimate;
use std::sync::Arc;

use crate::view::{self, Nav, View};
use crate::{App, Run, Trace, color, detail, export, identify};

fn chromatogram_help() -> [&'static str; 7] {
    [
        t(
            "Click: show the scan at that time; near a dot or triangle, that component.",
            "Klick: den Scan zu dieser Zeit zeigen; nahe einem Punkt oder Dreieck diese Komponente.",
        ),
        t(
            "← / →: previous or next scan.",
            "← / →: vorheriger oder nächster Scan.",
        ),
        t(
            "Drag, or mouse wheel: move the view.",
            "Ziehen oder Mausrad: die Ansicht verschieben.",
        ),
        t(
            "Ctrl + mouse wheel or pinch: zoom; drag on an axis: stretch it.",
            "Strg + Mausrad oder Zwei-Finger-Geste: zoomen; auf einer Achse ziehen: sie dehnen.",
        ),
        t(
            "Right mouse button, drag a box: zoom into the box.",
            "Rechte Maustaste, Rahmen ziehen: in den Rahmen zoomen.",
        ),
        t(
            "Double-click, or All: show all data.",
            "Doppelklick oder Alles: alle Daten zeigen.",
        ),
        t(
            "The dots and triangles mark the detected components; the triangles are easier to hit.",
            "Die Punkte und Dreiecke markieren die erkannten Komponenten; die Dreiecke sind leichter zu treffen.",
        ),
    ]
}

fn spectrum_help() -> [&'static str; 5] {
    [
        t(
            "Click a bar: the SIC shows that m/z.",
            "Klick auf einen Balken: das SIC zeigt dieses m/z.",
        ),
        t(
            "Drag, or mouse wheel: move the view.",
            "Ziehen oder Mausrad: die Ansicht verschieben.",
        ),
        t(
            "Ctrl + mouse wheel or pinch: zoom; drag on an axis: stretch it.",
            "Strg + Mausrad oder Zwei-Finger-Geste: zoomen; auf einer Achse ziehen: sie dehnen.",
        ),
        t(
            "Right mouse button, drag a box: zoom into the box.",
            "Rechte Maustaste, Rahmen ziehen: in den Rahmen zoomen.",
        ),
        t(
            "Double-click, or All: show all masses.",
            "Doppelklick oder Alles: alle Massen zeigen.",
        ),
    ]
}

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

/// One chromatogram plot with a line per run and `marks` as dots, unzoomed when `reset`, moved by `nav`
/// (only in one of two linked plots: the link hands the change to the other); returns
/// the clicked retention time and the snap radius `SNAP_PX`, both in minutes.
#[allow(clippy::too_many_arguments)]
fn trace_plot(
    ui: &mut egui::Ui,
    (id, reset, nav): (&str, bool, Option<Nav>),
    label: &str,
    height: f32,
    runs: &[Run],
    view: &View,
    cursor_min: f64,
    trace: fn(&Run) -> &[f64],
    (marks, apexes): (&[[f64; 2]], &[f64]),
    labels: &[identify::Label],
) -> Option<(f64, f64)> {
    // Neither a run colour nor the dots' red, so the triangles stand out from the lines.
    let marker = ui.visuals().strong_text_color();
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
        .show_crosshair(view.crosshair)
        .x_axis_label(tf!(
            "Retention time ({})",
            "Retentionszeit ({})",
            view.unit.label()
        ))
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
            if let Some(nav) = nav {
                nav.apply(p);
            }
            let b = p.plot_bounds();
            // The y range of this frame, so the triangles do not lag behind an auto-fitted axis.
            let (mut y_lo, mut y_hi) = (b.min()[1], b.max()[1]);
            if view.auto_y
                && let Some((l, h)) = y_span
            {
                let pad = ((h - l) * 0.05).max(f64::EPSILON);
                (y_lo, y_hi) = (l - pad, h + pad);
                p.set_plot_bounds_y(y_lo..=y_hi);
            }
            if !marks.is_empty() {
                let at = |m: &[f64; 2]| [view.time(m[0]), m[1]];
                p.points(
                    Points::new(
                        t("components", "Komponenten"),
                        marks.iter().map(at).collect::<Vec<_>>(),
                    )
                    .radius(3.0),
                );
            }
            if !apexes.is_empty() {
                // Triangles on the x axis: bigger click targets than the dots.
                // With Auto Y off and the plot fitting its items itself, a triangle below the data
                // would lower the range every frame, so it sits on the lowest data point.
                let y = match y_span {
                    Some((low, _)) if !view.auto_y && p.auto_bounds().y => low,
                    _ => y_lo + 8.0 * (y_hi - y_lo) / p.response().rect.height().max(1.0) as f64,
                };
                let tri = apexes
                    .iter()
                    .map(|&t| [view.time(t), y])
                    .collect::<Vec<_>>();
                p.points(
                    Points::new(t("components", "Komponenten"), tri)
                        .shape(egui_plot::MarkerShape::Up)
                        .radius(5.0)
                        .color(marker),
                );
            }
            // Hide overlapping names, best score first (3.3 does the same).
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
                p.text(
                    Text::new(t("hits", "Treffer"), at.into(), text)
                        .anchor(egui::Align2::LEFT_BOTTOM),
                );
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
            self.view.nav = view::nav_buttons(ui, &chromatogram_help());
            self.view.controls(ui);
            if ui.button(t("Run info", "Laufinfo")).clicked() {
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
                // Dropping the jobs cancels them.
                if ui.button(t("Cancel", "Abbrechen")).clicked() {
                    (self.loading, self.load_started) = (Vec::new(), None);
                }
                return;
            }
            ui.label(t(
                "Open GC-MS runs (NetCDF or old .bin) to see their chromatograms and spectra.",
                "GC-MS-Läufe öffnen (NetCDF oder altes .bin), um ihre Chromatogramme und Spektren zu sehen.",
            ));
            ui.add_space(8.0);
            self.open_chromatograms_button(ui);
            self.open_library_button(ui);
        });
    }

    fn trace_controls(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(&mut self.trace, Trace::Tic, "TIC")
                .on_hover_text(t(
                    "Total ion chromatogram: the sum of all masses per scan",
                    "Totalionenchromatogramm: die Summe aller Massen pro Scan",
                ));
            ui.selectable_value(&mut self.trace, Trace::Sic, "SIC")
                .on_hover_text(t(
                    "Selected ion chromatogram: only the chosen m/z",
                    "Ionenchromatogramm: nur das gewählte m/z",
                ));
            ui.selectable_value(&mut self.trace, Trace::Both, "TIC + SIC")
                .on_hover_text(t(
                    "Both, one above the other with a shared time axis",
                    "Beide, übereinander mit gemeinsamer Zeitachse",
                ));
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
            .on_hover_text(t(
                "Ion whose intensity the SIC trace shows",
                "Ion, dessen Intensität die SIC-Spur zeigt",
            ))
            .on_disabled_hover_text(t(
                "Used by the SIC trace; pick SIC or TIC + SIC",
                "Für die SIC-Spur; SIC oder TIC + SIC wählen",
            ));
            let label = ui.add_enabled(sic, egui::Label::new("±"));
            ui.add_enabled(
                sic,
                egui::DragValue::new(&mut self.tolerance)
                    .range(0.0..=10.0)
                    .speed(0.05)
                    .max_decimals(2),
            )
            .labelled_by(label.id)
            .on_hover_text(t(
                "Window around the m/z that the SIC trace adds up",
                "Fenster um das m/z, das die SIC-Spur aufsummiert",
            ))
            .on_disabled_hover_text(t(
                "Used by the SIC trace; pick SIC or TIC + SIC",
                "Für die SIC-Spur; SIC oder TIC + SIC wählen",
            ));
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
        let reset = !(self.plotted.as_ref()).is_some_and(|p| p.ptr_eq(&Arc::downgrade(&run.data)))
            || self.view.nav == Some(Nav::All);
        if reset {
            self.plotted = Some(Arc::downgrade(&run.data));
        }
        let plots = if self.trace == Trace::Both { 2.0 } else { 1.0 };
        let gap = ui.spacing().item_spacing.y * (plots - 1.0);
        let height = ((ui.available_height() - gap) / plots).max(64.0);
        let mut clicked = None;
        let top = ui.cursor().min;
        // The old "Axes labels" switch hides the y titles of the chromatograms.
        let axes = self.settings.appearance.axes_labels;
        let y_title = |label: String| if axes { label } else { String::new() };
        let apexes: Vec<f64> = run
            .components
            .iter()
            .map(|k| run.minutes[k.apex_scan as usize])
            .collect();
        if self.trace != Trace::Sic {
            // Components of the active run as dots on its TIC; not on the derivative, whose y differs.
            // The triangles on the time axis need no y, so they stay in every view.
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
                ("tic", reset, self.view.nav),
                &y_title(self.derivative_label("TIC".into())),
                height,
                &self.runs,
                &self.view,
                cursor,
                |r| &r.tic,
                (&marks, &apexes),
                &labels,
            );
        }
        if self.trace != Trace::Tic {
            let label = format!("SIC m/z {} ± {}", self.mz, self.tolerance);
            let label = y_title(self.derivative_label(label));
            let sic = trace_plot(
                ui,
                (
                    "sic",
                    reset,
                    self.view.nav.filter(|_| self.trace == Trace::Sic),
                ),
                &label,
                height,
                &self.runs,
                &self.view,
                cursor,
                |r| &r.sic,
                (&[], &apexes),
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
            let label = ui.label(t("Scan", "Scan"));
            let mut n = run.scan + 1;
            if ui
                .add(egui::DragValue::new(&mut n).range(1..=run.data.scan_count()))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Scan whose spectrum is shown; ←/→ step through the scans",
                    "Scan, dessen Spektrum gezeigt wird; ←/→ blättern durch die Scans",
                ))
                .changed()
            {
                scan = Some(n - 1);
            }
            let label = ui.label(format!("RT ({})", self.view.unit.label()));
            let mut rt = self.view.time(run.minutes[run.scan]);
            if ui
                .add(egui::DragValue::new(&mut rt).max_decimals(3))
                .labelled_by(label.id)
                .on_hover_text(t(
                    "Retention time of the shown scan; the nearest scan is taken",
                    "Retentionszeit des gezeigten Scans; der nächste Scan wird genommen",
                ))
                .changed()
            {
                scan = run.data.scan_at(self.view.minutes(rt) * 60.0);
            }
            self.view.spectrum_nav = view::nav_buttons(ui, &spectrum_help());
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
            md_core::spectrum::Normalization::None => t("Intensity", "Intensität"),
            _ => t("Normalised intensity", "Normierte Intensität"),
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
        let look = self.settings.appearance;
        let mut plot = Plot::new("spectrum");
        if look.axes_labels {
            plot = plot
                .x_axis_label("m/z")
                .y_axis_label(format!("{label} — {}", run.name));
        }
        // The old caption marks the found spectrum only, also when the library one is drawn below it.
        let labelled: Vec<(f64, f64)> = match (&cmp, mode) {
            (_, detail::Mode::Difference) => Vec::new(),
            (Some((c, _)), detail::Mode::Split) => {
                c.mz.iter().copied().zip(c.found.clone()).collect()
            }
            _ => mz.iter().copied().zip(y.iter().copied()).collect(),
        };
        if mode != detail::Mode::Normal {
            plot = plot.legend(egui_plot::Legend::default());
        }
        let nav = self.view.spectrum_nav;
        if nav == Some(Nav::All) {
            plot = plot.reset();
        }
        let shown = plot.show(ui, |p| {
            if let Some(nav) = nav {
                nav.apply(p);
            }
            for chart in charts {
                p.bar_chart(chart);
            }
            if look.peak_labels {
                let b = p.plot_bounds();
                for (m, v) in highest_in_view(&labelled, b.min()[0], b.max()[0]) {
                    // No name: a named item would get a legend row. The f32 masses print without widening noise.
                    let text = Text::new("", PlotPoint::new(m, v), format!("{}", m as f32))
                        .anchor(egui::Align2::CENTER_BOTTOM);
                    p.text(text);
                }
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

/// Old `SpectrumPlot::plotPeakLabels`: the `LABELLED_PEAKS` highest peaks with m/z in `lo..=hi`,
/// as (m/z, intensity).
fn highest_in_view(peaks: &[(f64, f64)], lo: f64, hi: f64) -> Vec<(f64, f64)> {
    let mut shown: Vec<(f64, f64)> = (peaks.iter().copied())
        .filter(|(m, v)| (lo..=hi).contains(m) && *v > 0.0)
        .collect();
    shown.sort_by(|a, b| b.1.total_cmp(&a.1));
    shown.truncate(crate::settings::LABELLED_PEAKS);
    shown
}

#[cfg(test)]
mod tests {
    use crate::App;

    #[test]
    fn labels_go_to_the_highest_peaks_in_view() {
        let peaks = [(73.0, 50.0), (147.0, 90.0), (200.0, 10.0), (300.0, 99.0)];
        assert_eq!(
            super::highest_in_view(&peaks, 0.0, 250.0),
            [(147.0, 90.0), (73.0, 50.0), (200.0, 10.0)]
        );
        let many: Vec<(f64, f64)> = (1..=8).map(|i| (f64::from(i), f64::from(i))).collect();
        let top = super::highest_in_view(&many, 0.0, 100.0);
        assert_eq!(
            top.iter().map(|p| p.0).collect::<Vec<_>>(),
            [8.0, 7.0, 6.0, 5.0, 4.0]
        );
    }

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
