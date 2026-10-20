//! Detected Peaks plot (old "Detected Peaks" dock): the ion traces of the selected component over its
//! scans, overlaid; a click on a curve picks that ion like a click in the peak table.

use eframe::egui;
use egui_plot::{Line, LineStyle, Plot};
use md_core::detail::ion_profiles;

use std::sync::Arc;

use crate::{App, Trace, color, detail};

/// A click this close to a curve, vertically, picks its ion.
const PICK_PX: f32 = 12.0;

impl App {
    pub fn detected_peaks(&mut self, ui: &mut egui::Ui) {
        let run = &self.runs[self.active];
        ui.horizontal_wrapped(|ui| {
            ui.strong("Detected peaks");
            ui.checkbox(&mut self.model_profile, "Model")
                .on_hover_text("Also draw the modelled peak profile of each ion, dashed");
            ui.weak("Click a curve to show its SIC");
        });
        let Some(k) = run.component else {
            ui.weak("Select a component to see its ion traces.");
            return;
        };
        let c = &run.components[k];
        let profiles = ion_profiles(c, &run.peak_masses, &run.data);
        let picked = |mass: u32| self.trace != Trace::Tic && self.mz == mass as f32;
        let first = f64::from(c.first_scan) + 1.0;
        // One plot memory per component, so a zoom on one never hides the next one's curves.
        let mut plot = Plot::new(("detected peaks", Arc::as_ptr(&run.data), run.rev, k));
        if self.settings.appearance.axes_labels {
            plot = plot.x_axis_label("Scan").y_axis_label("Intensity");
        }
        let shown = plot
            .y_axis_formatter(|m, _| crate::view::compact(m.value, m.step_size))
            .show(ui, |p| {
                for (i, ion) in profiles.iter().enumerate() {
                    let points = |ys: &[f64]| -> Vec<[f64; 2]> {
                        (ys.iter().enumerate())
                            .map(|(s, y)| [first + s as f64, *y])
                            .collect()
                    };
                    let name = ion.mass.to_string();
                    let line = Line::new(&name, points(&ion.measured))
                        .color(color(i))
                        .highlight(picked(ion.mass));
                    p.line(line);
                    if self.model_profile && !ion.modelled.is_empty() {
                        let model = Line::new(&name, points(&ion.modelled))
                            .color(color(i))
                            .style(LineStyle::dashed_dense());
                        p.line(model);
                    }
                }
                p.pointer_coordinate()
            });
        // The curve nearest to the click at that scan, if it passes within `PICK_PX`.
        let t = &shown.transform;
        if shown.response.clicked()
            && let Some(at) = shown.inner
            && let Some(ion) = (profiles.iter())
                .filter_map(|ion| {
                    let s = (at.x - first).round();
                    let y = ion.measured.get(usize::try_from(s as i64).ok()?)?;
                    let px = (t.position_from_point_y(*y) - t.position_from_point_y(at.y)).abs();
                    Some((ion.mass, px))
                })
                .filter(|(_, px)| *px <= PICK_PX)
                .min_by(|a, b| a.1.total_cmp(&b.1))
        {
            detail::pick_ion(&mut self.mz, &mut self.trace, f64::from(ion.0));
        }
    }
}
