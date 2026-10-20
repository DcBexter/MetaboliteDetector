//! Panel of the selected component (old component panel: S/N, sharpness, peak discrepancy, scores, ion
//! counts, peak list) and the comparison of its spectrum with the library hit (normal / split / difference).

use eframe::egui;
use egui_plot::{Bar, BarChart};
use md_core::detail::{Comparison, compare, component_detail};
use md_core::library::LibraryCompound;
use md_core::search::LibraryHit;
use md_core::spectrum::Normalization;

use crate::{Run, Trace, export, identify};

/// Old split plot: the found spectrum alone, found above and library mirrored below, or their difference.
#[derive(Clone, Copy, Default, PartialEq)]
pub enum Mode {
    #[default]
    Normal,
    Split,
    Difference,
}

impl Mode {
    const ALL: [Self; 3] = [Self::Normal, Self::Split, Self::Difference];

    fn label(self) -> &'static str {
        match self {
            Self::Normal => "Found only",
            Self::Split => "Split",
            Self::Difference => "Difference",
        }
    }
}

/// The library hit the selected component is compared with: the assigned one, else the best.
pub fn shown_hit(run: &Run) -> Option<(&LibraryCompound, &LibraryHit)> {
    let ((lib, ..), id, _) = run.ident.as_ref()?;
    let k = run.component?;
    let hit = id.assigned.get(k)?.as_ref().or(id.hits.get(k)?.first())?;
    Some((&lib.compounds[hit.compound], hit))
}

/// Found and library spectrum of the selected component on one axis, normalised as the plot is.
pub fn comparison(run: &Run, norm: Normalization) -> Option<(Comparison, &LibraryCompound)> {
    let (lib, _) = shown_hit(run)?;
    let found = &run.components[run.component?].spectrum;
    Some((compare(found, &lib.spectrum, norm), lib))
}

/// A click on an ion sets the SIC to its mass.
pub fn pick_ion(mz: &mut f32, trace: &mut Trace, mass: f64) {
    *mz = mass as f32;
    if *trace == Trace::Tic {
        *trace = Trace::Sic;
    }
}

/// Bars of the split and difference views; `Normal` is drawn by the caller.
pub fn charts(mode: Mode, c: &Comparison, library: &str) -> Vec<BarChart> {
    let bars = |ys: Vec<f64>| -> Vec<Bar> {
        (c.mz.iter().zip(ys))
            .map(|(m, y)| Bar::new(*m, y).width(0.3))
            .collect()
    };
    match mode {
        Mode::Normal => Vec::new(),
        Mode::Split => vec![
            BarChart::new("Found", bars(c.found.clone())),
            BarChart::new(
                format!("Library: {library}"),
                bars(c.library.iter().map(|v| -v).collect()),
            ),
        ],
        Mode::Difference => vec![BarChart::new("Found − library", bars(c.difference()))],
    }
}

/// Mode picker and export of both spectra above the spectrum plot; returns an export error.
pub fn controls(
    ui: &mut egui::Ui,
    mode: &mut Mode,
    title: &str,
    cmp: Option<&(Comparison, &LibraryCompound)>,
    file_stem: &str,
) -> Option<String> {
    let mut error = None;
    ui.horizontal_wrapped(|ui| {
        ui.strong(title);
        let none = "Select a component with a library hit to compare.";
        ui.add_enabled_ui(cmp.is_some(), |ui| {
            for m in Mode::ALL {
                ui.selectable_value(mode, m, m.label());
            }
            if ui
                .button("Copy both")
                .on_hover_text("Found and library spectrum as CSV text")
                .clicked()
                && let Some((c, _)) = cmp
            {
                ui.ctx().copy_text(c.csv_text());
            }
            if ui.button("Save both as CSV…").clicked()
                && let Some((c, _)) = cmp
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("CSV", &["csv"])
                    .set_file_name(format!("{file_stem}-comparison.csv"))
                    .save_file()
            {
                error = export::write_file(&path, |w| c.write_csv(w)).err();
            }
        })
        .response
        .on_disabled_hover_text(none);
        match cmp {
            Some((_, lib)) => ui.weak(format!("vs. {}", lib.name)),
            None => ui.weak("No library hit to compare"),
        };
    });
    error
}

/// Right panel: values of the selected component, its peaks and the library hits.
pub fn panel(ui: &mut egui::Ui, run: &Run, mz: &mut f32, trace: &mut Trace) {
    let Some(k) = run.component else {
        return;
    };
    crate::layout::side_panel(ui, egui::Panel::right("hits"), |ui| {
        body(ui, run, k, mz, trace);
    });
}

fn body(ui: &mut egui::Ui, run: &Run, k: usize, mz: &mut f32, trace: &mut Trace) {
    ui.strong(format!("Component {}", k + 1));
    let c = &run.components[k];
    let d = component_detail(c, &run.peaks, &run.peak_masses);
    let hit = shown_hit(run);
    egui::Grid::new("component values").show(ui, |ui| {
        let mut row = |name: &str, tip: &str, value: String| {
            ui.label(name).on_hover_text(tip);
            ui.label(value);
            ui.end_row();
        };
        row(
            "Mean S/N",
            "Mean signal-to-noise of the ion peaks",
            stat(&d, d.mean_signal_to_noise, 1),
        );
        row(
            "Peak sharpness",
            "Mean sharpness of the ion peaks",
            stat(&d, d.mean_sharpness, 2),
        );
        let disc = if d.discrepancy < 0.0 {
            "–".into()
        } else {
            format!("{:.1}", d.discrepancy)
        };
        row(
            "Peak discrepancy",
            "Mean peak quality (old \"Peak dis.indx\"); 0 is best",
            disc,
        );
        // The rewrite reads centroid spectra only, so there is no profile count (architecture.md).
        let ions = match hit {
            Some((l, _)) => format!(
                "{} found / {} library",
                c.spectrum.len(),
                l.spectrum.masses.len()
            ),
            None => format!("{} found", c.spectrum.len()),
        };
        row(
            "Ions (centroid)",
            "Number of ions in the spectrum; spectra are centroided",
            ions,
        );
        if let Some((_, h)) = hit {
            row(
                "Overall score",
                "Score of the compared library hit",
                format!("{:.3}", h.overall),
            );
            row(
                "Spectrum score",
                "Similarity of the spectra",
                format!("{:.3}", h.spec_sim),
            );
            row(
                "RI score",
                "Similarity of the retention indices",
                h.ri_sim.map_or("–".into(), |r| format!("{r:.3}")),
            );
        }
    });
    ui.separator();
    peak_list(ui, &d.peaks, mz, trace);
    ui.separator();
    identify::hits_list(ui, run);
}

/// A mean over the ion peaks, or a dash for a component that has none (saved by the old program).
fn stat(d: &md_core::detail::Detail, value: f64, decimals: usize) -> String {
    match d.peaks.is_empty() {
        true => "–".into(),
        false => format!("{value:.decimals$}"),
    }
}

/// A right-aligned number cell, so decimals line up down the column.
fn number(ui: &mut egui::Ui, text: impl Into<egui::WidgetText>) {
    let layout = egui::Layout::right_to_left(egui::Align::Center);
    ui.with_layout(layout, |ui| ui.label(text));
}

fn peak_list(
    ui: &mut egui::Ui,
    peaks: &[md_core::detail::PeakRow],
    mz: &mut f32,
    trace: &mut Trace,
) {
    ui.strong("Peaks");
    ui.weak("Click an ion to show its SIC");
    egui::ScrollArea::vertical()
        .id_salt("component peaks")
        .max_height(220.0)
        .show(ui, |ui| {
            egui::Grid::new("component peak list")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Ion");
                    for head in ["Quality", "Sharp.", "S/N", "Area"] {
                        number(ui, egui::RichText::new(head).strong());
                    }
                    ui.end_row();
                    for p in peaks {
                        let on = *trace != Trace::Tic && *mz == p.mass as f32;
                        if ui.selectable_label(on, p.mass.to_string()).clicked() {
                            pick_ion(mz, trace, f64::from(p.mass));
                        }
                        number(ui, format!("{:.1}", p.quality));
                        number(ui, format!("{:.2}", p.sharpness));
                        number(ui, format!("{:.1}", p.signal_to_noise));
                        number(ui, format!("{:.0}", p.integral));
                        ui.end_row();
                    }
                });
        });
}

#[cfg(test)]
mod tests {
    use md_core::peaks::Peak;

    use super::*;
    use crate::results::tests::run;

    #[test]
    fn panel_and_comparison_follow_the_selected_component() {
        let mut run = run();
        run.component = Some(0);
        run.peak_masses = vec![73];
        run.peaks = vec![Peak {
            row: 0,
            first_scan: 0,
            apex_scan: 1,
            last_scan: 2,
            precise_apex: None,
            signal_to_noise: 12.0,
            quality: 1.0,
            sharpness: 3.0,
        }];
        run.components[0].spectrum = vec![(73, 100), (147, 50)];
        assert!(comparison(&run, Normalization::BasePeak).is_some());
        let (mut mz, mut trace) = (147.0, Trace::Tic);
        let ctx = egui::Context::default();
        let mut out = ctx.run_ui(Default::default(), |ui| {
            panel(ui, &run, &mut mz, &mut trace);
        });
        out.textures_delta.clear();
        pick_ion(&mut mz, &mut trace, 73.0);
        assert_eq!((mz, trace == Trace::Sic), (73.0, true));
        // Without a selection or without a hit nothing can be compared.
        run.component = Some(1);
        assert!(comparison(&run, Normalization::None).is_none());
        run.component = None;
        assert!(comparison(&run, Normalization::None).is_none());
    }
}
