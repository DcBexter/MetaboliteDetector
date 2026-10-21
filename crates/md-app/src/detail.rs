//! Panel of the selected component (old component panel: S/N, sharpness, peak discrepancy, scores, ion
//! counts, peak list) and the comparison of its spectrum with the library hit (normal / split / difference).

use eframe::egui;
use egui_plot::{Bar, BarChart};
use md_core::detail::{Comparison, compare, component_detail};
use md_core::library::LibraryCompound;
use md_core::search::LibraryHit;
use md_core::spectrum::Normalization;

use crate::lang::{Msg, t, tf};
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
            Self::Normal => t("Found only", "Nur gefunden"),
            Self::Split => t("Split", "Geteilt"),
            Self::Difference => t("Difference", "Differenz"),
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
        // The ids are not the names: egui_plot hides a curve by its id, which must not change with the language (#341).
        Mode::Split => vec![
            BarChart::new(t("Found", "Gefunden"), bars(c.found.clone())).id("found"),
            BarChart::new(
                tf!("Library: {library}", "Bibliothek: {library}"),
                bars(c.library.iter().map(|v| -v).collect()),
            )
            .id("library"),
        ],
        Mode::Difference => vec![
            BarChart::new(
                t("Found − library", "Gefunden − Bibliothek"),
                bars(c.difference()),
            )
            .id("difference"),
        ],
    }
}

/// Mode picker and export of both spectra above the spectrum plot; returns an export error.
pub fn controls(
    ui: &mut egui::Ui,
    mode: &mut Mode,
    title: &str,
    cmp: Option<&(Comparison, &LibraryCompound)>,
    file_stem: &str,
) -> Option<Msg> {
    let mut error = None;
    // On its own line: beside a long title the buttons are squeezed and push the controls wider.
    ui.strong(title);
    ui.horizontal_wrapped(|ui| {
        let none = t(
            "Select a component with a library hit to compare.",
            "Zum Vergleichen eine Komponente mit Bibliothekstreffer wählen.",
        );
        ui.add_enabled_ui(cmp.is_some(), |ui| {
            for m in Mode::ALL {
                ui.selectable_value(mode, m, m.label());
            }
            if ui
                .button(t("Copy both", "Beide kopieren"))
                .on_hover_text(t(
                    "Found and library spectrum as CSV text",
                    "Gefundenes und Bibliotheksspektrum als CSV-Text",
                ))
                .clicked()
                && let Some((c, _)) = cmp
            {
                ui.ctx().copy_text(c.csv_text());
            }
            if ui
                .button(t("Save both as CSV…", "Als CSV speichern…"))
                .clicked()
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
            Some((_, lib)) => ui.weak(tf!("vs. {}", "vs. {}", lib.name)),
            None => ui.weak(t(
                "No library hit to compare",
                "Kein Bibliothekstreffer zum Vergleichen",
            )),
        };
    });
    error
}

/// Right panel: the selected component next to its library entry, its peaks, the library hits and
/// the NIST search; `cal` gives the RI while no library identified the run. True when the NIST
/// search was asked for.
pub fn panel(
    ui: &mut egui::Ui,
    run: &Run,
    nist: &identify::Nist,
    (mz, trace): (&mut f32, &mut Trace),
    cal: Option<&md_core::ri::RiCalibration>,
) -> bool {
    let mut search = false;
    crate::layout::dock(egui::Panel::right("info"), ui, 340.0).show(ui, |ui| {
        egui::ScrollArea::both().show(ui, |ui| {
            let Some(k) = run.component else {
                ui.strong(t("Component info", "Komponenteninfo"));
                ui.weak(t(
                    "Select a component in the list or on the chromatogram.",
                    "Eine Komponente in der Liste oder im Chromatogramm wählen.",
                ));
                return;
            };
            body(ui, run, k, (mz, trace), cal);
            if nist.lib.is_some() {
                ui.separator();
                search = identify::nist_list(ui, run, k, nist);
            }
        });
    });
    search
}

/// A value or a dash.
pub(crate) fn num(v: Option<f64>, decimals: usize) -> String {
    v.map_or("–".into(), |v| format!("{v:.decimals$}"))
}

fn body(
    ui: &mut egui::Ui,
    run: &Run,
    k: usize,
    (mz, trace): (&mut f32, &mut Trace),
    cal: Option<&md_core::ri::RiCalibration>,
) {
    ui.strong(tf!("Component {}", "Komponente {}", k + 1));
    let c = &run.components[k];
    let d = component_detail(c, &run.peaks, &run.peak_masses);
    let hit = shown_hit(run);
    if let Some((l, _)) = hit {
        ui.label(tf!("Library: {}", "Bibliothek: {}", l.name));
    }
    // The detected component beside its library entry, as in the old Compound Information.
    egui::Grid::new("component values").show(ui, |ui| {
        ui.label("");
        number(ui, egui::RichText::new(t("Detected", "Erkannt")).strong());
        if hit.is_some() {
            number(ui, egui::RichText::new(t("Library", "Bibliothek")).strong());
        }
        ui.end_row();
        let mut row = |name: &str, tip: &str, found: String, lib: Option<String>| {
            ui.label(name).on_hover_text(tip);
            number(ui, found);
            if hit.is_some() {
                number(ui, lib.unwrap_or_default());
            }
            ui.end_row();
        };
        let lib = |f: fn(&LibraryCompound, &LibraryHit) -> String| hit.map(|(l, h)| f(l, h));
        row(
            t("RT (min)", "RT (min)"),
            t("Retention time of the apex", "Retentionszeit des Apex"),
            format!("{:.2}", run.minutes[c.apex_scan as usize]),
            lib(|l, _| num(Some(l.spectrum.rt_ms / 60_000.0).filter(|t| *t > 0.0), 2)),
        );
        row(
            t("RI", "RI"),
            t(
                "Retention index; needs an RI calibration",
                "Retentionsindex; braucht eine RI-Kalibrierung",
            ),
            num(crate::results::component_ri(run, k, cal), 1),
            lib(|l, _| num(Some(l.spectrum.ri).filter(|ri| *ri >= 0.0), 1)),
        );
        // The rewrite reads centroid spectra only, so there is no profile count (architecture.md).
        row(
            t("Ions (centroid)", "Ionen (Zentroid)"),
            t(
                "Number of ions in the spectrum; spectra are centroided",
                "Anzahl der Ionen im Spektrum; Spektren sind zentroidiert",
            ),
            c.spectrum.len().to_string(),
            lib(|l, _| l.spectrum.masses.len().to_string()),
        );
        row(
            t("Mean S/N", "Mittleres S/N"),
            t(
                "Mean signal-to-noise of the ion peaks",
                "Mittleres Signal-Rausch-Verhältnis der Ionen-Peaks",
            ),
            stat(&d, d.mean_signal_to_noise, 1),
            None,
        );
        row(
            t("Peak sharpness", "Peakschärfe"),
            t(
                "Mean sharpness of the ion peaks",
                "Mittlere Schärfe der Ionen-Peaks",
            ),
            stat(&d, d.mean_sharpness, 2),
            None,
        );
        row(
            t("Peak discrepancy", "Peak-Abweichung"),
            t(
                "Mean peak quality (old \"Peak dis.indx\"); 0 is best",
                "Mittlere Peakqualität (alt „Peak dis.indx“); 0 ist am besten",
            ),
            num(Some(d.discrepancy).filter(|v| *v >= 0.0), 1),
            None,
        );
        if hit.is_some() {
            let mut score =
                |name: &str, tip: &str, f: fn(&LibraryCompound, &LibraryHit) -> String| {
                    row(name, tip, String::new(), lib(f))
                };
            score(
                t("Score", "Score"),
                t(
                    "Overall score of the compared library hit",
                    "Gesamt-Score des verglichenen Bibliothekstreffers",
                ),
                |_, h| format!("{:.3}", h.overall),
            );
            score(
                t("Spectrum score", "Spektrum-Score"),
                t("Similarity of the spectra", "Ähnlichkeit der Spektren"),
                |_, h| format!("{:.3}", h.spec_sim),
            );
            score(
                t("RI score", "RI-Score"),
                t(
                    "Similarity of the retention indices",
                    "Ähnlichkeit der Retentionsindizes",
                ),
                |_, h| num(h.ri_sim, 3),
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
    ui.strong(t("Peaks", "Peaks"));
    ui.weak(t(
        "Click an ion to show its SIC",
        "Ein Ion anklicken, um sein SIC zu zeigen",
    ));
    egui::ScrollArea::vertical()
        .id_salt("component peaks")
        .max_height(220.0)
        .show(ui, |ui| {
            egui::Grid::new("component peak list")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong(t("Ion", "Ion"));
                    for head in [
                        t("Quality", "Qualität"),
                        t("Sharp.", "Schärfe"),
                        "S/N",
                        t("Area", "Fläche"),
                    ] {
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
    fn chart_ids_do_not_depend_on_the_language() {
        use egui_plot::PlotItem;
        let c = Comparison {
            mz: vec![73.0],
            found: vec![1.0],
            library: vec![1.0],
        };
        let ids = |lang| {
            crate::lang::set(lang);
            [Mode::Split, Mode::Difference]
                .into_iter()
                .flat_map(|m| charts(m, &c, "lib"))
                .map(|chart| PlotItem::id(&chart))
                .collect::<Vec<_>>()
        };
        let en = ids(crate::lang::Lang::En);
        assert_eq!(en.len(), 3);
        assert_eq!(en, ids(crate::lang::Lang::De));
        crate::lang::set(crate::lang::Lang::En);
    }

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
            panel(ui, &run, &Default::default(), (&mut mz, &mut trace), None);
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
