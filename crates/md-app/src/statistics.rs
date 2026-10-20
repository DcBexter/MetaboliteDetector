//! Biomarker statistics window (F14, F15): two-group table, box plots and PCA of the batch table
//! the user sees. The numbers come from `md_core::stats`.

use std::sync::mpsc;

use eframe::egui;
use egui_extras::{Column, TableBuilder};
use egui_plot::{
    Bar, BarChart, BoxElem, BoxPlot, BoxSpread, Legend, Line, Plot, PlotPoints, Points, Text,
};
use md_core::batch::Table;
use md_core::normalize::POOL_GROUP;
use md_core::stats::{self, Biomarker, Pca, Samples, TwoGroups};

use crate::export::{write_file, write_png};

#[derive(Clone, Copy, PartialEq)]
enum Tab {
    Table,
    Box,
    Pca,
}

#[derive(Clone, Copy, PartialEq)]
enum PcaView {
    Scores,
    Loadings,
    Variance,
}

pub struct Statistics {
    pub open: bool,
    /// Replicate groups of the two sides, in the order the box plot shows them.
    sides: [Vec<String>; 2],
    paste: String,
    samples: Samples,
    missing_as_zero: bool,
    /// Show only substances with p at most this.
    p_max: f64,
    sort_by_p: bool,
    tab: Tab,
    /// Substance (table row) of the box plot.
    selected: Option<usize>,
    pca_view: PcaView,
    /// Zero-based components of the score and loading plots.
    pc: [usize; 2],
    cache: Option<Cache>,
    pca: Option<(u64, bool, Result<Pca, String>)>,
    note: String,
    /// Errors of screenshot callbacks, which run after the frame that asked for them.
    errors: (mpsc::Sender<String>, mpsc::Receiver<String>),
}

/// Comparison of the table generation `0` with the groups and options in `1`.
struct Cache(u64, TwoGroups, Vec<Biomarker>);

impl Default for Statistics {
    fn default() -> Self {
        Self {
            open: false,
            sides: Default::default(),
            paste: String::new(),
            samples: Samples::Independent,
            missing_as_zero: true,
            p_max: 1.0,
            sort_by_p: true,
            tab: Tab::Table,
            selected: None,
            pca_view: PcaView::Scores,
            pc: [0, 1],
            cache: None,
            pca: None,
            note: String::new(),
            errors: mpsc::channel(),
        }
    }
}

/// Replicate groups of the table in column order, without the pool.
fn groups_of(t: &Table) -> Vec<&str> {
    let mut out: Vec<&str> = Vec::new();
    for g in &t.groups {
        if g != POOL_GROUP && !out.contains(&g.as_str()) {
            out.push(g);
        }
    }
    out
}

/// Same colour per replicate group in every plot.
fn color(i: usize) -> egui::Color32 {
    egui::ecolor::Hsva::new(i as f32 * 0.618_034, 0.7, 0.8, 1.0).into()
}

impl Statistics {
    /// `table` is the batch table as shown (normalised), with a number that changes with it.
    pub fn window(&mut self, ctx: &egui::Context, table: Option<(u64, &Table)>) {
        if let Ok(e) = self.errors.1.try_recv() {
            self.note = e;
        }
        let mut open = self.open;
        egui::Window::new("Biomarker statistics")
            .open(&mut open)
            .default_size([900.0, 620.0])
            .show(ctx, |ui| match table {
                Some((generation, t)) => self.contents(ui, generation, t),
                None => {
                    ui.label("Run a batch quantification first (Analysis → Batch quantification).");
                }
            });
        self.open = open;
    }

    fn contents(&mut self, ui: &mut egui::Ui, generation: u64, t: &Table) {
        // Groups of an earlier table that the new one lacks.
        let all = groups_of(t);
        self.sides
            .iter_mut()
            .for_each(|s| s.retain(|g| all.contains(&g.as_str())));
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.tab, Tab::Table, "Two groups");
            ui.selectable_value(&mut self.tab, Tab::Box, "Box plot");
            ui.selectable_value(&mut self.tab, Tab::Pca, "PCA");
            ui.separator();
            ui.checkbox(&mut self.missing_as_zero, "Missing value counts as 0")
                .on_hover_text("The old program stores 0 for a substance that was not found.");
        });
        ui.separator();
        match self.tab {
            Tab::Pca => self.pca_tab(ui, generation, t),
            _ => {
                self.group_picker(ui, &all);
                let c = TwoGroups {
                    group1: self.sides[0].clone(),
                    group2: self.sides[1].clone(),
                    samples: self.samples,
                    missing_as_zero: self.missing_as_zero,
                };
                if self
                    .cache
                    .as_ref()
                    .is_none_or(|k| k.0 != generation || k.1 != c)
                {
                    let r = stats::compare(t, &c);
                    self.cache = Some(Cache(generation, c, r));
                }
                if self.tab == Tab::Table {
                    self.table_tab(ui, t);
                } else {
                    self.box_tab(ui, t);
                }
            }
        }
        if !self.note.is_empty() {
            ui.label(&self.note);
        }
    }

    /// Which replicate groups are in group 1 and 2, with a field to paste group names.
    fn group_picker(&mut self, ui: &mut egui::Ui, all: &[&str]) {
        egui::CollapsingHeader::new("Groups")
            .default_open(true)
            .show(ui, |ui| {
                ui.horizontal_top(|ui| {
                    egui::Grid::new("group picker")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.strong("Replicate group");
                            ui.strong("Group 1");
                            ui.strong("Group 2");
                            ui.end_row();
                            for g in all {
                                ui.label(*g);
                                for side in 0..2 {
                                    let mut on = self.sides[side].iter().any(|s| s == g);
                                    if ui.checkbox(&mut on, "").changed() {
                                        self.sides.iter_mut().for_each(|s| s.retain(|s| s != g));
                                        if on {
                                            self.sides[side].push((*g).into());
                                        }
                                    }
                                }
                                ui.end_row();
                            }
                        });
                    ui.vertical(|ui| {
                        ui.label("Paste group names (one per line), then fill a side:");
                        ui.add(egui::TextEdit::multiline(&mut self.paste).desired_rows(4));
                        ui.horizontal(|ui| {
                            for side in 0..2 {
                                if ui.button(format!("Fill group {}", side + 1)).clicked() {
                                    self.fill(side, all);
                                }
                            }
                        });
                    });
                });
                // The order of a side is the order of its boxes.
                for side in 0..2 {
                    if self.sides[side].len() > 1 {
                        ui.horizontal_wrapped(|ui| {
                            ui.label(format!("Order of group {}:", side + 1));
                            let mut swap = None;
                            for (i, g) in self.sides[side].iter().enumerate() {
                                ui.label(g);
                                if i > 0 && ui.small_button("◀").clicked() {
                                    swap = Some(i);
                                }
                            }
                            if let Some(i) = swap {
                                self.sides[side].swap(i - 1, i);
                            }
                        });
                    }
                }
            });
        ui.horizontal(|ui| {
            ui.label("Samples");
            ui.selectable_value(&mut self.samples, Samples::Independent, "independent runs")
                .on_hover_text("Welch's t-test over all runs of each side.");
            ui.selectable_value(&mut self.samples, Samples::GroupMeans, "group means")
                .on_hover_text(
                    "Welch's t-test over one mean per replicate group (old: not independent).",
                );
            ui.selectable_value(&mut self.samples, Samples::Paired, "paired")
                .on_hover_text(
                    "Paired t-test: the i-th run of group 1 with the i-th run of group 2.",
                );
        });
    }

    /// Puts the pasted names that are replicate groups into `side` (exact match, like the old
    /// program), taking them out of the other side.
    fn fill(&mut self, side: usize, all: &[&str]) {
        let names: Vec<&str> = self
            .paste
            .split(['\n', '\r', '\t', ';', ','])
            .map(str::trim)
            .filter(|n| all.contains(n))
            .collect();
        for n in names {
            self.sides.iter_mut().for_each(|s| s.retain(|s| s != n));
            self.sides[side].push(n.into());
        }
    }

    fn table_tab(&mut self, ui: &mut egui::Ui, t: &Table) {
        let Some(cache) = &self.cache else { return };
        if self.sides.iter().any(Vec::is_empty) {
            ui.label("Choose at least one replicate group on each side.");
            return;
        }
        ui.horizontal(|ui| {
            ui.label("Show only p ≤");
            ui.add(
                egui::DragValue::new(&mut self.p_max)
                    .range(0.0..=1.0)
                    .speed(0.001)
                    .max_decimals(4),
            );
            ui.checkbox(&mut self.sort_by_p, "sort by p");
            if self.samples == Samples::Paired && cache.1.group1.len() != cache.1.group2.len() {
                ui.colored_label(
                    ui.visuals().warn_fg_color,
                    "Paired runs: the sides need the same number of replicate groups.",
                );
            }
        });
        let mut shown: Vec<&Biomarker> = cache
            .2
            .iter()
            .filter(|r| self.p_max >= 1.0 || r.p.is_some_and(|p| p <= self.p_max))
            .collect();
        if self.sort_by_p {
            shown.sort_by(|a, b| a.p.unwrap_or(2.0).total_cmp(&b.p.unwrap_or(2.0)));
        }
        ui.horizontal(|ui| {
            ui.label(format!("{} of {} substances", shown.len(), cache.2.len()));
            if ui.button("Export CSV…").clicked()
                && let Some(p) = rfd::FileDialog::new()
                    .add_filter("CSV", &["csv"])
                    .set_file_name("biomarkers.csv")
                    .save_file()
            {
                let rows: Vec<Biomarker> = shown.iter().map(|r| (*r).clone()).collect();
                self.note = match write_file(&p, |w| stats::write_csv(t, &rows, w)) {
                    Ok(()) => format!("Saved {}", p.display()),
                    Err(e) => format!("Export failed: {e}"),
                };
            }
        });
        ui.label(match self.selected.and_then(|r| t.rows.get(r)) {
            Some(s) => format!("Selected: {} (see Box plot)", s.name),
            None => "Click a substance to see its box plot.".into(),
        });
        let mut clicked = None;
        TableBuilder::new(ui)
            .id_salt("biomarker table")
            .striped(true)
            .sense(egui::Sense::click())
            .column(Column::auto().at_least(160.0))
            .columns(Column::auto().at_least(60.0), 9)
            .header(20.0, |mut h| {
                for title in [
                    "Substance",
                    "RT (min)",
                    "Hits 1",
                    "Hits 2",
                    "Mean 1",
                    "Mean 2",
                    "SD 1",
                    "SD 2",
                    "Ratio",
                    "p-value",
                ] {
                    h.col(|ui| {
                        ui.strong(title);
                    });
                }
            })
            .body(|body| {
                body.rows(18.0, shown.len(), |mut row| {
                    let r = shown[row.index()];
                    let s = &t.rows[r.row];
                    row.set_selected(self.selected == Some(r.row));
                    let cells = [
                        s.name.clone(),
                        format!("{:.2}", s.rt_s / 60.0),
                        r.hits[0].to_string(),
                        r.hits[1].to_string(),
                        format!("{:.4}", r.mean[0]),
                        format!("{:.4}", r.mean[1]),
                        format!("{:.4}", r.sd[0]),
                        format!("{:.4}", r.sd[1]),
                        format!("{:.4}", r.ratio),
                        r.p.map_or("–".into(), |p| {
                            if p < 1e-4 {
                                format!("{p:.1e}")
                            } else {
                                format!("{p:.4}")
                            }
                        }),
                    ];
                    for (i, c) in cells.into_iter().enumerate() {
                        row.col(|ui| {
                            // Numbers right-aligned so the digits line up.
                            let layout = if i == 0 {
                                egui::Layout::left_to_right(egui::Align::Center)
                            } else {
                                egui::Layout::right_to_left(egui::Align::Center)
                            };
                            ui.with_layout(layout, |ui| {
                                ui.add(egui::Label::new(c).selectable(false));
                            });
                        });
                    }
                    if row.response().clicked() {
                        clicked = Some(r.row);
                    }
                });
            });
        if let Some(c) = clicked {
            self.selected = Some(c);
        }
    }

    fn box_tab(&mut self, ui: &mut egui::Ui, t: &Table) {
        let Some(row) = self.selected.filter(|r| *r < t.rows.len()) else {
            ui.label("Click a substance in the Two groups table.");
            return;
        };
        if self.sides.iter().all(Vec::is_empty) {
            ui.label("Choose replicate groups above to draw their boxes.");
            return;
        }
        ui.label(format!(
            "{} — summed quant ion integrals per run",
            t.rows[row].name
        ));
        let groups: Vec<(usize, &String)> = (0..2)
            .flat_map(|s| self.sides[s].iter().map(move |g| (s, g)))
            .collect();
        let names: Vec<String> = groups.iter().map(|(_, g)| (*g).clone()).collect();
        let boxes: Vec<BoxElem> = groups
            .iter()
            .enumerate()
            .filter_map(|(i, (side, g))| {
                let v = stats::group_values(t, row, std::slice::from_ref(*g), self.missing_as_zero);
                let f = stats::five_numbers(&v)?;
                Some(
                    BoxElem::new(i as f64, BoxSpread::new(f.min, f.q1, f.median, f.q3, f.max))
                        .name(*g)
                        .fill(color(*side * 3 + 1).gamma_multiply(0.5))
                        .stroke(egui::Stroke::new(1.5, color(*side * 3 + 1)))
                        .box_width(0.5)
                        .whisker_width(0.3),
                )
            })
            .collect();
        let plot = Plot::new("biomarker box plot")
            .x_axis_label("Replicate group")
            .y_axis_label("Intensity")
            .x_axis_formatter(move |m, _| {
                let i = m.value.round();
                if (m.value - i).abs() < 1e-6 && i >= 0.0 {
                    names.get(i as usize).cloned().unwrap_or_default()
                } else {
                    String::new()
                }
            })
            .show(ui, |p| p.box_plot(BoxPlot::new("boxes", boxes)));
        self.png_button(ui, plot.response.rect, "box-plot");
    }

    fn pca_tab(&mut self, ui: &mut egui::Ui, generation: u64, t: &Table) {
        if self
            .pca
            .as_ref()
            .is_none_or(|p| p.0 != generation || p.1 != self.missing_as_zero)
        {
            let r = stats::pca(t, POOL_GROUP, self.missing_as_zero);
            self.pca = Some((generation, self.missing_as_zero, r));
        }
        let Some((_, _, result)) = &self.pca else {
            return;
        };
        let p = match result {
            Ok(p) => p,
            Err(e) => {
                ui.colored_label(ui.visuals().error_fg_color, e);
                return;
            }
        };
        ui.label(format!(
            "{} runs (without {POOL_GROUP}), {} substances, z-scored per substance",
            p.samples.len(),
            p.rows.len()
        ));
        let k = p.proportion.len();
        self.pc.iter_mut().for_each(|c| *c = (*c).min(k - 1));
        ui.horizontal(|ui| {
            ui.selectable_value(&mut self.pca_view, PcaView::Scores, "Scores");
            ui.selectable_value(&mut self.pca_view, PcaView::Loadings, "Loadings");
            ui.selectable_value(&mut self.pca_view, PcaView::Variance, "Explained variance");
            if self.pca_view != PcaView::Variance {
                for (axis, c) in ["x", "y"].into_iter().zip(&mut self.pc) {
                    egui::ComboBox::from_label(axis)
                        .selected_text(format!("PC{}", *c + 1))
                        .show_ui(ui, |ui| {
                            for i in 0..k {
                                ui.selectable_value(
                                    c,
                                    i,
                                    format!("PC{} ({:.1} %)", i + 1, 100.0 * p.proportion[i]),
                                );
                            }
                        });
                }
            }
        });
        let rect = match self.pca_view {
            PcaView::Scores => scores_plot(ui, t, p, self.pc),
            PcaView::Loadings => loadings_plot(ui, t, p, self.pc),
            PcaView::Variance => variance_plot(ui, p),
        };
        let name = match self.pca_view {
            PcaView::Scores => "pca-scores",
            PcaView::Loadings => "pca-loadings",
            PcaView::Variance => "pca-variance",
        };
        self.png_button(ui, rect, name);
    }

    /// Saves the plot area as PNG: a crop of egui's own screenshot, like the chromatogram export.
    fn png_button(&mut self, ui: &mut egui::Ui, rect: egui::Rect, what: &str) {
        if ui.button("Save PNG…").clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("PNG image", &["png"])
                .set_file_name(format!("{what}.png"))
                .save_file()
        {
            // The screenshot is read after the next frame, so the dialog is gone from it.
            let (ctx, errors) = (ui.ctx().clone(), self.errors.0.clone());
            let ppp = ctx.pixels_per_point();
            ctx.request_screenshot(move |shot| {
                if let Err(e) = write_png(&path, &shot, rect, ppp) {
                    errors.send(e).ok();
                }
            });
        }
    }
}

/// Samples in the plane of two components, one colour per replicate group. Returns the plot area.
fn scores_plot(ui: &mut egui::Ui, t: &Table, p: &Pca, pc: [usize; 2]) -> egui::Rect {
    let label = |c: usize| format!("PC{} ({:.1} %)", c + 1, 100.0 * p.proportion[c]);
    let mut groups: Vec<&str> = Vec::new();
    for c in &p.samples {
        if !groups.contains(&t.groups[*c].as_str()) {
            groups.push(&t.groups[*c]);
        }
    }
    Plot::new("pca scores")
        .legend(Legend::default())
        .data_aspect(1.0)
        .x_axis_label(label(pc[0]))
        .y_axis_label(label(pc[1]))
        .show(ui, |plot| {
            for (gi, g) in groups.iter().enumerate() {
                let at = |i: usize| [p.scores[pc[0]][i], p.scores[pc[1]][i]];
                let mine: Vec<usize> = (0..p.samples.len())
                    .filter(|i| t.groups[p.samples[*i]] == *g)
                    .collect();
                plot.points(
                    Points::new(*g, mine.iter().map(|i| at(*i)).collect::<Vec<_>>())
                        .radius(4.0)
                        .color(color(gi)),
                );
                for i in mine {
                    let text = egui::RichText::new(&t.runs[p.samples[i]]).color(color(gi));
                    plot.text(Text::new("", at(i).into(), text).anchor(egui::Align2::LEFT_BOTTOM));
                }
            }
        })
        .response
        .rect
}

/// The substances with the largest weight in the two components ("impact of metabolites").
fn loadings_plot(ui: &mut egui::Ui, t: &Table, p: &Pca, pc: [usize; 2]) -> egui::Rect {
    const SHOWN: usize = 25;
    let weight = |j: usize| p.loadings[pc[0]][j].hypot(p.loadings[pc[1]][j]);
    let mut order: Vec<usize> = (0..p.rows.len()).collect();
    order.sort_by(|a, b| weight(*b).total_cmp(&weight(*a)));
    order.truncate(SHOWN);
    let names: Vec<String> = order
        .iter()
        .map(|j| t.rows[p.rows[*j]].name.clone())
        .collect();
    // Horizontal bars, largest on top: long substance names read on the y axis.
    let bars = |c: usize, shift: f64| -> Vec<Bar> {
        order
            .iter()
            .enumerate()
            .map(|(i, j)| Bar::new(-(i as f64) + shift, p.loadings[c][*j]).width(0.35))
            .collect()
    };
    Plot::new("pca loadings")
        .legend(Legend::default())
        .x_axis_label("Loading")
        .y_axis_formatter(move |m, _| {
            let i = (-m.value).round();
            if (-m.value - i).abs() < 1e-6 && i >= 0.0 {
                names.get(i as usize).cloned().unwrap_or_default()
            } else {
                String::new()
            }
        })
        .show(ui, |plot| {
            plot.bar_chart(
                BarChart::new(format!("PC{}", pc[0] + 1), bars(pc[0], 0.2))
                    .horizontal()
                    .color(color(0)),
            );
            plot.bar_chart(
                BarChart::new(format!("PC{}", pc[1] + 1), bars(pc[1], -0.2))
                    .horizontal()
                    .color(color(1)),
            );
        })
        .response
        .rect
}

/// Share of the variance per component and the running total.
fn variance_plot(ui: &mut egui::Ui, p: &Pca) -> egui::Rect {
    let mut sum = 0.0;
    let cumulative: PlotPoints = p
        .proportion
        .iter()
        .enumerate()
        .map(|(i, v)| {
            sum += v;
            [(i + 1) as f64, 100.0 * sum]
        })
        .collect();
    let bars: Vec<Bar> = p
        .proportion
        .iter()
        .enumerate()
        .map(|(i, v)| Bar::new((i + 1) as f64, 100.0 * v).width(0.6))
        .collect();
    Plot::new("pca variance")
        .legend(Legend::default())
        .x_axis_label("Component")
        .y_axis_label("Explained variance (%)")
        .show(ui, |plot| {
            plot.bar_chart(BarChart::new("single", bars).color(color(0)));
            plot.line(Line::new("cumulative", cumulative).color(color(1)));
        })
        .response
        .rect
}
