//! Shared layout, so no window or panel puts one control on top of another; the spacing itself is
//! in `style::apply`.

use eframe::egui;

use crate::lang::t;

/// A resizable list on the left and its details beside it, each with the panel margin, so the
/// details never start on the list's edge.
pub fn list_and_details<T>(
    ui: &mut egui::Ui,
    id: &str,
    state: &mut T,
    list: impl FnOnce(&mut T, &mut egui::Ui),
    details: impl FnOnce(&mut T, &mut egui::Ui),
) {
    egui::Panel::left(id)
        .resizable(true)
        .default_size(260.0)
        .min_size(160.0)
        .show(ui, |ui| list(state, ui));
    egui::CentralPanel::default().show(ui, |ui| details(state, ui));
}

/// The areas of the main window; the View menu shows or hides each, egui keeps their sizes.
/// `persist` saves this; eframe saves the sizes with egui's memory (S4, #253).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Panels {
    pub runs: bool,
    pub components: bool,
    pub chromatogram: bool,
    pub spectrum: bool,
    pub detected_peaks: bool,
    pub info: bool,
}

impl Default for Panels {
    fn default() -> Self {
        Self {
            runs: true,
            components: true,
            chromatogram: true,
            spectrum: true,
            detected_peaks: true,
            info: true,
        }
    }
}

impl Panels {
    /// Each area under its saved key; unlike the menu names these never change.
    pub fn saved(&mut self) -> [(&'static str, &mut bool); 6] {
        [
            ("runs", &mut self.runs),
            ("components", &mut self.components),
            ("chromatogram", &mut self.chromatogram),
            ("spectrum", &mut self.spectrum),
            ("detected_peaks", &mut self.detected_peaks),
            ("info", &mut self.info),
        ]
    }

    /// Each area with its View menu name and hover text, in menu order.
    pub fn toggles(&mut self) -> [(&mut bool, &'static str, &'static str); 6] {
        [
            (
                &mut self.runs,
                t("Runs", "Läufe"),
                t("The open runs, left", "Die offenen Läufe, links"),
            ),
            (
                &mut self.components,
                t("Component list", "Komponentenliste"),
                t(
                    "Every component of the active run with its library hit, left",
                    "Jede Komponente des aktiven Laufs mit ihrem Bibliothekstreffer, links",
                ),
            ),
            (
                &mut self.chromatogram,
                t("Chromatogram", "Chromatogramm"),
                t(
                    "TIC and SIC of the open runs, centre",
                    "TIC und SIC der offenen Läufe, Mitte",
                ),
            ),
            (
                &mut self.spectrum,
                t("Spectrum", "Spektrum"),
                t(
                    "Spectrum of the scan or component, below the chromatogram",
                    "Spektrum des Scans oder der Komponente, unter dem Chromatogramm",
                ),
            ),
            (
                &mut self.detected_peaks,
                t("Detected peaks", "Erkannte Peaks"),
                t(
                    "Ion traces of the selected component, beside the spectrum",
                    "Ionenspuren der gewählten Komponente, neben dem Spektrum",
                ),
            ),
            (
                &mut self.info,
                t("Component info", "Komponenteninfo"),
                t(
                    "Values, peaks and library hits of the selected component, right",
                    "Werte, Peaks und Bibliothekstreffer der gewählten Komponente, rechts",
                ),
            ),
        ]
    }
}

/// A resizable area at the side of the main window. It starts at `default` px but at most 30 % of
/// the width, where tool windows open (`ToolWindow`), and can be dragged to 40 %, so both sides
/// together leave the plots a fifth.
pub fn dock(panel: egui::Panel, ui: &egui::Ui, default: f32) -> egui::Panel {
    // Of the window, not of the room left: the second side panel starts as wide as the first.
    let w = ui.ctx().content_rect().width();
    panel
        .resizable(true)
        .default_size(default.min(w * 0.3))
        .min_size(120.0_f32.min(w * 0.3))
        .max_size(w * 0.4)
}

/// Where a tool window opens the first time; egui keeps where the user moves it later.
pub trait ToolWindow {
    /// Opens window `n` over the plots, right of the run list (a side panel, so at most 30 % of the
    /// width), each a step down and right so none opens on another. The starting `size` is capped to
    /// the room right of and below that spot (above the status bar), so egui does not push the window
    /// back over the list or the menu bar to keep it on screen.
    /// shortcut: a window whose content needs more than that room still grows over the list, as most
    /// do at 520 px; narrower content would fix it if small windows matter.
    /// A window without its own size passes egui's default `[340.0, 420.0]`: unset, egui guesses 600
    /// px wide.
    fn tool_window(self, ctx: &egui::Context, n: u8, size: [f32; 2]) -> Self;
}

impl ToolWindow for egui::Window<'_> {
    fn tool_window(self, ctx: &egui::Context, n: u8, [w, h]: [f32; 2]) -> Self {
        let (screen, step) = (ctx.content_rect(), 20.0 * f32::from(n));
        let x = screen.min.x + screen.width() * 0.3 + 16.0 + step;
        let y = screen.min.y + 40.0 + step;
        self.default_pos([x, y]).default_size([
            w.min(screen.max.x - x - 40.0),
            h.min(screen.max.y - y - 28.0),
        ])
    }
}

#[cfg(test)]
mod tests {
    use crate::*;
    use egui::accesskit::NodeId;
    use std::collections::{HashMap, HashSet};

    /// Containers hold other widgets on purpose; only the controls inside them must not overlap.
    const CONTAINERS: [&str; 5] = ["GenericContainer", "Pane", "Window", "Unknown", "Menu"];

    /// The app with every window, dialog and panel open, each with something in it.
    fn app_with_every_window() -> App {
        let mut app = App::default();
        (
            app.ri.open,
            app.settings.open,
            app.batch.open,
            app.isotopes.open,
        ) = (true, true, true, true);
        (app.editor.open, app.view.info_open, app.updater.open) = (true, true, true);
        app.tof = crate::tof_calibration::TofWindow::demo();
        (app.statistics.open, app.mid.open) = (true, true);
        app.export_all.open = true;
        let mut c = md_core::library::LibraryCompound::new("Phenylalanine 2TMS");
        (c.spectrum.masses, c.spectrum.values) = (vec![73.0, 218.0], vec![1000, 800]);
        app.editor.compounds = vec![c];
        app.editor.select(Some(0));
        // A run with an identified component selected, so the hit list shows.
        let mut run = results::tests::run();
        (run.name, run.component) = ("a_long_run_name_from_the_lab_01".into(), Some(0));
        run.analysis = settings::Analysis::new(app.settings.key());
        let settings = md_core::search::SearchSettings::default();
        app.runs.push(run);
        // Enough runs that the run list scrolls.
        for i in 0..80 {
            let mut more = results::tests::run();
            more.name = format!("run_{i:02}");
            // Analysed with the settings in use, so no background re-analysis runs (#200).
            more.analysis = settings::Analysis::new(app.settings.key());
            app.runs.push(more);
        }
        // Enough substances that the statistics table scrolls too.
        let compounds: Vec<_> = (0..40)
            .map(|i| md_core::library::LibraryCompound::new(&format!("Substance {i}")))
            .collect();
        let lib = Arc::new(identify::Library {
            path: "lib.mdlib".into(),
            name: "lib.mdlib".into(),
            prepared: md_core::search::PreparedLibrary::new(&compounds, &settings),
            compounds,
            settings,
        });
        let found = |rt_s| md_core::batch::RunResult {
            path: format!("a_long_run_name_from_the_lab_{rt_s}.cdf").into(),
            found: (0..40)
                .map(|compound| md_core::batch::Found {
                    compound,
                    score: 0.9,
                    rt_s,
                    ions: [(73, 900.0 + compound as f64), (147, 400.0)]
                        .map(|(mass, integral)| md_core::batch::IonPeak {
                            mass,
                            integral,
                            quality: 0.1,
                            signal_to_noise: 10.0,
                        })
                        .into(),
                    rescanned: false,
                    rescued: Vec::new(),
                    alternatives: Vec::new(),
                })
                .collect(),
            ..Default::default()
        };
        // Two replicate groups, by the last digit of the file name, for the statistics table.
        app.batch.settings.replicates = r"_\d(\d)\.cdf$".into();
        app.statistics.compare("0", "1");
        app.batch
            .finished(lib, [60.0, 61.0, 70.0, 71.0].map(found).into());
        app
    }

    /// What AccessKit says about the drawn widgets.
    #[derive(Default)]
    struct Nodes {
        /// Role and label or value of each node.
        names: HashMap<NodeId, (String, String)>,
        /// (combo box, its label): egui draws a labelled combo box inside its label widget.
        combo_labels: HashSet<(NodeId, NodeId)>,
    }

    /// Draws `app` at `size` until windows and side panel widths settle (a side panel takes a few
    /// passes to grow to its content), and collects the AccessKit nodes.
    fn draw(app: &mut App, size: [f32; 2]) -> (egui::Context, Nodes) {
        let ctx = egui::Context::default();
        style::apply(&ctx);
        ctx.enable_accesskit();
        let mut frame = eframe::Frame::_new_kittest();
        app.open_editor_dialogs(&ctx);
        // A stale run starts a background analysis; one finishing between frames changes the
        // status line and other text late, so the frames would depend on timing (#200).
        let key = Some(app.settings.key());
        assert!(
            app.runs.iter().all(|r| r.analysis.key == key),
            "a run would be re-analysed in the background"
        );
        let mut nodes = Nodes::default();
        for _ in 0..8 {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
                ..Default::default()
            };
            let mut out = ctx.run_ui(input, |ui| eframe::App::ui(app, ui, &mut frame));
            out.textures_delta.clear();
            for (id, node) in
                (out.platform_output.accesskit_update.into_iter()).flat_map(|u| u.nodes)
            {
                let text = node.label().or(node.value()).unwrap_or_default().to_owned();
                nodes.names.insert(id, (format!("{:?}", node.role()), text));
                if node.role() == egui::accesskit::Role::ComboBox {
                    (nodes.combo_labels).extend(node.labelled_by().iter().map(|l| (id, *l)));
                }
            }
        }
        (ctx, nodes)
    }

    /// Controls of one layer cut off at the side; clipped by a scroll area vertically is fine (also
    /// out of view entirely), and so is cut off at the side above a horizontal scroll bar.
    fn cut_off(widgets: &[egui::WidgetRect], nodes: &Nodes) -> Vec<String> {
        let name = |w: &egui::WidgetRect| nodes.names.get(&w.id.accesskit_id());
        let bars: Vec<egui::Rect> = (widgets.iter())
            .filter(|w| name(w).is_some_and(|n| n.0 == "ScrollBar"))
            .map(|w| w.interact_rect)
            .filter(|r| r.width() > r.height())
            .collect();
        let mut found = Vec::new();
        for w in widgets {
            let r = (w.rect, w.interact_rect);
            let (left, right) = (r.1.min.x > r.0.min.x + 0.5, r.1.max.x < r.0.max.x - 0.5);
            // The clip edge is the edge of the scroll area, so it must meet the bar's end; on the
            // right a vertical bar may sit between them.
            let bar = egui::style::ScrollStyle::solid().allocated_width() + 1.0;
            let scrolls = |b: &egui::Rect| {
                b.min.y >= r.0.min.y
                    && (!left || (b.min.x - r.1.min.x).abs() <= 1.0)
                    && (!right || (b.max.x - r.1.max.x).abs() <= bar)
            };
            if let Some(n) = name(w)
                && !CONTAINERS.contains(&n.0.as_str())
                && w.rect.width() > 0.0
                // Scrolled out of view vertically: hidden, so not cut off.
                // shortcut: assumes a scroll area did the hiding, so a row unique to the hidden
                // part, or content hidden by a clip that does not scroll, goes unchecked; check
                // the clip owner if such content appears.
                && r.1.height() > 0.0
                && (left || right)
                && !bars.iter().any(scrolls)
            {
                found.push(format!("{n:?} {:?} is cut off at {:?}", r.0, r.1));
            }
        }
        found
    }

    /// Controls that reach past the screen edge where nothing scrolls them back: there the screen
    /// clips them, not a scroll area inside it (a window wider or taller than the screen).
    /// shortcut: a frameless scroll area right on the screen edge would count as the screen, and a
    /// window with its own `constrain_to` rect escapes the check; neither exists in the app yet.
    fn off_screen(widgets: &[egui::WidgetRect], nodes: &Nodes, screen: egui::Rect) -> Vec<String> {
        let mut found = Vec::new();
        for w in widgets {
            let (r, i) = (w.rect, w.interact_rect);
            let at_edge = |past: bool, clip: f32, edge: f32| past && (clip - edge).abs() <= 0.5;
            let off = at_edge(r.min.x < screen.min.x - 0.5, i.min.x, screen.min.x)
                || at_edge(r.max.x > screen.max.x + 0.5, i.max.x, screen.max.x)
                || at_edge(r.min.y < screen.min.y - 0.5, i.min.y, screen.min.y)
                || at_edge(r.max.y > screen.max.y + 0.5, i.max.y, screen.max.y);
            if let Some(n) = nodes.names.get(&w.id.accesskit_id())
                && !CONTAINERS.contains(&n.0.as_str())
                && off
            {
                found.push(format!("{n:?} {r:?} reaches past the screen"));
            }
        }
        found
    }

    /// Pairs of controls of one layer that cover each other, except a combo box and its own label.
    fn covered(widgets: &[egui::WidgetRect], nodes: &Nodes) -> Vec<String> {
        let controls: Vec<_> = (widgets.iter())
            .filter_map(|w| {
                let id = w.id.accesskit_id();
                Some((id, nodes.names.get(&id)?, w.interact_rect, w.rect))
            })
            .filter(|(_, (role, _), r, _)| !CONTAINERS.contains(&role.as_str()) && r.is_positive())
            .collect();
        let mut found = Vec::new();
        for (i, (ia, a, ra, fa)) in controls.iter().enumerate() {
            for (ib, b, rb, fb) in &controls[i + 1..] {
                let own_label = nodes.combo_labels.contains(&(*ia, *ib))
                    || nodes.combo_labels.contains(&(*ib, *ia));
                // A bar floats over the edge of the scroll area, so it lies over a row that the
                // area clips at that edge; that row is half hidden by design, not covered by mistake.
                let slack = egui::style::ScrollStyle::solid().allocated_width() + 1.0;
                let clipped = |bar: &(String, String),
                               r: &egui::Rect,
                               other: &egui::Rect,
                               full: &egui::Rect| {
                    bar.0 == "ScrollBar"
                        && ((r.width() > r.height()
                            && other.max.y < full.max.y - 0.5
                            && (other.max.y - r.min.y).abs() <= slack)
                            || (r.width() < r.height()
                                && other.max.x < full.max.x - 0.5
                                && (other.max.x - r.min.x).abs() <= slack))
                };
                let under_bar = clipped(a, ra, rb, fb) || clipped(b, rb, ra, fa);
                if !own_label && !under_bar && ra.shrink(0.5).intersects(rb.shrink(0.5)) {
                    found.push(format!("{a:?} {ra:?} overlaps {b:?} {rb:?}"));
                }
            }
        }
        found
    }

    /// Draws the app with every window open and returns each control that another one in the same
    /// layer covers, that is cut off at the side or that reaches past the screen.
    fn overlaps(size: [f32; 2]) -> Vec<String> {
        let mut app = app_with_every_window();
        let (ctx, nodes) = draw(&mut app, size);
        let shown = [
            "RI calibration",
            "Settings",
            "Batch quantification",
            "Natural isotopes",
            "Run info",
            "Library editor",
            "Add component to library",
            "Copy retention times",
            "Component list",
            "Detected peaks",
            "Biomarker statistics",
            "Substance 0",
            "Mass isotopomer distributions",
            "Library hits",
            "a_long_run_name_from_the_lab_60.cdf",
        ];
        if crate::lang::de() {
            // The titles are German; at least the windows of `windows_open_beside_the_run_list_and_apart`.
            let windows = (nodes.names.values()).filter(|(role, _)| role == "Window");
            assert!(windows.count() >= 12, "German: windows not drawn");
        }
        for title in shown.into_iter().filter(|_| !crate::lang::de()) {
            let drawn = nodes
                .names
                .values()
                .any(|(_, text)| text.starts_with(title));
            assert!(
                drawn,
                "{title} was not drawn; the test would check nothing there"
            );
        }
        let screen = ctx.content_rect();
        ctx.viewport(|v| {
            (v.prev_pass.widgets.layers())
                .flat_map(|(_, widgets)| {
                    [
                        cut_off(widgets, &nodes),
                        covered(widgets, &nodes),
                        off_screen(widgets, &nodes, screen),
                    ]
                })
                .flatten()
                .collect()
        })
    }

    #[test]
    fn windows_open_beside_the_run_list_and_apart() {
        for size in [[1200.0, 800.0], [1000.0, 700.0], [800.0, 600.0]] {
            let mut app = app_with_every_window();
            let (ctx, nodes) = draw(&mut app, size);
            // The right edge of the left panel: where its run names are clipped.
            let list = ctx.viewport(|v| {
                (v.prev_pass.widgets.layers())
                    .flat_map(|(_, ws)| ws)
                    .filter(|w| {
                        (nodes.names.get(&w.id.accesskit_id())).is_some_and(|n| n.1 == "run_00")
                    })
                    .map(|w| w.interact_rect.max.x)
                    .fold(0.0, f32::max)
            });
            assert!(list > 0.0, "{size:?}: no run list");
            let windows: Vec<egui::Rect> = ctx.memory(|m| {
                (m.layer_ids().filter(|l| l.order == egui::Order::Middle))
                    .filter_map(|l| m.area_rect(l.id))
                    .collect()
            });
            assert!(
                windows.len() >= 12,
                "{size:?}: only {} windows",
                windows.len()
            );
            let over_list: Vec<_> = windows.iter().filter(|w| w.min.x < list).collect();
            assert!(
                over_list.is_empty(),
                "{size:?}: cover the run list (to {list}): {over_list:?}"
            );
            for (i, w) in windows.iter().enumerate() {
                assert!(
                    !windows[..i].iter().any(|o| o.min == w.min),
                    "{size:?}: two windows open at {:?}",
                    w.min
                );
            }
        }
    }

    #[test]
    fn no_control_covers_another() {
        for size in [[1200.0, 800.0], [800.0, 600.0], [520.0, 700.0]] {
            let found = overlaps(size);
            assert!(found.is_empty(), "{size:?}:\n{}", found.join("\n"));
        }
    }

    /// German texts are longer (S8, #257); they must fit as well.
    #[test]
    fn no_control_covers_another_in_german() {
        crate::lang::set(crate::lang::Lang::De);
        for size in [[1200.0, 800.0], [800.0, 600.0], [520.0, 700.0]] {
            let found = overlaps(size);
            assert!(found.is_empty(), "{size:?}:\n{}", found.join("\n"));
        }
    }

    /// The statistics table has more rows than fit and is wider than its window: dragging the
    /// vertical scroll bar at the window's right edge must scroll the rows.
    #[test]
    fn the_statistics_table_scrolls_by_its_bar() {
        let size = [1280.0, 800.0];
        let mut app = app_with_every_window();
        let (ctx, nodes) = draw(&mut app, size);
        // Once the table exists, close the other windows: they sit on top of the bar and take
        // the mouse.
        (
            app.ri.open,
            app.settings.open,
            app.batch.open,
            app.isotopes.open,
        ) = (false, false, false, false);
        (app.editor.open, app.view.info_open, app.mid.open) = (false, false, false);
        app.tof.open = false;
        let mut frame = eframe::Frame::_new_kittest();
        let mut step = |app: &mut App, events: Vec<egui::Event>| {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
                events,
                ..Default::default()
            };
            let mut out = ctx.run_ui(input, |ui| eframe::App::ui(app, ui, &mut frame));
            out.textures_delta.clear();
        };
        for _ in 0..3 {
            step(&mut app, vec![]);
        }
        let named = |w: &egui::WidgetRect, text: &str| {
            (nodes.names.get(&w.id.accesskit_id())).is_some_and(|n| n.1 == text)
        };
        let rows_drawn = |ctx: &egui::Context| {
            ctx.viewport(|v| {
                (v.prev_pass.widgets.layers())
                    .flat_map(|(_, ws)| ws)
                    .any(|w| named(w, "Substance 0") && w.interact_rect.is_positive())
            })
        };
        // The bar is the tall thin scroll bar in the window with the title.
        let bar = ctx
            .viewport(|v| {
                let (_, widgets) = (v.prev_pass.widgets.layers())
                    .find(|(_, ws)| ws.iter().any(|w| named(w, "Biomarker statistics")))
                    .unwrap();
                let bars = widgets.iter().filter(|w| {
                    (nodes.names.get(&w.id.accesskit_id())).is_some_and(|n| n.0 == "ScrollBar")
                });
                bars.map(|w| w.interact_rect)
                    .find(|r| r.width() < r.height() && r.is_positive())
            })
            .expect("no vertical scroll bar in the statistics window");
        assert!(rows_drawn(&ctx), "the first row is not shown to begin with");
        let (from, to) = (
            bar.center_top() + egui::vec2(0.0, 10.0),
            bar.center_bottom(),
        );
        let button = |pos, pressed| egui::Event::PointerButton {
            pos,
            button: egui::PointerButton::Primary,
            pressed,
            modifiers: egui::Modifiers::NONE,
        };
        step(&mut app, vec![egui::Event::PointerMoved(from)]);
        step(&mut app, vec![button(from, true)]);
        // The drag needs several moves before egui counts it as one.
        for k in 1..=8 {
            let y = from.y + (to.y - from.y) * k as f32 / 8.0;
            step(
                &mut app,
                vec![egui::Event::PointerMoved(egui::pos2(from.x, y))],
            );
        }
        step(&mut app, vec![button(to, false)]);
        step(&mut app, vec![]);
        assert!(
            !rows_drawn(&ctx),
            "dragging the bar did not scroll the table"
        );
    }
}
