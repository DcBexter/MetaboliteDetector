//! Shared layout, so no window or panel puts one control on top of another; the spacing itself is
//! in `style::apply`.

use eframe::egui;

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

/// A panel beside the plots that grows to its content like a plain panel, but to at most 30 % of
/// the width, so the plots and their controls keep room in a narrow window; the rest scrolls.
pub fn side_panel(
    ui: &mut egui::Ui,
    panel: fn(&'static str) -> egui::Panel,
    id: &'static str,
    add: impl FnOnce(&mut egui::Ui),
) {
    // Width the content asked for in the last pass: a scroll area takes no more than it is given.
    let key = egui::Id::unique(("side panel wanted width", id));
    let wanted: f32 = ui.data(|d| d.get_temp(key)).unwrap_or_default();
    let max = ui.available_width() * 0.3;
    let shown = panel(id)
        .max_size(max)
        .min_size(wanted.min(max))
        .show(ui, |ui| {
            let content = egui::ScrollArea::both().show(ui, add).content_size.x;
            (ui.max_rect().width(), content)
        });
    let (inner, content) = shown.inner;
    // The panel is its content plus the frame margin.
    let now = (content + shown.response.rect.width() - inner).ceil();
    if (now - wanted).abs() > 0.5 {
        ui.data_mut(|d| d.insert_temp(key, now));
        ui.ctx().request_repaint();
    }
}

/// Where a tool window opens the first time; egui keeps where the user moves it later.
pub trait ToolWindow {
    /// Opens window `n` over the plots, right of the run list (a side panel, so at most 30 % of the
    /// width), each a step down and right so none opens on another. The starting `size` is capped to
    /// the room there, so egui does not push the window back over the list to keep it on screen.
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
        self.default_pos([x, screen.min.y + 40.0 + step])
            .default_size([w.min(screen.max.x - x - 40.0), h])
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
        (app.editor.open, app.view.info_open) = (true, true);
        (app.results.open, app.statistics.open, app.mid.open) = (true, true, true);
        let mut c = md_core::library::LibraryCompound::new("Phenylalanine 2TMS");
        (c.spectrum.masses, c.spectrum.values) = (vec![73.0, 218.0], vec![1000, 800]);
        app.editor.compounds = vec![c];
        app.editor.select(Some(0));
        // A run with an identified component selected, so the hit list shows.
        let mut run = results::tests::run();
        (run.name, run.component) = ("a_long_run_name_from_the_lab_01".into(), Some(0));
        run.analysis = settings::Analysis::new(app.settings.key());
        let lib = run.ident.as_ref().unwrap().0.0.clone();
        app.runs.push(run);
        // Enough runs that the run list scrolls.
        for i in 0..80 {
            let mut more = results::tests::run();
            more.name = format!("run_{i:02}");
            app.runs.push(more);
        }
        let found = |rt_s| md_core::batch::RunResult {
            path: format!("a_long_run_name_from_the_lab_{rt_s}.cdf").into(),
            found: vec![md_core::batch::Found {
                compound: 0,
                score: 0.9,
                rt_s,
                ions: [(73, 900.0), (147, 400.0)]
                    .map(|(mass, integral)| md_core::batch::IonPeak {
                        mass,
                        integral,
                        quality: 0.1,
                        signal_to_noise: 10.0,
                    })
                    .into(),
                rescanned: false,
                rescued: Vec::new(),
            }],
            ..Default::default()
        };
        app.batch.finished(lib, vec![found(60.0), found(61.0)]);
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
                Some((id, nodes.names.get(&id)?, w.interact_rect))
            })
            .filter(|(_, (role, _), r)| !CONTAINERS.contains(&role.as_str()) && r.is_positive())
            .collect();
        let mut found = Vec::new();
        for (i, (ia, a, ra)) in controls.iter().enumerate() {
            for (ib, b, rb) in &controls[i + 1..] {
                let own_label = nodes.combo_labels.contains(&(*ia, *ib))
                    || nodes.combo_labels.contains(&(*ib, *ia));
                if !own_label && ra.shrink(0.5).intersects(rb.shrink(0.5)) {
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
            "Results",
            "Biomarker statistics",
            "Mass isotopomer distributions",
            "Library hits",
            "a_long_run_name_from_the_lab_60.cdf",
        ];
        for title in shown {
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
            let (ctx, _) = draw(&mut app, size);
            let key = egui::Id::unique(("side panel wanted width", "runs"));
            let list = ctx
                .data(|d| d.get_temp::<f32>(key))
                .unwrap()
                .min(size[0] * 0.3);
            let windows: Vec<egui::Rect> = ctx.memory(|m| {
                (m.layer_ids().filter(|l| l.order == egui::Order::Middle))
                    .filter_map(|l| m.area_rect(l.id))
                    .collect()
            });
            assert!(
                windows.len() >= 11,
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
}
