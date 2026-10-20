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
        (app.results.open, app.statistics.open) = (true, true);
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

    /// Controls of one layer cut off at the side; clipped by a scroll area vertically is fine,
    /// and so is cut off at the side above a horizontal scroll bar.
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
                && (left || right)
                && !bars.iter().any(scrolls)
            {
                found.push(format!("{n:?} {:?} is cut off at {:?}", r.0, r.1));
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
    /// layer covers or that is cut off at the side.
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
        ctx.viewport(|v| {
            (v.prev_pass.widgets.layers())
                .flat_map(|(_, widgets)| [cut_off(widgets, &nodes), covered(widgets, &nodes)])
                .flatten()
                .collect()
        })
    }

    #[test]
    fn no_control_covers_another() {
        for size in [[1200.0, 800.0], [800.0, 600.0], [520.0, 700.0]] {
            let found = overlaps(size);
            assert!(found.is_empty(), "{size:?}:\n{}", found.join("\n"));
        }
    }
}
