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

/// A panel beside the plots: at most 30 % of the width, so the plots and their controls keep room
/// in a narrow window; content wider or longer than the panel scrolls.
pub fn side_panel(ui: &mut egui::Ui, panel: egui::Panel, add: impl FnOnce(&mut egui::Ui)) {
    panel
        .max_size(ui.available_width() * 0.3)
        .show(ui, |ui| egui::ScrollArea::both().show(ui, add));
}

#[cfg(test)]
mod tests {
    use crate::*;
    use std::collections::HashMap;

    /// Containers hold other widgets on purpose; only the controls inside them must not overlap.
    const CONTAINERS: [&str; 5] = ["GenericContainer", "Pane", "Window", "Unknown", "Menu"];

    /// Draws the app with its windows open and returns each control that another one in the same
    /// layer covers or that is cut off at the side.
    fn overlaps(size: [f32; 2]) -> Vec<String> {
        let mut app = App::default();
        (
            app.ri.open,
            app.settings.open,
            app.batch.open,
            app.isotopes.open,
        ) = (true, true, true, true);
        app.editor.open = true;
        app.view.info_open = true;
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
        app.results.open = true;
        app.statistics.open = true;
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
            }],
            ..Default::default()
        };
        app.batch.finished(lib, vec![found(60.0), found(61.0)]);
        let ctx = egui::Context::default();
        style::apply(&ctx);
        ctx.enable_accesskit();
        let mut frame = eframe::Frame::_new_kittest();
        app.open_editor_dialogs(&ctx);
        let mut names = HashMap::new();
        // (control, its label): egui draws a labelled combo box inside its label widget.
        let mut labelled = std::collections::HashSet::new();
        for _ in 0..4 {
            let input = egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, size.into())),
                ..Default::default()
            };
            let mut out = ctx.run_ui(input, |ui| eframe::App::ui(&mut app, ui, &mut frame));
            out.textures_delta.clear();
            for (id, node) in out
                .platform_output
                .accesskit_update
                .into_iter()
                .flat_map(|u| u.nodes)
            {
                let text = node.label().or(node.value()).unwrap_or_default().to_owned();
                names.insert(id, (format!("{:?}", node.role()), text));
                labelled.extend(node.labelled_by().iter().map(|l| (id, *l)));
            }
        }
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
            let drawn = names.values().any(|(_, text)| text.starts_with(title));
            assert!(
                drawn,
                "{title} was not drawn; the test would check nothing there"
            );
        }
        ctx.viewport(|v| {
            let mut found = Vec::new();
            for (_, widgets) in v.prev_pass.widgets.layers() {
                // Horizontal scroll bars: what they scroll to may be cut off at the side.
                let bars: Vec<egui::Rect> = (widgets.iter())
                    .filter(|w| {
                        let r = w.interact_rect;
                        names
                            .get(&w.id.accesskit_id())
                            .is_some_and(|n| n.0 == "ScrollBar")
                            && r.width() > r.height()
                    })
                    .map(|w| w.interact_rect)
                    .collect();
                // Clipped by a scroll area vertically is fine; cut off at the side only above a
                // horizontal scroll bar.
                for w in widgets {
                    let r = (w.rect, w.interact_rect);
                    let scrolls = |b: &egui::Rect| {
                        // The clip edge, which is all that is left of a control clipped entirely.
                        let x = b.x_range().expand(1.0);
                        b.min.y >= r.0.min.y && (x.contains(r.1.min.x) || x.contains(r.1.max.x))
                    };
                    if let Some(name) = names.get(&w.id.accesskit_id())
                        && !CONTAINERS.contains(&name.0.as_str())
                        && w.rect.width() > 0.0
                        && (r.1.min.x > r.0.min.x + 0.5 || r.1.max.x < r.0.max.x - 0.5)
                        && !bars.iter().any(scrolls)
                    {
                        found.push(format!("{name:?} {:?} is cut off at {:?}", r.0, r.1));
                    }
                }
                let controls: Vec<_> = widgets
                    .iter()
                    .filter_map(|w| {
                        let id = w.id.accesskit_id();
                        Some((id, names.get(&id)?, w.interact_rect))
                    })
                    .filter(|(_, (role, _), r)| {
                        !CONTAINERS.contains(&role.as_str()) && r.is_positive()
                    })
                    .collect();
                for (i, (ia, a, ra)) in controls.iter().enumerate() {
                    for (ib, b, rb) in &controls[i + 1..] {
                        let pair = labelled.contains(&(*ia, *ib)) || labelled.contains(&(*ib, *ia));
                        if !pair && ra.shrink(0.5).intersects(rb.shrink(0.5)) {
                            found.push(format!("{a:?} {ra:?} overlaps {b:?} {rb:?}"));
                        }
                    }
                }
            }
            found
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
