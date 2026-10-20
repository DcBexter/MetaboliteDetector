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
        let data = md_core::Chromatogram {
            scan_times_s: vec![60.0, 61.0],
            offsets: vec![0, 1, 2],
            mz: vec![73.0, 73.0],
            intensity: vec![1.0, 2.0],
        };
        app.runs.push(Run {
            path: "a_long_run_name_from_the_lab_01.cdf".into(),
            name: "a_long_run_name_from_the_lab_01".into(),
            minutes: vec![1.0, 61.0 / 60.0],
            tic: data.tic(),
            data: std::sync::Arc::new(data),
            sic: Vec::new(),
            sic_for: None,
            sic_job: None,
            components: Vec::new(),
            peak_masses: Vec::new(),
            peaks: Vec::new(),
            analysis: settings::Analysis::new(app.settings.key()),
            scan: 0,
            component: None,
            visible: true,
            ident: None,
            ident_job: None,
            saved_calibration: None,
            rev: 0,
        });
        let ctx = egui::Context::default();
        style::apply(&ctx);
        ctx.enable_accesskit();
        let mut frame = eframe::Frame::_new_kittest();
        let mut names = HashMap::new();
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
            }
        }
        ctx.viewport(|v| {
            let mut found = Vec::new();
            for (_, widgets) in v.prev_pass.widgets.layers() {
                // Clipped by a scroll area vertically is fine; cut off at the side is not.
                for w in widgets {
                    let r = (w.rect, w.interact_rect);
                    if let Some(name) = names.get(&w.id.accesskit_id())
                        && !CONTAINERS.contains(&name.0.as_str())
                        && w.rect.width() > 0.0
                        && (r.1.min.x > r.0.min.x + 0.5 || r.1.max.x < r.0.max.x - 0.5)
                    {
                        found.push(format!("{name:?} {:?} is cut off at {:?}", r.0, r.1));
                    }
                }
                let controls: Vec<_> = widgets
                    .iter()
                    .filter_map(|w| Some((names.get(&w.id.accesskit_id())?, w.interact_rect)))
                    .filter(|((role, _), r)| {
                        !CONTAINERS.contains(&role.as_str()) && r.is_positive()
                    })
                    .collect();
                for (i, (a, ra)) in controls.iter().enumerate() {
                    for (b, rb) in &controls[i + 1..] {
                        if ra.shrink(0.5).intersects(rb.shrink(0.5)) {
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
