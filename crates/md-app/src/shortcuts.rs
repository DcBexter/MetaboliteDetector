//! App-wide keyboard shortcuts (S3). Menu items and keys run the same `Action`, so they cannot differ.

use eframe::egui::{self, Key, KeyboardShortcut, Modifiers};

use crate::App;
use crate::results::stepped;

pub const REPO: &str = "https://github.com/DcBexter/MetaboliteDetector";
pub const MANUAL: &str = "https://github.com/DcBexter/MetaboliteDetector#readme";

#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Action {
    Open,
    SaveComponents,
    ExportCsv,
    AddToLibrary,
    Quit,
    Manual,
}

impl Action {
    pub const ALL: [Action; 6] = [
        Action::Open,
        Action::SaveComponents,
        Action::ExportCsv,
        Action::AddToLibrary,
        Action::Quit,
        Action::Manual,
    ];

    pub fn shortcut(self) -> KeyboardShortcut {
        let (mods, key) = match self {
            Action::Open => (Modifiers::COMMAND, Key::O),
            Action::SaveComponents => (Modifiers::COMMAND, Key::S),
            Action::ExportCsv => (Modifiers::COMMAND, Key::E),
            Action::AddToLibrary => (Modifiers::COMMAND, Key::L),
            Action::Quit => (Modifiers::COMMAND, Key::Q),
            Action::Manual => (Modifiers::NONE, Key::F1),
        };
        KeyboardShortcut::new(mods, key)
    }
}

impl App {
    /// Runs the action of a pressed shortcut, then steps components with ↑/↓ while the list is hidden.
    pub fn shortcuts(&mut self, ctx: &egui::Context) {
        for action in Action::ALL {
            if ctx.input_mut(|i| i.consume_shortcut(&action.shortcut())) {
                self.run_action(ctx, action);
            }
        }
        self.step_component(ctx);
    }

    pub fn run_action(&mut self, ctx: &egui::Context, action: Action) {
        // The status bar shows an old error before the result of these two.
        if matches!(action, Action::SaveComponents | Action::ExportCsv) {
            self.error = None;
        }
        match action {
            Action::Open => self.pick_chromatograms(ctx),
            Action::SaveComponents => {
                if let Some(run) = self.runs.get(self.active) {
                    self.results.save(run, self.ri.for_run(run));
                }
            }
            Action::ExportCsv => {
                if let Some(run) = self.runs.get(self.active) {
                    let cutoff = self.settings.applied().ls_cutoff;
                    self.results.export_csv(run, self.ri.for_run(run), cutoff);
                }
            }
            Action::AddToLibrary => self.add_component_to_library(),
            Action::Quit => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
            Action::Manual => ctx.open_url(egui::OpenUrl::new_tab(MANUAL)),
        }
    }

    /// ↑/↓ with the component list hidden: previous or next component by retention time. With the
    /// list shown it steps in the order of the list (filter, sort) itself.
    fn step_component(&mut self, ctx: &egui::Context) {
        if self.panels.components || ctx.text_edit_focused() {
            return;
        }
        let Some(run) = self.runs.get_mut(self.active) else {
            return;
        };
        let step = ctx.input(|i| {
            i32::from(i.key_pressed(Key::ArrowDown)) - i32::from(i.key_pressed(Key::ArrowUp))
        });
        let to = (step != 0).then(|| stepped(run.components.len(), run.component, step));
        if let Some(k) = to.flatten() {
            run.component = Some(k);
            run.scan = run.components[k].apex_scan as usize;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn press(app: &mut App, key: Key) {
        let event = egui::Event::Key {
            key,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: Modifiers::NONE,
        };
        let ctx = egui::Context::default();
        let input = egui::RawInput {
            events: vec![event],
            ..Default::default()
        };
        ctx.run_ui(input, |ui| app.shortcuts(ui.ctx()))
            .textures_delta
            .clear();
    }

    #[test]
    fn arrows_step_components_only_while_the_list_is_hidden() {
        let mut app = App::default();
        app.runs.push(crate::results::tests::run());
        // The list steps itself when shown.
        press(&mut app, Key::ArrowDown);
        assert_eq!(app.runs[0].component, None);
        app.panels.components = false;
        press(&mut app, Key::ArrowDown);
        assert_eq!((app.runs[0].component, app.runs[0].scan), (Some(0), 1));
        press(&mut app, Key::ArrowDown);
        press(&mut app, Key::ArrowDown);
        assert_eq!(app.runs[0].component, Some(1));
        press(&mut app, Key::ArrowUp);
        assert_eq!(app.runs[0].component, Some(0));
    }

    #[test]
    fn saving_and_exporting_replace_an_old_error() {
        let mut app = App::default();
        let ctx = egui::Context::default();
        app.error = Some("Could not load x".into());
        app.run_action(&ctx, Action::SaveComponents);
        assert_eq!(app.error, None);
    }
}
