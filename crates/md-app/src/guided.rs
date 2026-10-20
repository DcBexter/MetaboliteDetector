//! Guided mode (S7, #256): a global switch "Guided / Expert", asked at the first start, and four step
//! flows (RI calibration, component detection, batch quantification, MID). Each step draws a building
//! block of the expert window on the same state, so a flow ends in that window with its values set.
//! In Guided mode the expert options of the windows sit in a closed "Advanced" section.

use eframe::egui;

use crate::App;
use crate::layout::ToolWindow;

/// Every text this module shows, in one place for the translation (S8, #257).
pub mod text {
    use super::{Flow, Step};

    pub const ADVANCED: &str = "Advanced";
    pub const SWITCH: &str = "Guided mode";
    pub const SWITCH_HINT: &str = "Guided: the main tasks step by step, expert options under a \
        closed \"Advanced\". Expert: every option open. Every function stays reachable in both.";
    pub const MODE: &str = "Mode";
    pub const ASK_TITLE: &str = "How do you want to work?";
    pub const GUIDED: &str = "Guided";
    pub const GUIDED_HINT: &str = "Recommended to start: the main tasks step by step \
        (Analysis → Step by step), expert options under \"Advanced\".";
    pub const EXPERT: &str = "Expert";
    pub const EXPERT_HINT: &str = "Every window with all its options open, as in the old program.";
    pub const ASK_LATER: &str = "You can switch at any time in Help and in Settings.";
    pub const CONTINUE: &str = "Continue";
    pub const MENU: &str = "Step by step";
    pub const BACK: &str = "Back";
    pub const NEXT: &str = "Next";
    pub const SKIP: &str = "Skip";
    pub const SKIP_HINT: &str = "Go on without this step; you can come back with Back";
    pub const FINISH: &str = "Finish";
    pub const TO_EXPERT: &str = "Open the window";
    pub const TO_EXPERT_HINT: &str =
        "Leave the steps and go on in the full window, with the values set so far";
    pub const LIBRARY_LOADING: &str = "Library: loading…";
    pub const BACK_HINT: &str = "The step before; what you set stays";
    pub const NEXT_HINT: &str = "The next step";
    pub const CONTINUE_HINT: &str = "Start working in the chosen mode";

    pub fn active_run(name: Option<&str>) -> String {
        match name {
            Some(n) => format!("Active run: {n}"),
            None => "No run open yet".into(),
        }
    }

    pub fn library(name: Option<&str>) -> String {
        match name {
            Some(n) => format!("Library: {n}"),
            None => "No library open yet".into(),
        }
    }

    pub fn step_of(x: usize, n: usize) -> String {
        format!("Step {x} of {n}")
    }

    pub fn result((runs, errors): (usize, usize)) -> String {
        format!("Last run: {runs} runs, errors: {errors}. Finish shows the details.")
    }

    pub fn finish_hint(window: &str) -> String {
        format!("Close the steps and open the window \"{window}\" with the values set")
    }

    /// Menu entry, its hover text and the title of the steps' window.
    pub fn flow(f: Flow) -> (&'static str, &'static str, &'static str) {
        match f {
            Flow::Ri => (
                "RI calibration…",
                "Retention indices from an alkane mix or a .calib file",
                "RI calibration step by step",
            ),
            Flow::Detection => (
                "Component detection…",
                "Detect the components of many runs and save them as .cmp files",
                "Component detection step by step",
            ),
            Flow::Batch => (
                "Batch quantification…",
                "Quantify the library's substances across many runs",
                "Batch quantification step by step",
            ),
            Flow::Mid => (
                "Mass isotopomer distributions…",
                "Isotopomer distributions of the library's fragments in labelling experiments",
                "MIDs step by step",
            ),
        }
    }

    /// Title, what to do, and why Next waits.
    pub fn step(s: Step) -> (&'static str, &'static str, &'static str) {
        match s {
            Step::OpenRun => (
                "Open the alkane run",
                "Open the run of your n-alkane mix. It becomes the active run.",
                "Open a run first",
            ),
            Step::FindAlkanes => (
                "Find the alkanes",
                "Set the marker ion and the first and last alkane of your mix, then find them \
                 in the active run. A saved .calib file works too.",
                "Find the alkanes or load a .calib file first",
            ),
            Step::CheckPoints => (
                "Check the points",
                "Each alkane gives one point from retention time to retention index. Switch \
                 off a point that was found at the wrong peak.",
                "",
            ),
            Step::Library => (
                "Open the library",
                "The substances are searched in this library. It is opened again at the next start.",
                "Open a library first",
            ),
            Step::BatchRuns => (
                "Add the runs",
                "Add every run of the experiment.",
                "Add runs first",
            ),
            Step::SaveCmp => (
                "Where the components go",
                "The components of each run are saved as run.cmp, next to the run or in a \
                 folder you choose. They load faster later.",
                "",
            ),
            Step::QuantSettings => (
                "Groups and filters",
                "Name the replicate groups by a pattern on the file name or a sample table, and \
                 how reproducible a substance must be.",
                "",
            ),
            Step::StartBatch => (
                "Start",
                "Start the scan. It goes on in the background; the window shows the result.",
                "",
            ),
            Step::MidMethod => (
                "Correction",
                "Choose how the labelled spectra are corrected and the tracer element.",
                "",
            ),
            Step::MidRuns => (
                "Add the runs",
                "Add the labelled runs and, for the reference correction, unlabelled runs of the \
                 same substances.",
                "Add the runs first, or fix what the step shows",
            ),
            Step::StartMid => (
                "Start",
                "Start the calculation. It goes on in the background; the window shows the result.",
                "",
            ),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum Mode {
    Guided,
    /// The default where nothing asked yet (tests), so every option is open as before.
    #[default]
    Expert,
}

impl Mode {
    /// The value saved under `mode` (S4).
    pub fn key(self) -> &'static str {
        match self {
            Mode::Guided => "guided",
            Mode::Expert => "expert",
        }
    }

    pub fn from_key(v: &str) -> Option<Mode> {
        [Mode::Guided, Mode::Expert]
            .into_iter()
            .find(|m| m.key() == v)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Flow {
    Ri,
    Detection,
    Batch,
    Mid,
}

impl Flow {
    pub const ALL: [Flow; 4] = [Flow::Ri, Flow::Detection, Flow::Batch, Flow::Mid];

    fn steps(self) -> &'static [Step] {
        use Step::*;
        match self {
            Flow::Ri => &[OpenRun, FindAlkanes, CheckPoints],
            Flow::Detection => &[Library, BatchRuns, SaveCmp, StartBatch],
            Flow::Batch => &[Library, BatchRuns, QuantSettings, StartBatch],
            Flow::Mid => &[Library, MidMethod, MidRuns, StartMid],
        }
    }

    /// The expert window the flow ends in: its title and its open flag.
    fn window(self, app: &mut App) -> (&'static str, &mut bool) {
        match self {
            Flow::Ri => ("RI calibration", &mut app.ri.open),
            Flow::Detection | Flow::Batch => ("Batch quantification", &mut app.batch.open),
            Flow::Mid => ("Mass isotopomer distributions", &mut app.mid.open),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Step {
    OpenRun,
    FindAlkanes,
    CheckPoints,
    Library,
    BatchRuns,
    SaveCmp,
    QuantSettings,
    StartBatch,
    MidMethod,
    MidRuns,
    StartMid,
}

#[derive(Default)]
pub struct Guided {
    pub mode: Mode,
    /// The first start asks for the mode (nothing saved yet).
    pub ask: bool,
    /// The running flow and its step (from 0).
    pub flow: Option<(Flow, usize)>,
    /// The detection flow switched ".cmp" saving on; Finish, cancel or another flow switch it off again.
    cmp_switched_on: bool,
}

fn mode_id() -> egui::Id {
    egui::Id::unique("guided mode")
}

/// Makes the mode known to the windows' "Advanced" sections for this frame.
pub fn publish(ctx: &egui::Context, mode: Mode) {
    ctx.data_mut(|d| d.insert_temp(mode_id(), mode));
}

/// Expert options: closed in Guided mode, open in Expert mode, and reachable in both.
pub fn advanced(ui: &mut egui::Ui, id: &str, add: impl FnOnce(&mut egui::Ui)) {
    let mode = ui
        .ctx()
        .data(|d| d.get_temp::<Mode>(mode_id()))
        .unwrap_or_default();
    // Salted with the mode, so a switch opens or closes every section anew.
    egui::CollapsingHeader::new(text::ADVANCED)
        .id_salt((id, mode))
        .default_open(mode == Mode::Expert)
        .show(ui, add);
}

/// The switch in Help and in Settings.
pub fn switch(ui: &mut egui::Ui, mode: &mut Mode) {
    let mut on = *mode == Mode::Guided;
    if ui
        .checkbox(&mut on, text::SWITCH)
        .on_hover_text(text::SWITCH_HINT)
        .changed()
    {
        *mode = if on { Mode::Guided } else { Mode::Expert };
    }
}

impl App {
    /// Analysis → Step by step: the four flows, in both modes.
    pub fn flows_menu(&mut self, ui: &mut egui::Ui) {
        ui.menu_button(text::MENU, |ui| {
            for f in Flow::ALL {
                let (name, hint, _) = text::flow(f);
                if ui.button(name).on_hover_text(hint).clicked() {
                    self.start_flow(f);
                }
            }
        });
    }

    pub fn start_flow(&mut self, f: Flow) {
        self.cancel_flow();
        self.guided.cmp_switched_on = f == Flow::Detection && self.batch.save_cmp_by_default();
        self.guided.flow = Some((f, 0));
    }

    /// The first-start question and the window of the running flow.
    pub fn guided_windows(&mut self, ctx: &egui::Context) {
        if self.guided.ask {
            let g = &mut self.guided;
            let modal = egui::Modal::new(egui::Id::unique("ask mode")).show(ctx, |ui| {
                ui.set_max_width(420.0);
                ui.heading(text::ASK_TITLE);
                ui.add_space(8.0);
                for (mode, name, hint) in [
                    (Mode::Guided, text::GUIDED, text::GUIDED_HINT),
                    (Mode::Expert, text::EXPERT, text::EXPERT_HINT),
                ] {
                    ui.radio_value(&mut g.mode, mode, name);
                    ui.indent(name, |ui| ui.weak(hint));
                }
                ui.add_space(8.0);
                ui.label(text::ASK_LATER);
                ui.add_space(8.0);
                ui.button(text::CONTINUE)
                    .on_hover_text(text::CONTINUE_HINT)
                    .clicked()
            });
            // Escape or a click outside keeps the preselected choice.
            if modal.inner || modal.should_close() {
                g.ask = false;
            }
        }
        self.flow_window(ctx);
    }

    fn flow_window(&mut self, ctx: &egui::Context) {
        let Some((flow, i)) = self.guided.flow else {
            return;
        };
        let steps = flow.steps();
        let step = steps[i];
        let (title, what, waiting) = text::step(step);
        let window = flow.window(self).0;
        let mut open = true;
        let mut go = None;
        egui::Window::new(text::flow(flow).2)
            .tool_window(ctx, 6, [520.0, 480.0])
            .open(&mut open)
            .show(ctx, |ui| {
                ui.weak(text::step_of(i + 1, steps.len()));
                ui.heading(title);
                ui.label(what);
                ui.separator();
                // Per step, so a step does not open at the scroll offset of the one before.
                egui::ScrollArea::vertical()
                    .id_salt(step)
                    .max_height(ctx.content_rect().height() * 0.5)
                    .show(ui, |ui| self.step_ui(ui, step));
                ui.separator();
                ui.horizontal(|ui| {
                    if i > 0
                        && ui
                            .button(text::BACK)
                            .on_hover_text(text::BACK_HINT)
                            .clicked()
                    {
                        go = Some(i - 1);
                    }
                    if i + 1 == steps.len() {
                        if ui
                            .button(text::FINISH)
                            .on_hover_text(text::finish_hint(window))
                            .clicked()
                        {
                            go = Some(steps.len());
                        }
                        return;
                    }
                    if ui
                        .add_enabled(self.done(step), egui::Button::new(text::NEXT))
                        .on_hover_text(text::NEXT_HINT)
                        .on_disabled_hover_text(waiting)
                        .clicked()
                    {
                        go = Some(i + 1);
                    }
                    if ui
                        .button(text::SKIP)
                        .on_hover_text(text::SKIP_HINT)
                        .clicked()
                    {
                        go = Some(i + 1);
                    }
                    if ui
                        .button(text::TO_EXPERT)
                        .on_hover_text(text::TO_EXPERT_HINT)
                        .clicked()
                    {
                        go = Some(steps.len());
                    }
                });
            });
        self.guided.flow = match go {
            _ if !open => {
                self.cancel_flow();
                None
            }
            Some(n) if n == steps.len() => {
                self.restore_cmp();
                *flow.window(self).1 = true;
                None
            }
            Some(n) => Some((flow, n)),
            None => Some((flow, i)),
        };
    }

    /// Ends the running flow unfinished; a detection takes back the ".cmp" saving it switched on.
    pub fn cancel_flow(&mut self) {
        if self.guided.flow.take().is_some() {
            self.restore_cmp();
        }
    }

    /// The batch started in the last step has its options already; the shared setting goes back.
    fn restore_cmp(&mut self) {
        if std::mem::take(&mut self.guided.cmp_switched_on) {
            self.batch.stop_saving_cmp();
        }
    }

    /// Whether Next may go on; Skip always can.
    fn done(&self, step: Step) -> bool {
        match step {
            Step::OpenRun => !self.runs.is_empty(),
            Step::FindAlkanes => self.ri.cal.is_some() || self.ri.rt_only,
            Step::Library => self.library.is_some(),
            Step::BatchRuns => self.batch.has_runs(),
            Step::MidRuns => self.mid.problem().is_none(),
            _ => true,
        }
    }

    /// The expert window's building block for the step.
    fn step_ui(&mut self, ui: &mut egui::Ui, step: Step) {
        match step {
            Step::OpenRun => {
                self.open_chromatograms_button(ui);
                let name = self.runs.get(self.active).map(|r| r.name.as_str());
                ui.label(text::active_run(name));
            }
            Step::FindAlkanes => self.ri.contents(ui, self.runs.get(self.active)),
            Step::CheckPoints => {
                // A table edit that fails (times or RIs not rising) keeps the old calibration.
                self.ri.error_label(ui);
                self.ri.table(ui);
            }
            Step::Library => {
                self.open_library_button(ui);
                if self.library_job.is_some() {
                    ui.label(text::LIBRARY_LOADING);
                } else {
                    let name = self.library.as_ref().map(|l| l.name.as_str());
                    ui.label(text::library(name));
                }
            }
            Step::BatchRuns => self.batch.runs_row(ui),
            Step::SaveCmp => self.batch.cmp_grid(ui),
            Step::QuantSettings => self.batch.settings_grid(ui),
            Step::StartBatch => {
                let (cal, reference) = self.ri.for_search(self.library.as_deref());
                (self.batch).job_row(ui, self.library.as_ref(), cal, reference);
                if let Some(r) = self.batch.summary() {
                    ui.label(text::result(r));
                }
            }
            Step::MidMethod => self.mid.settings(ui),
            Step::MidRuns => {
                self.mid.run_lists(ui);
                // Also a tracer problem of the step before, which Next waits for too.
                if let Some(p) = self.mid.problem() {
                    ui.weak(p);
                }
            }
            Step::StartMid => {
                let (cal, reference) = self.ri.for_search(self.library.as_deref());
                (self.mid).start_row(ui, self.library.as_ref(), cal, reference);
                if let Some(r) = self.mid.summary() {
                    ui.label(text::result(r));
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mode_keys_round_trip() {
        for m in [Mode::Guided, Mode::Expert] {
            assert_eq!(Mode::from_key(m.key()), Some(m));
        }
        assert_eq!(Mode::from_key("wizard"), None);
    }

    #[test]
    fn result_names_runs_and_errors() {
        assert_eq!(
            text::result((3, 1)),
            "Last run: 3 runs, errors: 1. Finish shows the details."
        );
    }
}
