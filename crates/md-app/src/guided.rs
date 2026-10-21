//! Guided mode (S7, #256): a global switch "Guided / Expert", asked at the first start, and four step
//! flows (RI calibration, component detection, batch quantification, MID). Each step draws a building
//! block of the expert window on the same state, so a flow ends in that window with its values set.
//! In Guided mode the expert options of the windows sit in a closed "Advanced" section.

use eframe::egui;

use crate::App;
use crate::lang::t;
use crate::layout::ToolWindow;

/// Every text this module shows, in one place for the translation (S8, #257).
pub mod text {
    use super::{Flow, Step};
    use crate::lang::{t, tf};

    pub fn advanced() -> &'static str {
        t("Advanced", "Erweitert")
    }
    pub fn switch() -> &'static str {
        t("Guided mode", "Geführter Modus")
    }
    pub fn switch_hint() -> &'static str {
        t(
            "Guided: the main tasks step by step, expert options under a \
             closed \"Advanced\". Expert: every option open. Every function stays reachable in both.",
            "Geführt: die Hauptaufgaben Schritt für Schritt, Expertenoptionen unter einem \
             geschlossenen „Erweitert“. Experte: alle Optionen offen. Jede Funktion bleibt in beiden \
             erreichbar.",
        )
    }
    pub fn mode() -> &'static str {
        t("Mode", "Modus")
    }
    pub fn ask_title() -> &'static str {
        t("How do you want to work?", "Wie möchten Sie arbeiten?")
    }
    pub fn guided() -> &'static str {
        t("Guided", "Geführt")
    }
    pub fn guided_hint() -> &'static str {
        t(
            "Recommended to start: the main tasks step by step \
             (Analysis → Step by step), expert options under \"Advanced\".",
            "Für den Einstieg empfohlen: die Hauptaufgaben Schritt für Schritt \
             (Analyse → Schritt für Schritt), Expertenoptionen unter „Erweitert“.",
        )
    }
    pub fn expert() -> &'static str {
        t("Expert", "Experte")
    }
    pub fn expert_hint() -> &'static str {
        t(
            "Every window with all its options open, as in the old program.",
            "Jedes Fenster mit allen Optionen offen, wie im alten Programm.",
        )
    }
    pub fn ask_later() -> &'static str {
        t(
            "You can switch at any time in Help and in Settings.",
            "Sie können jederzeit unter Hilfe und in den Einstellungen wechseln.",
        )
    }
    pub fn continue_() -> &'static str {
        t("Continue", "Weiter")
    }
    pub fn menu() -> &'static str {
        t("Step by step", "Schritt für Schritt")
    }
    pub fn back() -> &'static str {
        t("Back", "Zurück")
    }
    pub fn next() -> &'static str {
        t("Next", "Weiter")
    }
    pub fn skip() -> &'static str {
        t("Skip", "Überspringen")
    }
    pub fn skip_hint() -> &'static str {
        t(
            "Go on without this step; you can come back with Back",
            "Ohne diesen Schritt weitermachen; mit Zurück kommen Sie wieder hierher",
        )
    }
    pub fn finish() -> &'static str {
        t("Finish", "Fertigstellen")
    }
    pub fn to_expert() -> &'static str {
        t("Open the window", "Fenster öffnen")
    }
    pub fn to_expert_hint() -> &'static str {
        t(
            "Leave the steps and go on in the full window, with the values set so far",
            "Die Schritte verlassen und im vollen Fenster weitermachen, mit den bisher gesetzten Werten",
        )
    }
    pub fn library_loading() -> &'static str {
        t("Library: loading…", "Bibliothek: wird geladen…")
    }
    pub fn back_hint() -> &'static str {
        t(
            "The step before; what you set stays",
            "Der Schritt davor; Ihre Eingaben bleiben",
        )
    }
    pub fn next_hint() -> &'static str {
        t("The next step", "Der nächste Schritt")
    }
    pub fn continue_hint() -> &'static str {
        t(
            "Start working in the chosen mode",
            "Im gewählten Modus beginnen",
        )
    }

    pub fn active_run(name: Option<&str>) -> String {
        match name {
            Some(n) => tf!("Active run: {n}", "Aktiver Lauf: {n}"),
            None => t("No run open yet", "Noch kein Lauf geöffnet").into(),
        }
    }

    pub fn library(name: Option<&str>) -> String {
        match name {
            Some(n) => tf!("Library: {n}", "Bibliothek: {n}"),
            None => t("No library open yet", "Noch keine Bibliothek geöffnet").into(),
        }
    }

    pub fn step_of(x: usize, n: usize) -> String {
        tf!("Step {x} of {n}", "Schritt {x} von {n}")
    }

    pub fn result((runs, errors): (usize, usize)) -> String {
        tf!(
            "Last run: {runs} runs, errors: {errors}. Finish shows the details.",
            "Letzter Durchgang: {runs} Läufe, Fehler: {errors}. Fertigstellen zeigt die Details."
        )
    }

    pub fn finish_hint(window: &str) -> String {
        tf!(
            "Close the steps and open the window \"{window}\" with the values set",
            "Die Schritte schließen und das Fenster „{window}“ mit den gesetzten Werten öffnen"
        )
    }

    /// Menu entry, its hover text and the title of the steps' window.
    pub fn flow(f: Flow) -> (&'static str, &'static str, &'static str) {
        match f {
            Flow::Ri => (
                t("RI calibration…", "RI-Kalibrierung…"),
                t(
                    "Retention indices from an alkane mix or a .calib file",
                    "Retentionsindizes aus einem Alkan-Mix oder einer .calib-Datei",
                ),
                t(
                    "RI calibration step by step",
                    "RI-Kalibrierung Schritt für Schritt",
                ),
            ),
            Flow::Detection => (
                t("Component detection…", "Komponentenerkennung…"),
                t(
                    "Detect the components of many runs and save them as .cmp files",
                    "Die Komponenten vieler Läufe erkennen und als .cmp-Dateien speichern",
                ),
                t(
                    "Component detection step by step",
                    "Komponentenerkennung Schritt für Schritt",
                ),
            ),
            Flow::Batch => (
                t("Batch quantification…", "Batch-Quantifizierung…"),
                t(
                    "Quantify the library's substances across many runs",
                    "Die Substanzen der Bibliothek über viele Läufe quantifizieren",
                ),
                t(
                    "Batch quantification step by step",
                    "Batch-Quantifizierung Schritt für Schritt",
                ),
            ),
            Flow::Mid => (
                t("Mass isotopomer distributions…", "Isotopomere (MID)…"),
                t(
                    "Isotopomer distributions of the library's fragments in labelling experiments",
                    "Isotopomerverteilungen der Fragmente der Bibliothek in Markierungsexperimenten",
                ),
                t("MIDs step by step", "MIDs Schritt für Schritt"),
            ),
        }
    }

    /// The English title of the steps' window: its id, stable across a language switch.
    pub fn flow_id(f: Flow) -> &'static str {
        match f {
            Flow::Ri => "RI calibration step by step",
            Flow::Detection => "Component detection step by step",
            Flow::Batch => "Batch quantification step by step",
            Flow::Mid => "MIDs step by step",
        }
    }

    /// Title, what to do, and why Next waits.
    pub fn step(s: Step) -> (&'static str, &'static str, &'static str) {
        match s {
            Step::OpenRun => (
                t("Open the alkane run", "Den Alkan-Lauf öffnen"),
                t(
                    "Open the run of your n-alkane mix. It becomes the active run.",
                    "Öffnen Sie den Lauf Ihres n-Alkan-Mixes. Er wird der aktive Lauf.",
                ),
                t("Open a run first", "Erst einen Lauf öffnen"),
            ),
            Step::FindAlkanes => (
                t("Find the alkanes", "Die Alkane finden"),
                t(
                    "Set the marker ion and the first and last alkane of your mix, then find them \
                     in the active run. A saved .calib file works too.",
                    "Setzen Sie das Marker-Ion und das erste und letzte Alkan Ihres Mixes und \
                     finden Sie sie dann im aktiven Lauf. Eine gespeicherte .calib-Datei geht auch.",
                ),
                t(
                    "Find the alkanes or load a .calib file first",
                    "Erst die Alkane finden oder eine .calib-Datei laden",
                ),
            ),
            Step::CheckPoints => (
                t("Check the points", "Die Punkte prüfen"),
                t(
                    "Each alkane gives one point from retention time to retention index. Switch \
                     off a point that was found at the wrong peak.",
                    "Jedes Alkan ergibt einen Punkt von Retentionszeit zu Retentionsindex. \
                     Schalten Sie einen Punkt ab, der am falschen Peak gefunden wurde.",
                ),
                "",
            ),
            Step::Library => (
                t("Open the library", "Die Bibliothek öffnen"),
                t(
                    "The substances are searched in this library. It is opened again at the next start.",
                    "Die Substanzen werden in dieser Bibliothek gesucht. Sie wird beim nächsten \
                     Start wieder geöffnet.",
                ),
                t("Open a library first", "Erst eine Bibliothek öffnen"),
            ),
            Step::BatchRuns => (
                t("Add the runs", "Die Läufe hinzufügen"),
                t(
                    "Add every run of the experiment.",
                    "Fügen Sie jeden Lauf des Experiments hinzu.",
                ),
                t("Add runs first", "Erst Läufe hinzufügen"),
            ),
            Step::SaveCmp => (
                t("Where the components go", "Wohin die Komponenten kommen"),
                t(
                    "The components of each run are saved as run.cmp, next to the run or in a \
                     folder you choose. They load faster later.",
                    "Die Komponenten jedes Laufs werden als run.cmp gespeichert, neben dem Lauf \
                     oder in einem Ordner Ihrer Wahl. Sie laden später schneller.",
                ),
                "",
            ),
            Step::QuantSettings => (
                t("Groups and filters", "Gruppen und Filter"),
                t(
                    "Name the replicate groups by a pattern on the file name or a sample table, and \
                     how reproducible a substance must be.",
                    "Benennen Sie die Replikatgruppen über ein Muster im Dateinamen oder eine \
                     Probentabelle, und legen Sie fest, wie reproduzierbar eine Substanz sein muss.",
                ),
                "",
            ),
            Step::StartBatch => (
                t("Start", "Starten"),
                t(
                    "Start the scan. It goes on in the background; the window shows the result.",
                    "Den Durchgang starten. Er läuft im Hintergrund; das Fenster zeigt das Ergebnis.",
                ),
                "",
            ),
            Step::MidMethod => (
                t("Correction", "Korrektur"),
                t(
                    "Choose how the labelled spectra are corrected and the tracer element.",
                    "Wählen Sie, wie die markierten Spektren korrigiert werden, und das \
                     Tracer-Element.",
                ),
                "",
            ),
            Step::MidRuns => (
                t("Add the runs", "Die Läufe hinzufügen"),
                t(
                    "Add the labelled runs and, for the reference correction, unlabelled runs of the \
                     same substances.",
                    "Fügen Sie die markierten Läufe hinzu und, für die Referenzkorrektur, \
                     unmarkierte Läufe derselben Substanzen.",
                ),
                t(
                    "Add the runs first, or fix what the step shows",
                    "Erst die Läufe hinzufügen oder beheben, was der Schritt zeigt",
                ),
            ),
            Step::StartMid => (
                t("Start", "Starten"),
                t(
                    "Start the calculation. It goes on in the background; the window shows the result.",
                    "Die Berechnung starten. Sie läuft im Hintergrund; das Fenster zeigt das \
                     Ergebnis.",
                ),
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
            Flow::Ri => (t("RI calibration", "RI-Kalibrierung"), &mut app.ri.open),
            Flow::Detection | Flow::Batch => (
                t("Batch quantification", "Batch-Quantifizierung"),
                &mut app.batch.open,
            ),
            Flow::Mid => (
                t(
                    "Mass isotopomer distributions",
                    "Massenisotopomerverteilungen",
                ),
                &mut app.mid.open,
            ),
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
    advanced_open(ui, id, false, add);
}

/// Like [`advanced`]; `force` keeps the section open, so a value to fix is never hidden.
pub fn advanced_open(ui: &mut egui::Ui, id: &str, force: bool, add: impl FnOnce(&mut egui::Ui)) {
    let mode = ui
        .ctx()
        .data(|d| d.get_temp::<Mode>(mode_id()))
        .unwrap_or_default();
    // Salted with the mode, so a switch opens or closes every section anew.
    egui::CollapsingHeader::new(text::advanced())
        .id_salt((id, mode))
        .default_open(mode == Mode::Expert)
        .open(force.then_some(true))
        .show(ui, add);
}

/// The switch in Help and in Settings.
pub fn switch(ui: &mut egui::Ui, mode: &mut Mode) {
    let mut on = *mode == Mode::Guided;
    if ui
        .checkbox(&mut on, text::switch())
        .on_hover_text(text::switch_hint())
        .changed()
    {
        *mode = if on { Mode::Guided } else { Mode::Expert };
    }
}

impl App {
    /// Analysis → Step by step: the four flows, in both modes.
    pub fn flows_menu(&mut self, ui: &mut egui::Ui) {
        ui.menu_button(text::menu(), |ui| {
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
                ui.heading(text::ask_title());
                ui.add_space(8.0);
                for (mode, name, hint) in [
                    (Mode::Guided, text::guided(), text::guided_hint()),
                    (Mode::Expert, text::expert(), text::expert_hint()),
                ] {
                    ui.radio_value(&mut g.mode, mode, name);
                    ui.indent(name, |ui| ui.weak(hint));
                }
                ui.add_space(8.0);
                ui.label(text::ask_later());
                ui.add_space(8.0);
                ui.button(text::continue_())
                    .on_hover_text(text::continue_hint())
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
            .id(egui::Id::unique(text::flow_id(flow)))
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
                            .button(text::back())
                            .on_hover_text(text::back_hint())
                            .clicked()
                    {
                        go = Some(i - 1);
                    }
                    if i + 1 == steps.len() {
                        if ui
                            .button(text::finish())
                            .on_hover_text(text::finish_hint(window))
                            .clicked()
                        {
                            go = Some(steps.len());
                        }
                        return;
                    }
                    if ui
                        .add_enabled(self.done(step), egui::Button::new(text::next()))
                        .on_hover_text(text::next_hint())
                        .on_disabled_hover_text(waiting)
                        .clicked()
                    {
                        go = Some(i + 1);
                    }
                    if ui
                        .button(text::skip())
                        .on_hover_text(text::skip_hint())
                        .clicked()
                    {
                        go = Some(i + 1);
                    }
                    if ui
                        .button(text::to_expert())
                        .on_hover_text(text::to_expert_hint())
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
                    ui.label(text::library_loading());
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
