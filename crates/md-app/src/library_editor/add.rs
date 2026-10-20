//! "Add component" and "Copy retention times" of the library editor (F18; old `DialogAddLib` and
//! `MainWindow::copyRT`, docs/reverse-engineering/algorithms.md §7e). Both work on the library being
//! edited; saving stays with the editor.

use std::sync::{Arc, mpsc};

use crate::layout::ToolWindow;
use eframe::egui;
use md_core::Chromatogram;
use md_core::batch::{IonPeak, QuantSettings, ion_peaks_by};
use md_core::deconv::Component;
use md_core::library::LibraryCompound;
use md_core::library_add::{
    AddMode, RtCopy, add, apply_rt_copy, compound_from_component, find, plan_rt_copy,
};
use md_core::ri::RiCalibration;
use md_core::search::SearchSettings;

use super::LibraryEditor;
use crate::{Run, poll_job, spawn};

/// What the editor takes from the active run and the settings.
pub struct Source<'a> {
    pub run: &'a Run,
    /// The run's RI calibration, if any.
    pub ri: Option<&'a RiCalibration>,
    pub quant: QuantSettings,
    pub search: SearchSettings,
}

impl Source<'_> {
    pub fn component(&self) -> Option<&Component> {
        self.run.components.get(self.run.component?)
    }
}

/// The old radio buttons for a name that is already in the library.
#[derive(Clone, Copy, PartialEq)]
enum Choice {
    Replace,
    KeepFeatures,
    QuantIons,
}

pub struct AddDialog {
    comp: Component,
    data: Arc<Chromatogram>,
    ri: Option<f64>,
    quant: QuantSettings,
    name: String,
    comment: String,
    choice: Choice,
    /// The question "Do you really want to …?" is showing.
    confirm: bool,
    /// The ion peaks behind the component, largest integral first.
    ions: Vec<IonPeak>,
}

pub struct RtDialog {
    /// Searching the components against the library.
    job: Option<mpsc::Receiver<Vec<RtCopy>>>,
    plan: Vec<RtCopy>,
}

impl AddDialog {
    /// Name, the old radio buttons for a name already in the library, and the comment.
    fn choices(&mut self, ui: &mut egui::Ui, library: &[LibraryCompound], exists: bool) {
        let before = (self.name.clone(), self.choice);
        ui.horizontal(|ui| {
            ui.label("Name");
            ui.text_edit_singleline(&mut self.name);
        });
        egui::ComboBox::from_label("Existing entry")
            .selected_text("choose…")
            .height(300.0)
            .show_ui(ui, |ui| {
                for c in library {
                    if ui.selectable_label(c.name == self.name, &c.name).clicked() {
                        self.name.clone_from(&c.name);
                    }
                }
            });
        if exists {
            ui.separator();
            ui.label("The library has an entry with this name:");
            ui.radio_value(&mut self.choice, Choice::Replace, "Replace entry")
                .on_hover_text("Remove the entry and add the found compound.");
            ui.radio_value(&mut self.choice, Choice::KeepFeatures, "Replace RT/Spec")
                .on_hover_text(
                    "Replace spectrum, retention time and retention index but keep quantification ions and comments.",
                );
            ui.radio_value(&mut self.choice, Choice::QuantIons, "Exchange quant. ions")
                .on_hover_text("Set only the automatic quantification ions of the entry.");
        } else {
            ui.label("A new entry is added.");
        }
        ui.label("Comment");
        ui.add(egui::TextEdit::multiline(&mut self.comment).desired_rows(2));
        // The question was for the entry and choice it showed.
        if (self.name.as_str(), self.choice) != (before.0.as_str(), before.1) {
            self.confirm = false;
        }
    }

    /// RT and RI of the component and the question.
    fn status(&self, ui: &mut egui::Ui) {
        ui.label(format!(
            "RT {:.3} min, RI {}",
            self.data.scan_times_s[self.comp.apex_scan as usize] / 60.0,
            self.ri.map_or("unknown".into(), |r| format!("{r:.1}")),
        ));
        if self.confirm {
            ui.colored_label(
                ui.visuals().warn_fg_color,
                match self.choice {
                    Choice::Replace => "Do you really want to replace the existing entry?",
                    Choice::KeepFeatures => {
                        "Do you really want to replace spectrum, RT and RI of the entry?"
                    }
                    Choice::QuantIons => {
                        "Do you really want to set the detected quantification ions?"
                    }
                },
            );
        }
    }

    /// OK (asks first for an existing name) and Cancel; returns (apply, cancel).
    fn buttons(&mut self, ui: &mut egui::Ui, exists: bool) -> (bool, bool) {
        let (mut apply, mut cancel) = (false, false);
        ui.horizontal(|ui| {
            let ready = !self.name.trim().is_empty();
            let label = if self.confirm { "Yes" } else { "OK" };
            let ok = ui
                .add_enabled(ready, egui::Button::new(label))
                .on_disabled_hover_text("Needs a name.");
            if ok.clicked() {
                if exists && !self.confirm {
                    self.confirm = true;
                } else {
                    apply = true;
                }
            }
            cancel = ui.button("Cancel").clicked();
        });
        (apply, cancel)
    }
}

impl LibraryEditor {
    /// Opens "Add component" for the selected component of the active run.
    pub(super) fn start_add(&mut self, src: &Source) {
        let Some(comp) = src.component() else { return };
        let run = src.run;
        let rt_s = run.data.scan_times_s[comp.apex_scan as usize];
        // From the peaks the components were built from: also right for a `.cmp` made with another
        // deconvolution width, which a fresh deconvolution would not find again (#312).
        let ions = ion_peaks_by(comp, &run.peaks, |i| run.peak_masses[i]);
        self.add = Some(AddDialog {
            comp: comp.clone(),
            data: run.data.clone(),
            ri: src.ri.map(|cal| cal.ri(rt_s)),
            quant: src.quant.clone(),
            name: format!("{} {:.3} min", run.name, rt_s / 60.0),
            comment: String::new(),
            choice: Choice::Replace,
            confirm: false,
            ions,
        });
    }

    /// Opens "Copy retention times": the components of the active run against the edited library.
    pub(super) fn start_rt_copy(&mut self, ctx: &egui::Context, src: &Source) {
        let run = src.run;
        let (comps, data) = (run.components.clone(), run.data.clone());
        let ris = src
            .ri
            .map(|cal| cal.component_ris(&data.scan_times_s, &comps));
        let (target, s) = (self.compounds.clone(), src.search.clone());
        let job = spawn(ctx, move || {
            plan_rt_copy(&comps, &data.scan_times_s, ris.as_deref(), &target, &s)
        });
        self.rt = Some(RtDialog {
            job: Some(job),
            plan: Vec::new(),
        });
    }

    pub(super) fn add_windows(&mut self, ctx: &egui::Context) {
        self.add_window(ctx);
        self.rt_window(ctx);
    }

    fn add_window(&mut self, ctx: &egui::Context) {
        let Some(d) = &mut self.add else { return };
        let existing = find(&self.compounds, &d.name);
        let (mut open, mut apply, mut cancel) = (true, false, false);
        egui::Window::new("Add component to library")
            .tool_window(ctx, 9, [340.0, 420.0])
            .open(&mut open)
            .collapsible(false)
            .show(ctx, |ui| {
                d.choices(ui, &self.compounds, existing.is_some());
                d.status(ui);
                ui.separator();
                (apply, cancel) = d.buttons(ui, existing.is_some());
            });
        if !open || cancel {
            self.add = None;
        } else if apply {
            self.apply_add(existing);
        }
    }

    fn apply_add(&mut self, existing: Option<usize>) {
        let Some(d) = self.add.take() else { return };
        let c = compound_from_component(
            &d.name,
            &d.comp,
            &d.data.scan_times_s,
            d.ri,
            &d.ions,
            &d.quant,
        );
        let mode = match (existing, d.choice) {
            (None, _) => AddMode::New,
            (Some(i), Choice::Replace) => AddMode::Replace(i),
            (Some(i), Choice::KeepFeatures) => AddMode::KeepFeatures(i),
            (Some(i), Choice::QuantIons) => AddMode::QuantIons(i),
        };
        match add(&mut self.compounds, c, mode, &d.comment) {
            Ok(i) => {
                self.dirty = true;
                self.select(Some(i));
            }
            Err(e) => self.error = Some(e.to_string()),
        }
    }

    fn rt_window(&mut self, ctx: &egui::Context) {
        let Some(d) = &mut self.rt else { return };
        if let Some(rx) = &d.job
            && let Some(r) = poll_job(rx, "retention time search")
        {
            d.job = None;
            match r {
                Ok(plan) => d.plan = plan,
                Err(e) => self.error = Some(e),
            }
        }
        let (mut open, mut apply) = (true, false);
        egui::Window::new("Copy retention times")
            .tool_window(ctx, 10, [340.0, 420.0])
            .open(&mut open)
            .show(ctx, |ui| {
                ui.label(format!(
                    "Assigned components of the active run → entries of {}.",
                    if self.name.is_empty() {
                        "the new library"
                    } else {
                        &self.name
                    }
                ));
                if d.job.is_some() {
                    ui.horizontal(|ui| {
                        ui.spinner();
                        ui.label("Searching the library…");
                    });
                } else if d.plan.is_empty() {
                    ui.label("No component was assigned to an entry.");
                }
                egui::ScrollArea::vertical()
                    .max_height(360.0)
                    .show(ui, |ui| {
                        egui::Grid::new("rt plan").striped(true).show(ui, |ui| {
                            for h in [
                                "Compound",
                                "RT (min)",
                                "Lib RT (min)",
                                "RI",
                                "Lib RI",
                                "Score",
                            ] {
                                ui.strong(h);
                            }
                            ui.end_row();
                            for p in &d.plan {
                                ui.label(&p.name);
                                ui.label(rt_minutes(p.rt_ms));
                                ui.label(rt_minutes(p.old_rt_ms));
                                ui.label(p.ri.map_or("–".into(), ri_text));
                                ui.label(ri_text(p.old_ri));
                                ui.label(format!("{:.3}", p.score));
                                ui.end_row();
                            }
                        });
                    });
                if ui
                    .add_enabled(
                        !d.plan.is_empty(),
                        egui::Button::new(format!("Copy {} RTs", d.plan.len())),
                    )
                    .on_hover_text("Replaces RT and RI of the entries; Save writes the library.")
                    .clicked()
                {
                    apply = true;
                }
            });
        if apply {
            let plan = std::mem::take(&mut d.plan);
            match apply_rt_copy(&mut self.compounds, &plan) {
                Ok(()) => self.dirty = true,
                Err(e) => self.error = Some(e.to_string()),
            }
            open = false;
        }
        if !open {
            self.rt = None;
        }
    }
}

/// RI as shown in the copy dialog; a negative value is unknown.
fn ri_text(v: f64) -> String {
    if v >= 0.0 {
        format!("{v:.1}")
    } else {
        "–".into()
    }
}

/// RT in minutes from milliseconds; a negative value is unknown.
fn rt_minutes(v: f64) -> String {
    if v >= 0.0 {
        format!("{:.3}", v / 60_000.0)
    } else {
        "–".into()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The ions come from the peaks behind the components, not from a new deconvolution with the
    /// current settings: a `.cmp` made with another width has to work too (#312).
    #[test]
    fn add_takes_the_ions_of_the_stored_peaks() {
        let mut run = crate::results::tests::run();
        run.component = Some(1);
        run.analysis = crate::settings::Analysis::default();
        let src = Source {
            run: &run,
            ri: None,
            quant: QuantSettings::default(),
            search: SearchSettings::default(),
        };
        let mut ed = LibraryEditor::default();
        ed.start_add(&src);
        assert!(ed.error.is_none());
        let ions = &ed.add.as_ref().expect("dialog").ions;
        assert_eq!(
            ions,
            &[IonPeak {
                mass: 73,
                integral: 5.0,
                quality: 0.5,
                signal_to_noise: 12.0
            }]
        );
    }
}
