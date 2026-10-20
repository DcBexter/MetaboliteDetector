//! Library search in the app: load a library, identify the components of every run on worker threads
//! (with RI once calibrated, shifted per run onto a reference substance if one is chosen), label hits on
//! the TIC and list the hits of the selected component.

use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};

use eframe::egui;
use md_core::library::{LibraryCompound, import_msl, read_lbr};
use md_core::ri::RiCalibration;
use md_core::search::{
    Identification, PreparedLibrary, SearchSettings, map_compounds_to_library, shift_to_reference,
};

use crate::{App, Run, poll_job, spawn};

pub struct Library {
    pub name: String,
    pub compounds: Vec<LibraryCompound>,
    pub prepared: PreparedLibrary,
    /// What `prepared` was scaled and filtered with; searches use the same.
    pub settings: SearchSettings,
}

/// Library, RI calibration and RT shift reference (library index, cutoff) an identification was
/// computed with.
pub type Key = (
    Arc<Library>,
    Option<Arc<RiCalibration>>,
    Option<(usize, f64)>,
);

/// Identification of one run, what it was computed with, and the calibration shifted onto the
/// reference when it was found.
pub type Identified = (Key, Identification, Option<RiCalibration>);

/// A text on the TIC: position, name, colour.
pub type Label = ([f64; 2], String, egui::Color32);

/// Old label colours: above 0.9, above 0.7, the rest.
fn score_color(score: f64) -> egui::Color32 {
    match score {
        s if s > 0.9 => egui::Color32::from_rgb(0, 150, 0),
        s if s > 0.7 => egui::Color32::from_rgb(220, 140, 0),
        _ => egui::Color32::from_rgb(210, 40, 40),
    }
}

/// Libraries are recognised by extension; everything else opens as a run.
pub fn is_library(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("lbr") || e.eq_ignore_ascii_case("msl"))
}

fn load(path: &Path, s: SearchSettings) -> Result<Library, String> {
    let lbr = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("lbr"));
    let compounds = if lbr {
        read_lbr(path)
    } else {
        import_msl(path)
    }
    .map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Library {
        name: path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into(),
        prepared: PreparedLibrary::new(&compounds, &s),
        settings: s,
        compounds,
    })
}

impl App {
    pub fn library_menu(&mut self, ui: &mut egui::Ui) {
        if ui.button("Open library (.lbr, .msl)…").clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("Library", &["lbr", "msl", "MSL"])
                .pick_file()
        {
            self.open_library(ui.ctx(), path);
        }
    }

    pub fn open_library(&mut self, ctx: &egui::Context, path: PathBuf) {
        let s = self.settings.applied.search();
        self.library_job = Some(spawn(ctx, move || load(&path, s).map(Arc::new)));
    }

    /// Takes a loaded library and (re)identifies every run that was not searched with the current
    /// library, RI calibration and reference.
    pub fn update_identification(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.library_job
            && let Some(result) = poll_job(rx, "library search")
        {
            match result {
                Ok(Ok(lib)) => self.library = Some(lib),
                Ok(Err(e)) => self.error = Some(e),
                Err(e) => {
                    self.error = Some(e);
                    self.library_failed = Some(self.settings.applied.search());
                }
            }
            self.library_job = None;
        }
        let Some(lib) = &self.library else {
            return;
        };
        // New search settings: prepare the library again, then every run is searched again.
        let wanted = self.settings.applied.search();
        if lib.settings != wanted {
            if self.library_job.is_none() && self.library_failed.as_ref() != Some(&wanted) {
                let lib = lib.clone();
                self.library_job = Some(spawn(ctx, move || {
                    Ok(Arc::new(Library {
                        name: lib.name.clone(),
                        prepared: PreparedLibrary::new(&lib.compounds, &wanted),
                        compounds: lib.compounds.clone(),
                        settings: wanted,
                    }))
                }));
            }
            return;
        }
        let cal = &self.ri.cal;
        let reference = (self.ri.reference.as_ref())
            .filter(|_| cal.is_some())
            // The window offers only compounds with RI; a namesake without RI must not shadow them.
            .and_then(|n| (lib.compounds.iter()).position(|c| &c.name == n && c.spectrum.ri >= 0.0))
            .map(|k| (k, self.ri.ref_cutoff));
        let current = |k: &Key| {
            Arc::ptr_eq(&k.0, lib)
                && k.1.as_ref().map(Arc::as_ptr) == cal.as_ref().map(Arc::as_ptr)
                && k.2 == reference
        };
        for run in &mut self.runs {
            if let Some((k, rx)) = &run.ident_job
                && let Some(result) = poll_job(rx, &format!("{}: identification", run.name))
            {
                // A crashed search is not retried for this key; the run shows no hits.
                let (id, shifted) = result.unwrap_or_else(|e| {
                    self.error = Some(e);
                    let n = run.components.len();
                    let none = Identification {
                        hits: vec![Vec::new(); n],
                        assigned: vec![None; n],
                    };
                    (none, None)
                });
                run.ident = Some((k.clone(), id, shifted));
                run.ident_job = None;
            }
            if !run.ident.as_ref().is_some_and(|i| current(&i.0))
                && !run.ident_job.as_ref().is_some_and(|j| current(&j.0))
            {
                let key: Key = (lib.clone(), cal.clone(), reference);
                let (k, comps, times) = (key.clone(), run.components.clone(), run.data.clone());
                let rx = spawn(ctx, move || {
                    let (lib, times, s) =
                        (&k.0.prepared, &times.scan_times_s, k.0.settings.clone());
                    let shifted = k.1.as_ref().zip(k.2).and_then(|(c, (r, cutoff))| {
                        shift_to_reference(c, times, &comps, lib, r, cutoff, &s)
                    });
                    let cal = shifted.as_ref().or(k.1.as_deref());
                    let ris = cal.map(|c| c.component_ris(times, &comps));
                    let id = map_compounds_to_library(&comps, ris.as_deref(), lib, &s);
                    (id, shifted)
                });
                run.ident_job = Some((key, rx));
            }
        }
    }
}

/// Names of the assigned library compounds at the component apexes on the TIC.
pub fn labels(run: &Run) -> Vec<Label> {
    let Some(((lib, ..), id, _)) = &run.ident else {
        return Vec::new();
    };
    run.components
        .iter()
        .zip(&id.assigned)
        .filter_map(|(k, hit)| {
            let h = hit.as_ref()?;
            let s = k.apex_scan as usize;
            Some((
                [run.minutes[s], run.tic[s]],
                lib.compounds[h.compound].name.clone(),
                score_color(h.overall),
            ))
        })
        .collect()
}

/// Best hits of the selected component, the assigned one marked.
pub fn hits_panel(ui: &mut egui::Ui, run: &Run) {
    let (Some(((lib, ..), id, _)), Some(k)) = (&run.ident, run.component) else {
        return;
    };
    egui::Panel::right("hits").show(ui, |ui| {
        ui.strong(format!("Library hits ({})", lib.name));
        // A new search runs after a library or calibration change; until then these scores are old.
        if run.ident_job.is_some() {
            ui.horizontal(|ui| {
                ui.spinner();
                ui.label("Searching again…");
            });
        }
        let assigned = id.assigned[k].map(|h| h.compound);
        egui::Grid::new("hit table").striped(true).show(ui, |ui| {
            ui.label("Score");
            ui.label("Spectrum");
            ui.label("RI");
            ui.label("Name");
            ui.end_row();
            for h in id.hits[k].iter().take(15) {
                ui.colored_label(score_color(h.overall), format!("{:.3}", h.overall));
                ui.label(format!("{:.3}", h.spec_sim));
                // Without calibration (or library RI) the overall score is the spectrum score alone.
                ui.label(h.ri_sim.map_or("–".into(), |r| format!("{r:.3}")));
                let name = &lib.compounds[h.compound].name;
                if assigned == Some(h.compound) {
                    ui.strong(name).on_hover_text("Assigned to this component");
                } else {
                    ui.label(name);
                }
                ui.end_row();
            }
        });
        if id.hits[k].is_empty() {
            ui.label("No hit.");
        }
    });
}

pub type LibraryJob = mpsc::Receiver<Result<Arc<Library>, String>>;
pub type IdentJob = (Key, mpsc::Receiver<(Identification, Option<RiCalibration>)>);
