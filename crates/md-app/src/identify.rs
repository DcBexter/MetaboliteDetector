//! Library search in the app: load a library, identify the components of every run on worker threads
//! (with RI once calibrated, shifted per run onto a reference substance if one is chosen), label hits on
//! the TIC and list the hits of the selected component. A second, NIST library is searched for the
//! selected component on request only (old "NIST search...").

use std::path::{Path, PathBuf};
use std::sync::{Arc, mpsc};

use eframe::egui;
use md_core::library::LibraryCompound;
use md_core::ri::RiCalibration;
use md_core::search::{
    Identification, LibraryHit, PreparedLibrary, SearchSettings, library_hits,
    map_compounds_to_library, shift_to_reference,
};

use crate::{App, Run, poll_job, spawn};

pub struct Library {
    pub path: PathBuf,
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

/// A text on the TIC: position, name, colour, library score.
pub type Label = ([f64; 2], String, egui::Color32, f64);

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
    path.extension().is_some_and(|e| {
        ["lbr", "msl", "msp", "sdf", "mdlib"]
            .iter()
            .any(|x| e.eq_ignore_ascii_case(x))
    })
}

fn load(path: &Path, s: SearchSettings, nist: bool) -> Result<Library, String> {
    let compounds = md_core::library::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    Ok(Library {
        path: path.to_path_buf(),
        name: path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into(),
        prepared: prepare(
            &compounds,
            &s,
            false,
            nist || md_core::library::is_nist(path),
        ),
        settings: s,
        compounds,
    })
}

fn prepare(
    c: &[LibraryCompound],
    s: &SearchSettings,
    rt_only: bool,
    filter: bool,
) -> PreparedLibrary {
    match rt_only {
        true => PreparedLibrary::new_rt_only(c, s),
        false => PreparedLibrary::new(c, s),
    }
    .with_fragment_filter(filter)
}

/// The same library prepared again with new search settings.
fn reprepare(lib: &Library, s: SearchSettings, rt_only: bool, filter: bool) -> Library {
    Library {
        path: lib.path.clone(),
        name: lib.name.clone(),
        prepared: prepare(&lib.compounds, &s, rt_only, filter),
        compounds: lib.compounds.clone(),
        settings: s,
    }
}

/// Run path, component, its spectrum, hits.
pub type NistHits = (PathBuf, usize, Vec<(u32, u32)>, Vec<LibraryHit>);

/// Old "NIST search...": a second library next to the normal one (old setting `rc_lib_nist_file`),
/// searched only for the selected component and only on request. Always with the old fragment
/// pre-filter, like the old `SQLLibrary`.
#[derive(Default)]
pub struct Nist {
    pub lib: Option<Arc<Library>>,
    job: Option<LibraryJob>,
    /// Search settings whose preparation crashed; not retried until they change.
    failed: Option<SearchSettings>,
    /// The last search: run path, component, its spectrum, hits.
    pub hits: Option<NistHits>,
}

impl Nist {
    pub fn busy(&self) -> bool {
        self.job.is_some()
    }

    /// Hits of the last search if it was for this component of this run as it is now (a new
    /// analysis replaces the components).
    pub fn hits_for(&self, run: &Run, k: usize) -> Option<&[LibraryHit]> {
        let (path, c, spectrum, hits) = self.hits.as_ref()?;
        let same = run
            .components
            .get(k)
            .is_some_and(|x| &x.spectrum == spectrum);
        (path == &run.path && *c == k && same).then_some(hits)
    }
}

impl App {
    pub fn open_library_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .button("Open library (.mdlib, .lbr, .msl, NIST .msp/.sdf)…")
            .on_hover_text(
                "Load a mass spectral library; components are then identified against it",
            )
            .clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter(
                    "Library",
                    &["mdlib", "lbr", "msl", "MSL", "msp", "MSP", "sdf", "SDF"],
                )
                .pick_file()
        {
            self.open_library(ui.ctx(), path);
        }
    }

    pub fn open_library(&mut self, ctx: &egui::Context, path: PathBuf) {
        let s = self.settings.applied.search();
        self.library_job = Some(spawn(ctx, move || load(&path, s, false).map(Arc::new)));
    }

    pub fn open_nist_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .button("Open NIST library (.sdf, .msp)…")
            .on_hover_text(
                "Load a NIST library for \"NIST search\" on the selected component; the normal library stays",
            )
            .clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("NIST text", &["sdf", "SDF", "msp", "MSP"])
                .pick_file()
        {
            let s = self.settings.applied.search();
            self.nist.hits = None;
            self.nist.job = Some(spawn(ui.ctx(), move || load(&path, s, true).map(Arc::new)));
        }
    }

    /// Takes a loaded NIST library and prepares it again when the search settings change.
    pub fn update_nist(&mut self, ctx: &egui::Context) {
        let wanted = self.settings.applied.search();
        if let Some(rx) = &self.nist.job
            && let Some(result) = poll_job(rx, "NIST library")
        {
            match result {
                Ok(Ok(lib)) => self.nist.lib = Some(lib),
                Ok(Err(e)) => self.error = Some(e),
                Err(e) => {
                    self.error = Some(e);
                    self.nist.failed = Some(wanted.clone());
                }
            }
            self.nist.job = None;
        }
        if let Some(lib) = &self.nist.lib
            && lib.settings != wanted
            && self.nist.job.is_none()
            && self.nist.failed.as_ref() != Some(&wanted)
        {
            let lib = lib.clone();
            self.nist.hits = None;
            self.nist.job = Some(spawn(ctx, move || {
                Ok(Arc::new(reprepare(&lib, wanted, false, true)))
            }));
        }
    }

    /// Old `nistSearch`: the selected component of the active run against the NIST library, with RI
    /// once the run is calibrated.
    pub fn nist_search(&mut self) {
        let (Some(run), Some(lib)) = (self.runs.get(self.active), &self.nist.lib) else {
            return;
        };
        let Some(k) = run.component else {
            return;
        };
        let c = &run.components[k];
        let ri =
            (self.ri.for_run(run)).map(|cal| cal.ri(run.data.scan_times_s[c.apex_scan as usize]));
        // shortcut: on the UI thread; with the pre-filter a 300 000 entry library takes ~10 ms.
        let hits = library_hits(&c.spectrum, ri, &lib.prepared, &lib.settings);
        self.nist.hits = Some((run.path.clone(), k, c.spectrum.clone(), hits));
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
        let rt_only = self.ri.rt_only;
        if lib.settings != wanted || lib.prepared.is_rt_only() != rt_only {
            if self.library_job.is_none() && self.library_failed.as_ref() != Some(&wanted) {
                let lib = lib.clone();
                let filter = md_core::library::is_nist(&lib.path);
                self.library_job = Some(spawn(ctx, move || {
                    Ok(Arc::new(reprepare(&lib, wanted, rt_only, filter)))
                }));
            }
            return;
        }
        let reference = (self.ri.reference.as_ref())
            // The window offers only compounds with RI; a namesake without RI must not shadow them.
            .and_then(|n| {
                (lib.compounds.iter()).position(|c| &c.name == n && self.ri.anchor(c) >= 0.0)
            })
            .map(|k| (k, self.ri.ref_cutoff));
        for run in &mut self.runs {
            let cal = self.ri.base(run).cloned();
            let reference = reference.filter(|_| cal.is_some());
            let current = |k: &Key| {
                Arc::ptr_eq(&k.0, lib)
                    && k.1.as_ref().map(Arc::as_ptr) == cal.as_ref().map(Arc::as_ptr)
                    && k.2 == reference
            };
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
                run.rev += 1;
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
pub fn labels(run: &Run, show: impl Fn(f64) -> bool) -> Vec<Label> {
    let Some(((lib, ..), id, _)) = &run.ident else {
        return Vec::new();
    };
    run.components
        .iter()
        .zip(&id.assigned)
        .filter_map(|(k, hit)| {
            let h = hit.as_ref().filter(|h| show(h.overall))?;
            let s = k.apex_scan as usize;
            Some((
                [run.minutes[s], run.tic[s]],
                lib.compounds[h.compound].name.clone(),
                score_color(h.overall),
                h.overall,
            ))
        })
        .collect()
}

/// Best hits of the selected component, the assigned one marked.
pub fn hits_list(ui: &mut egui::Ui, run: &Run) {
    let (Some(((lib, ..), id, _)), Some(k)) = (&run.ident, run.component) else {
        return;
    };
    ui.strong(format!("Library hits ({})", lib.name));
    // A new search runs after a library or calibration change; until then these scores are old.
    if run.ident_job.is_some() {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label("Searching again…");
        });
    }
    let assigned = id.assigned[k].map(|h| h.compound);
    egui::Grid::new("hit table")
        .striped(true)
        .show(ui, |ui| hit_rows(ui, &id.hits[k], lib, assigned));
    if id.hits[k].is_empty() {
        ui.label("No hit.");
    }
}

/// Hit table rows: score, spectrum score, RI score, name (bold for the assigned one).
fn hit_rows(ui: &mut egui::Ui, hits: &[LibraryHit], lib: &Library, assigned: Option<usize>) {
    for head in ["Score", "Spectrum", "RI", "Name"] {
        ui.strong(head);
    }
    ui.end_row();
    for h in hits.iter().take(15) {
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
}

/// "NIST search" for the selected component and its hits; true when the search was asked for.
pub fn nist_list(ui: &mut egui::Ui, run: &Run, k: usize, nist: &Nist) -> bool {
    let Some(lib) = &nist.lib else {
        return false;
    };
    let clicked = ui
        .add_enabled(!nist.busy(), egui::Button::new("NIST search"))
        .on_hover_text(format!(
            "Search this component in the NIST library {}; the library hits above stay",
            lib.name
        ))
        .on_disabled_hover_text("The NIST library is being prepared")
        .clicked();
    if let Some(hits) = nist.hits_for(run, k) {
        ui.strong(format!("NIST hits ({})", lib.name));
        egui::Grid::new("nist hit table")
            .striped(true)
            .show(ui, |ui| hit_rows(ui, hits, lib, None));
        if hits.is_empty() {
            ui.label("No hit.");
        }
    }
    clicked
}

pub type LibraryJob = mpsc::Receiver<Result<Arc<Library>, String>>;
pub type IdentJob = (Key, mpsc::Receiver<(Identification, Option<RiCalibration>)>);

/// Which labels stay visible when they are drawn left-aligned at `starts` (pixels) with the
/// given `widths`: the highest score wins, a label overlapping a kept one is hidden.
pub fn visible_labels(starts: &[f32], widths: &[f32], scores: &[f64]) -> Vec<bool> {
    let mut order: Vec<usize> = (0..starts.len()).collect();
    order.sort_by(|a, b| scores[*b].total_cmp(&scores[*a]));
    let mut keep = vec![false; starts.len()];
    let mut kept: Vec<usize> = Vec::new();
    for i in order {
        if kept
            .iter()
            .all(|k| starts[i] >= starts[*k] + widths[*k] || starts[*k] >= starts[i] + widths[i])
        {
            keep[i] = true;
            kept.push(i);
        }
    }
    keep
}

#[cfg(test)]
mod tests {
    use super::*;

    /// NIST search fills its own hit list for the selected component; the library hits stay.
    #[test]
    fn nist_search_keeps_the_library_hits() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../md-core/tests/data/synthetic-nist.sdf");
        let s = SearchSettings {
            cutoff: 0.0,
            ..Default::default()
        };
        let nist = load(&path, s, true).unwrap();
        let alpha = &nist.compounds[0].spectrum;
        let mut run = crate::results::tests::run();
        run.component = Some(0);
        run.components[0].spectrum = (alpha.masses.iter().zip(&alpha.values))
            .map(|(m, v)| (*m as u32, *v as u32))
            .collect();
        let before = run.ident.as_ref().map(|i| i.1.hits.clone());
        let mut app = App::default();
        app.runs.push(run);
        app.nist.lib = Some(Arc::new(nist));
        app.nist_search();
        let hits = app.nist.hits_for(&app.runs[0], 0).unwrap();
        let lib = app.nist.lib.as_ref().unwrap();
        assert_eq!(lib.compounds[hits[0].compound].name, "Synthetic alpha");
        assert_eq!(app.runs[0].ident.as_ref().map(|i| i.1.hits.clone()), before);
        // Another component, or new components after a re-analysis, do not show these hits.
        assert!(app.nist.hits_for(&app.runs[0], 1).is_none());
        app.runs[0].components[0].spectrum.pop();
        assert!(app.nist.hits_for(&app.runs[0], 0).is_none());
    }

    #[test]
    fn overlapping_labels_keep_the_best_score() {
        // 0 and 1 overlap (1 scores higher), 2 is clear of both.
        let keep = visible_labels(&[0.0, 30.0, 100.0], &[50.0, 50.0, 40.0], &[0.8, 0.95, 0.5]);
        assert_eq!(keep, [false, true, true]);
        // Touching labels do not overlap.
        assert_eq!(
            visible_labels(&[0.0, 50.0], &[50.0, 50.0], &[0.9, 0.8]),
            [true, true]
        );
    }
}
