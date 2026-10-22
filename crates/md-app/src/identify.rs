//! Library search in the app: load a library, identify the components of every run on worker threads
//! (with RI once calibrated, shifted per run onto a reference substance if one is chosen), label hits on
//! the TIC and list the hits of the selected component. A second, NIST library is searched for the
//! selected component on request only (old "NIST search...").

use std::path::{Path, PathBuf};
use std::sync::Arc;

use eframe::egui;
use md_core::job::Cancel;
use md_core::library::LibraryCompound;
use md_core::ri::RiCalibration;
use md_core::search::{
    Identification, LibraryHit, PreparedLibrary, SearchSettings, library_hits,
    map_compounds_to_library, shift_to_reference,
};

use crate::lang::{t, tf};
use crate::settings::Pick;
use crate::{App, Job, Run};

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
pub fn score_color(score: f64) -> egui::Color32 {
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
    }) || md_core::library::is_targetsearch(path)
}

/// Reads and prepares a library. The read cannot be interrupted; `cancel` stops before the preparation.
fn load(path: &Path, s: SearchSettings, nist: bool, cancel: &Cancel) -> Result<Library, String> {
    let compounds = md_core::library::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if cancel.is_cancelled() {
        return Err("cancelled".into());
    }
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

/// Run path, component, its spectrum, the RI searched with, hits.
pub type NistHits = (
    PathBuf,
    usize,
    Vec<(u32, u32)>,
    Option<f64>,
    Vec<LibraryHit>,
);

/// Old "NIST search...": a second library next to the normal one (old setting `rc_lib_nist_file`),
/// searched only for the selected component and only on request. Always with the old fragment
/// pre-filter, like the old `SQLLibrary`.
#[derive(Default)]
pub struct Nist {
    pub lib: Option<Arc<Library>>,
    job: Option<LibraryJob>,
    /// Search settings whose preparation crashed; not retried until they change.
    failed: Option<SearchSettings>,
    /// The running search of the selected component.
    search: Option<Job<NistHits>>,
    /// The last search: run path, component, its spectrum, the RI searched with, hits.
    pub hits: Option<NistHits>,
}

impl Nist {
    pub fn busy(&self) -> bool {
        self.job.is_some() || self.search.is_some()
    }

    /// Hits of the last search if it was for this component of this run as it is now (a new
    /// analysis replaces the components).
    pub fn hits_for(&self, run: &Run, k: usize) -> Option<&[LibraryHit]> {
        let (path, c, spectrum, _, hits) = self.hits.as_ref()?;
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
            .button(t(
                "Open library (.mdlib, .lbr, .msl, NIST .msp/.sdf, TargetSearch .txt)…",
                "Bibliothek öffnen (.mdlib, .lbr, .msl, NIST .msp/.sdf, TargetSearch .txt)…",
            ))
            .on_hover_text(t(
                "Load a mass spectral library; components are then identified against it",
                "Lädt eine Massenspektren-Bibliothek; die Komponenten werden dann mit ihr identifiziert",
            ))
            .clicked()
        {
            self.pick_library(ui.ctx(), Pick::Library);
        }
    }

    /// Asks for a compound or NIST library file and opens it (menu buttons and Settings).
    pub fn pick_library(&mut self, ctx: &egui::Context, which: Pick) {
        let dialog = rfd::FileDialog::new();
        match which {
            Pick::Library => {
                let types = [
                    "mdlib", "lbr", "msl", "MSL", "msp", "MSP", "sdf", "SDF", "txt",
                ];
                if let Some(path) = dialog
                    .add_filter(t("Library", "Bibliothek"), &types)
                    .pick_file()
                {
                    self.open_library(ctx, path);
                }
            }
            Pick::Nist => {
                let types = ["sdf", "SDF", "msp", "MSP"];
                if let Some(path) = dialog
                    .add_filter(t("NIST text", "NIST-Text"), &types)
                    .pick_file()
                {
                    self.open_nist(ctx, path);
                }
            }
        }
    }

    pub fn open_library(&mut self, ctx: &egui::Context, path: PathBuf) {
        let s = self.settings.analysis.search.clone();
        self.library_job = Some(Job::spawn(ctx, move |_, c| {
            load(&path, s, false, c).map(Arc::new)
        }));
    }

    pub fn open_nist_button(&mut self, ui: &mut egui::Ui) {
        if ui
            .button(t("Open NIST library (.sdf, .msp)…", "NIST-Bibliothek öffnen (.sdf, .msp)…"))
            .on_hover_text(t(
                "Load a NIST library for \"NIST search\" on the selected component; the normal library stays",
                "Lädt eine NIST-Bibliothek für die \"NIST-Suche\" zur ausgewählten Komponente; die normale Bibliothek bleibt",
            ))
            .clicked()
        {
            self.pick_library(ui.ctx(), Pick::Nist);
        }
    }

    pub fn open_nist(&mut self, ctx: &egui::Context, path: PathBuf) {
        let s = self.settings.analysis.search.clone();
        self.nist.hits = None;
        self.nist.search = None;
        self.nist.job = Some(Job::spawn(ctx, move |_, c| {
            load(&path, s, true, c).map(Arc::new)
        }));
    }

    /// Takes a loaded NIST library and prepares it again when the search settings change.
    pub fn update_nist(&mut self, ctx: &egui::Context) {
        let wanted = &self.settings.analysis.search;
        if let Some(job) = &self.nist.job
            && let Some(result) = job.poll(|| t("NIST library", "NIST-Bibliothek").into())
        {
            match result {
                Ok(Ok(lib)) => {
                    self.nist.lib = Some(lib);
                    self.nist.failed = None;
                }
                Ok(Err(e)) => self.error = Some(e.into()),
                Err(e) => {
                    self.error = Some(e);
                    self.nist.failed = Some(wanted.clone());
                }
            }
            self.nist.job = None;
        }
        if let Some(lib) = &self.nist.lib
            && lib.settings != *wanted
            && self.nist.job.is_none()
            && self.nist.failed.as_ref() != Some(wanted)
        {
            let (lib, wanted) = (lib.clone(), wanted.clone());
            self.nist.hits = None;
            self.nist.search = None;
            self.nist.job = Some(Job::spawn(ctx, move |_, _| {
                Ok(Arc::new(reprepare(&lib, wanted, false, true)))
            }));
        }
        if let Some(job) = &self.nist.search
            && let Some(result) = job.poll(|| t("NIST search", "NIST-Suche").into())
        {
            match result {
                Ok(hits) => self.nist.hits = Some(hits),
                Err(e) => self.error = Some(e),
            }
            self.nist.search = None;
        }
        // A new RI calibration of that run makes the hits stale.
        if let Some((path, k, _, ri, _)) = &self.nist.hits
            && let Some(run) = self.runs.iter().find(|r| &r.path == path)
            && run.components.get(*k).is_some()
            && crate::results::component_ri(run, *k, self.ri.for_library(run)) != *ri
        {
            self.nist.hits = None;
        }
    }

    /// Old `nistSearch`: the selected component of the active run against the NIST library, with RI
    /// once the run is calibrated.
    pub fn nist_search(&mut self, ctx: &egui::Context) {
        let (Some(run), Some(lib)) = (self.runs.get(self.active), &self.nist.lib) else {
            return;
        };
        let Some(k) = run.component else {
            return;
        };
        let c = &run.components[k];
        let ri = crate::results::component_ri(run, k, self.ri.for_library(run));
        let (path, spectrum, lib) = (run.path.clone(), c.spectrum.clone(), lib.clone());
        self.nist.search = Some(Job::spawn(ctx, move |_, _| {
            let hits = library_hits(&spectrum, ri, &lib.prepared, &lib.settings);
            (path, k, spectrum, ri, hits)
        }));
    }

    /// The settings the library is prepared with. The fragment values steer only the NIST
    /// pre-filter, so changing them does not search a library without it again.
    fn wanted_search(&self) -> SearchSettings {
        let mut s = self.settings.analysis.search.clone();
        if let Some(lib) = &self.library
            && !md_core::library::is_nist(&lib.path)
        {
            s.max_fragments = lib.settings.max_fragments;
            s.min_equal_fragments = lib.settings.min_equal_fragments;
        }
        s
    }

    /// Takes a loaded library and (re)identifies every run that was not searched with the current
    /// library, RI calibration and reference.
    pub fn update_identification(&mut self, ctx: &egui::Context) {
        if let Some(job) = &self.library_job
            && let Some(result) = job.poll(|| t("library search", "Bibliothekssuche").into())
        {
            match result {
                Ok(Ok(lib)) => self.library = Some(lib),
                Ok(Err(e)) => self.error = Some(e.into()),
                Err(e) => {
                    self.error = Some(e);
                    self.library_failed = Some(self.wanted_search());
                }
            }
            self.library_job = None;
        }
        let Some(lib) = &self.library else {
            return;
        };
        // New search settings: prepare the library again, then every run is searched again.
        let wanted = self.wanted_search();
        let rt_only = self.ri.rt_only;
        if lib.settings != wanted || lib.prepared.is_rt_only() != rt_only {
            if self.library_job.is_none() && self.library_failed.as_ref() != Some(&wanted) {
                let lib = lib.clone();
                let filter = md_core::library::is_nist(&lib.path);
                self.library_job = Some(Job::spawn(ctx, move |_, _| {
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
            if let Some((k, job)) = &run.ident_job
                && let Some(result) = job.poll(|| {
                    let name = &run.name;
                    tf!("{name}: identification", "{name}: Identifizierung")
                })
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
                let job = Job::spawn(ctx, move |_, c| {
                    let (lib, times, s) =
                        (&k.0.prepared, &times.scan_times_s, k.0.settings.clone());
                    let shifted = k.1.as_ref().zip(k.2).and_then(|(c, (r, cutoff))| {
                        shift_to_reference(c, times, &comps, lib, r, cutoff, &s)
                    });
                    let cal = shifted.as_ref().or(k.1.as_deref());
                    let ris = cal.map(|c| c.component_ris(times, &comps));
                    let id = map_compounds_to_library(&comps, ris.as_deref(), lib, &s, c)
                        .unwrap_or_default();
                    (id, shifted)
                });
                run.ident_job = Some((key, job));
            }
        }
    }
}

/// Names of the assigned library compounds at the component apexes on the TIC, coloured by score.
pub fn labels(run: &Run, show: impl Fn(f64) -> bool) -> Vec<Label> {
    let Some(((lib, ..), id, _)) = &run.ident else {
        return Vec::new();
    };
    md_core::view::hit_labels(
        &run.components,
        &id.assigned,
        &lib.compounds,
        &run.minutes,
        &run.tic,
        show,
    )
    .into_iter()
    .map(|l| (l.at, l.name, score_color(l.score), l.score))
    .collect()
}

/// Best hits of the selected component, the assigned one marked.
pub fn hits_list(ui: &mut egui::Ui, run: &Run) {
    let (Some(((lib, ..), id, _)), Some(k)) = (&run.ident, run.component) else {
        return;
    };
    let name = &lib.name;
    ui.strong(tf!("Library hits ({name})", "Bibliothekstreffer ({name})"));
    // A new search runs after a library or calibration change; until then these scores are old.
    if run.ident_job.is_some() {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label(t("Searching again…", "Erneute Suche…"));
        });
    }
    let assigned = id.assigned[k].map(|h| h.compound);
    egui::Grid::new("hit table")
        .striped(true)
        .show(ui, |ui| hit_rows(ui, &id.hits[k], lib, assigned));
    if id.hits[k].is_empty() {
        ui.label(t("No hit.", "Kein Treffer."));
    } else if let Some(note) = assignment_note(&id.hits[k], assigned) {
        ui.weak(note);
    }
}

/// Why the peak has no name, or a name other than the first hit (#242): each library compound goes
/// to one component only, the one with the highest score (ties: the lower index).
fn assignment_note(hits: &[LibraryHit], assigned: Option<usize>) -> Option<&'static str> {
    match assigned {
        None => Some(t(
            "Not assigned: each compound went to a component that scores as high or higher.",
            "Nicht zugeordnet: Jede Substanz ging an eine Komponente mit gleich hohem oder höherem Score.",
        )),
        Some(a) if hits.first().is_some_and(|h| h.compound != a) => Some(t(
            "Assigned to a lower hit: the better ones went to components that score as high or higher.",
            "Einem niedrigeren Treffer zugeordnet: Die besseren gingen an Komponenten mit gleich hohem oder höherem Score.",
        )),
        Some(_) => None,
    }
}

/// Hit table rows: score, spectrum score, RI score, name (bold for the assigned one).
fn hit_rows(ui: &mut egui::Ui, hits: &[LibraryHit], lib: &Library, assigned: Option<usize>) {
    for head in ["Score", t("Spectrum", "Spektrum"), "RI", "Name"] {
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
            ui.strong(name).on_hover_text(t(
                "Assigned to this component",
                "Dieser Komponente zugeordnet",
            ));
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
        .add_enabled(!nist.busy(), egui::Button::new(t("NIST search", "NIST-Suche")))
        .on_hover_text({
            let name = &lib.name;
            tf!(
                "Search this component in the NIST library {name}; the library hits above stay",
                "Sucht diese Komponente in der NIST-Bibliothek {name}; die Bibliothekstreffer oben bleiben",
            )
        })
        .on_disabled_hover_text(t(
            "The NIST library is being prepared or searched",
            "Die NIST-Bibliothek wird vorbereitet oder durchsucht",
        ))
        .clicked();
    if let Some(hits) = nist.hits_for(run, k) {
        let name = &lib.name;
        ui.strong(tf!("NIST hits ({name})", "NIST-Treffer ({name})"));
        egui::Grid::new("nist hit table")
            .striped(true)
            .show(ui, |ui| hit_rows(ui, hits, lib, None));
        if hits.is_empty() {
            ui.label(t("No hit.", "Kein Treffer."));
        }
    }
    clicked
}

pub type LibraryJob = Job<Result<Arc<Library>, String>>;
pub type IdentJob = (Key, Job<(Identification, Option<RiCalibration>)>);

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cancelled_library_load_stops_before_preparing() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../md-core/tests/data/synthetic-nist.sdf");
        let s = App::default().settings.applied.analysis().search;
        let cancel = Cancel::default();
        cancel.cancel();
        assert_eq!(
            load(&path, s, true, &cancel).err().as_deref(),
            Some("cancelled")
        );
    }

    /// NIST search fills its own hit list for the selected component; the library hits stay.
    #[test]
    fn nist_search_keeps_the_library_hits() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../md-core/tests/data/synthetic-nist.sdf");
        let mut app = App::default();
        let nist = load(
            &path,
            app.settings.analysis.search.clone(),
            true,
            &Cancel::default(),
        )
        .unwrap();
        let alpha = &nist.compounds[0].spectrum;
        let mut run = crate::results::tests::run();
        run.component = Some(0);
        run.components[0].spectrum = (alpha.masses.iter().zip(&alpha.values))
            .map(|(m, v)| (*m as u32, *v as u32))
            .collect();
        let before = run.ident.as_ref().map(|i| i.1.hits.clone());
        app.runs.push(run);
        app.nist.lib = Some(Arc::new(nist));
        let ctx = egui::Context::default();
        app.nist_search(&ctx);
        while app.nist.busy() {
            app.update_nist(&ctx);
            std::thread::yield_now();
        }
        let hits = app.nist.hits_for(&app.runs[0], 0).unwrap();
        let lib = app.nist.lib.as_ref().unwrap();
        assert_eq!(lib.compounds[hits[0].compound].name, "Synthetic alpha");
        assert_eq!(app.runs[0].ident.as_ref().map(|i| i.1.hits.clone()), before);
        // Another component, or new components after a re-analysis, do not show these hits.
        assert!(app.nist.hits_for(&app.runs[0], 1).is_none());
        app.runs[0].components[0].spectrum.pop();
        assert!(app.nist.hits_for(&app.runs[0], 0).is_none());

        // A new calibration alone keeps the hits: the list shows the stored identification's RI
        // until a new identification ends. Once it does, hits searched with the old RI are dropped.
        app.runs[0].components[0]
            .spectrum
            .clone_from(&app.nist.hits.as_ref().unwrap().2);
        app.update_nist(&ctx);
        assert!(app.nist.hits_for(&app.runs[0], 0).is_some());
        app.ri.cal = Some(Arc::new(RiCalibration::identity()));
        app.update_nist(&ctx);
        assert!(app.nist.hits.is_some());
        app.runs[0].ident.as_mut().unwrap().0.1 = Some(Arc::new(RiCalibration::identity()));
        app.update_nist(&ctx);
        assert!(app.nist.hits.is_none());

        // Settings that change while a search runs drop its hits too.
        app.nist_search(&ctx);
        app.settings
            .edit_applied(|s| s.sql_minimal_number_equal_fragments = 2);
        app.update_nist(&ctx);
        while app.nist.busy() {
            app.update_nist(&ctx);
            std::thread::yield_now();
        }
        assert!(app.nist.hits.is_none());

        // The fragment values matter only for a library with the pre-filter.
        let lib = app.nist.lib.clone().unwrap();
        app.settings
            .edit_applied(|s| s.sql_minimal_number_equal_fragments = 1);
        app.library = Some(Arc::new(Library {
            path: "other.msl".into(),
            ..reprepare(&lib, lib.settings.clone(), false, false)
        }));
        assert_eq!(app.wanted_search(), lib.settings);
        app.library = Some(lib.clone());
        assert_eq!(app.wanted_search().min_equal_fragments, 1);
    }

    /// Changing LS_CUTOFF prepares the library again and identifies the runs again (#399): a hit
    /// below the new cutoff is gone, one above it is back.
    #[test]
    fn a_new_ls_cutoff_identifies_again() {
        let path =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../md-core/tests/data/synthetic-nist.sdf");
        let mut app = App::default();
        app.settings.edit_applied(|s| s.ls_cutoff = 0.0);
        // The pre-filter of a NIST library needs shared fragments; any one will do.
        app.settings
            .edit_applied(|s| s.sql_minimal_number_equal_fragments = 0);
        let lib = load(
            &path,
            app.settings.analysis.search.clone(),
            false,
            &Cancel::default(),
        )
        .unwrap();
        let alpha = &lib.compounds[0].spectrum;
        let mut run = crate::results::tests::run();
        // One ion missing: a good hit, but below a cutoff of 1.
        run.components[0].spectrum = (alpha.masses.iter().zip(&alpha.values))
            .map(|(m, v)| (*m as u32, *v as u32))
            .skip(1)
            .collect();
        app.runs.push(run);
        app.library = Some(Arc::new(lib));
        let ctx = egui::Context::default();
        // Wait for both jobs to finish, with a deadline instead of a poll count.
        let settle = |app: &mut App| {
            let deadline = std::time::Instant::now() + std::time::Duration::from_secs(60);
            loop {
                app.update_identification(&ctx);
                if app.library_job.is_none() && app.runs[0].ident_job.is_none() {
                    break;
                }
                assert!(
                    std::time::Instant::now() < deadline,
                    "identification jobs did not finish"
                );
                std::thread::yield_now();
            }
        };
        settle(&mut app);
        let assigned = |app: &App| app.runs[0].ident.as_ref().unwrap().1.assigned[0];
        let score = assigned(&app).expect("assigned at cutoff 0").overall;
        assert!((0.0..1.0).contains(&score));
        app.settings
            .edit_applied(|s| s.ls_cutoff = (score + 1.0) / 2.0);
        settle(&mut app);
        assert!(assigned(&app).is_none() && app.error.is_none());
        app.settings.edit_applied(|s| s.ls_cutoff = 0.0);
        settle(&mut app);
        assert_eq!(assigned(&app).map(|h| h.overall), Some(score));
    }

    #[test]
    fn note_explains_an_assignment_other_than_the_first_hit() {
        let hit = |compound| LibraryHit {
            compound,
            spec_sim: 0.5,
            ri_sim: None,
            overall: 0.5,
        };
        let hits = [hit(3), hit(7)];
        assert_eq!(assignment_note(&hits, Some(3)), None);
        assert!(
            assignment_note(&hits, Some(7))
                .unwrap()
                .starts_with("Assigned")
        );
        assert!(
            assignment_note(&hits, None)
                .unwrap()
                .starts_with("Not assigned")
        );
    }
}
