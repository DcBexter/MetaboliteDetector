//! Growing a library from measurements: add a found component to a library (old `DialogAddLib`) and
//! copy the retention times of found substances into a target library (old `MainWindow::copyRT`).
//! Behaviour and what is not ported: `docs/reverse-engineering/algorithms.md` §7e.
//!
//! Pure functions on `Vec<LibraryCompound>`; saving goes through `library::save_mdlib`.

use crate::batch::{IonPeak, QuantSettings, quant_ions};
use crate::deconv::Component;
use crate::library::{LibraryCompound, LibrarySpectrum};
use crate::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};

/// Source name the old dialog gives a compound it adds (`addCompound(c, "no source")`).
const NO_SOURCE: &str = "no source";

/// The four choices of the old dialog (radio buttons `rb_replace`, `rb_keep`, `rb_quant_ions`, and
/// adding under a name not in the library). The index is the entry of the library.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AddMode {
    /// A new entry; the name must not be in the library.
    New,
    /// "Replace entry": the entry is removed and the found compound takes its place.
    Replace(usize),
    /// "Replace RT/Spec": spectrum, RT and RI are replaced; quant ions, features and comment stay.
    KeepFeatures(usize),
    /// "Exchange quant. ions": only the automatic quant ions of the entry are set.
    QuantIons(usize),
}

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum AddError {
    #[error("the library has no entry {0}")]
    NoEntry(usize),
    #[error("the library already has an entry named {0:?}")]
    NameTaken(String),
    #[error("the compound needs a name")]
    NoName,
}

/// Index of the entry named `name` (the old `LibrarySearch::contains` compares names exactly).
pub fn find(library: &[LibraryCompound], name: &str) -> Option<usize> {
    library.iter().position(|c| c.name == name)
}

/// The old `LibraryCompound(Compound, source)` for a found component: its deconvoluted spectrum,
/// retention time at the apex scan, RI (`None` = not calibrated, stored as -1 like the old program)
/// and the automatic quant ions of §7.2. `ions` are the component's ion peaks, largest integral first.
pub fn compound_from_component(
    name: &str,
    comp: &Component,
    scan_times_s: &[f64],
    ri: Option<f64>,
    ions: &[IonPeak],
    quant: &QuantSettings,
) -> LibraryCompound {
    let mut c = LibraryCompound::new(name);
    c.spectrum = LibrarySpectrum {
        scan: comp.apex_scan,
        rt_ms: scan_times_s[comp.apex_scan as usize] * 1000.0,
        ri: ri.unwrap_or(-1.0),
        bin_width: 1.0,
        masses: comp.spectrum.iter().map(|&(m, _)| m as f32).collect(),
        values: comp
            .spectrum
            .iter()
            .map(|&(_, h)| i32::try_from(h).unwrap_or(i32::MAX))
            .collect(),
    };
    c.auto_quant_ions = quant_ions(ions, quant)
        .into_iter()
        .map(|m| m as f32)
        .collect();
    set_combined(&mut c);
    c.sources = vec![NO_SOURCE.into()];
    c
}

/// Combined spectrum and ion frequency of a single spectrum, like the MSL importer leaves them.
fn set_combined(c: &mut LibraryCompound) {
    let s = &c.spectrum;
    c.comb_spec = s
        .masses
        .iter()
        .copied()
        .zip(s.values.iter().copied())
        .collect();
    c.ion_freq = s.masses.iter().map(|&m| (m, 1)).collect();
    c.spec_count = 1;
}

/// Puts `found` into `library` as `mode` says and returns the index of the changed entry. `comment`
/// goes into the feature `COMMENT` when it is not empty. In the modes that name an entry, the entry
/// keeps its name. The library is only changed when the call succeeds.
pub fn add(
    library: &mut Vec<LibraryCompound>,
    mut found: LibraryCompound,
    mode: AddMode,
    comment: &str,
) -> Result<usize, AddError> {
    let entry = |i: usize, library: &[LibraryCompound]| {
        library
            .get(i)
            .map(|c| c.name.clone())
            .ok_or(AddError::NoEntry(i))
    };
    match mode {
        AddMode::New => {
            if found.name.trim().is_empty() {
                return Err(AddError::NoName);
            }
            if find(library, &found.name).is_some() {
                return Err(AddError::NameTaken(found.name));
            }
            set_comment(&mut found, comment);
            library.push(found);
            Ok(library.len() - 1)
        }
        AddMode::Replace(i) => {
            found.name = entry(i, library)?;
            set_comment(&mut found, comment);
            library[i] = found;
            Ok(i)
        }
        AddMode::KeepFeatures(i) => {
            found.name = entry(i, library)?;
            let old = &library[i];
            // The old dialog copies the entry's features over the new ones, which drops a comment typed
            // in the dialog; here the typed comment is appended to the old one.
            found.quant_ions.clone_from(&old.quant_ions);
            found.auto_quant_ions.clone_from(&old.auto_quant_ions);
            found.features.clone_from(&old.features);
            set_comment(&mut found, comment);
            library[i] = found;
            Ok(i)
        }
        AddMode::QuantIons(i) => {
            entry(i, library)?;
            library[i].auto_quant_ions = found.auto_quant_ions;
            Ok(i)
        }
    }
}

fn set_comment(c: &mut LibraryCompound, comment: &str) {
    if !comment.is_empty() {
        c.features
            .entry("COMMENT".into())
            .and_modify(|old| {
                old.push('\n');
                old.push_str(comment);
            })
            .or_insert_with(|| comment.into());
    }
}

/// One line of the preview of "Copy retention times to target lib": a found component and the entry
/// of the target library it was assigned to.
#[derive(Debug, Clone, PartialEq)]
pub struct RtCopy {
    /// Index into the components given to [`plan_rt_copy`].
    pub component: usize,
    /// Index into the target library.
    pub target: usize,
    pub name: String,
    pub score: f64,
    /// Retention time of the component at its apex, in ms.
    pub rt_ms: f64,
    /// RI of the component; `None` while the run is not calibrated (the entry keeps its RI then).
    pub ri: Option<f64>,
    /// What the target entry holds now.
    pub old_rt_ms: f64,
    pub old_ri: f64,
}

/// Searches every component against `target` (old `CompoundMapper::mapCompsToLib`, one-to-one by best
/// score) and lists the assigned ones, ordered like the components. `ris[i]` is the RI of component `i`.
pub fn plan_rt_copy(
    comps: &[Component],
    scan_times_s: &[f64],
    ris: Option<&[f64]>,
    target: &[LibraryCompound],
    s: &SearchSettings,
) -> Vec<RtCopy> {
    let id = map_compounds_to_library(comps, ris, &PreparedLibrary::new(target, s), s);
    id.assigned
        .iter()
        .enumerate()
        .filter_map(|(i, hit)| {
            let hit = hit.as_ref()?;
            let t = &target[hit.compound];
            Some(RtCopy {
                component: i,
                target: hit.compound,
                name: t.name.clone(),
                score: hit.overall,
                rt_ms: scan_times_s[comps[i].apex_scan as usize] * 1000.0,
                ri: ris.and_then(|r| r.get(i).copied()),
                old_rt_ms: t.spectrum.rt_ms,
                old_ri: t.spectrum.ri,
            })
        })
        .collect()
}

/// Writes the retention times (and RIs, where known) of `plan` into `target`. Fails without changing
/// anything when a line points outside the library or at an entry of another name (the library
/// changed since the plan).
pub fn apply_rt_copy(target: &mut [LibraryCompound], plan: &[RtCopy]) -> Result<(), AddError> {
    if let Some(p) = plan
        .iter()
        .find(|p| target.get(p.target).is_none_or(|t| t.name != p.name))
    {
        return Err(AddError::NoEntry(p.target));
    }
    for p in plan {
        let s = &mut target[p.target].spectrum;
        s.rt_ms = p.rt_ms;
        if let Some(ri) = p.ri {
            s.ri = ri;
        }
    }
    Ok(())
}

/// Why a library entry is not usable as it stands. The editor still saves such entries (the user may
/// be half way through), but names them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Issue {
    /// No spectrum peak with an intensity above zero.
    NoSpectrum,
    /// A quant ion outside the mass range of the spectrum.
    QuantIonOutOfRange(f32),
    /// A quant ion inside the mass range that has no peak in the spectrum.
    QuantIonNotInSpectrum(f32),
}

impl std::fmt::Display for Issue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NoSpectrum => write!(f, "the spectrum has no peaks"),
            Self::QuantIonOutOfRange(m) => write!(f, "quant ion {m} is outside the mass range"),
            Self::QuantIonNotInSpectrum(m) => write!(f, "quant ion {m} is not in the spectrum"),
        }
    }
}

/// The problems of one entry; empty when it is fine. Without a spectrum the quant ions are not
/// checked, the missing spectrum is the one problem.
pub fn check(c: &LibraryCompound) -> Vec<Issue> {
    let s = &c.spectrum;
    let peaks: Vec<f32> = s
        .masses
        .iter()
        .zip(&s.values)
        .filter(|(_, v)| **v > 0)
        .map(|(m, _)| *m)
        .collect();
    let (Some(&lo), Some(&hi)) = (
        peaks.iter().min_by(|a, b| a.total_cmp(b)),
        peaks.iter().max_by(|a, b| a.total_cmp(b)),
    ) else {
        return vec![Issue::NoSpectrum];
    };
    let half = s.bin_width.max(1.0) / 2.0;
    c.quant_ions
        .iter()
        .filter_map(|&q| {
            if !q.is_finite() || q < lo - half || q > hi + half {
                Some(Issue::QuantIonOutOfRange(q))
            } else if peaks.iter().all(|m| (m - q).abs() > half) {
                Some(Issue::QuantIonNotInSpectrum(q))
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deconv::{DeconvSettings, deconvolute};
    use crate::library::{parse_mdlib, write_mdlib};
    use crate::peaks::PeakSettings;

    fn comp(spectrum: &[(u32, u32)], apex: u32) -> Component {
        Component {
            first_scan: 0,
            apex_scan: apex,
            last_scan: 4,
            model: Vec::new(),
            spectrum: spectrum.to_vec(),
            peaks: Vec::new(),
            exact_spectrum: Vec::new(),
        }
    }

    fn ion(mass: u32, integral: f64) -> IonPeak {
        IonPeak {
            mass,
            integral,
            quality: 0.1,
            signal_to_noise: 10.0,
        }
    }

    const TIMES: [f64; 5] = [1.0, 2.0, 3.0, 4.0, 5.0];

    fn found(name: &str) -> LibraryCompound {
        let ions = [ion(73, 90.0), ion(147, 50.0), ion(45, 5.0)];
        compound_from_component(
            name,
            &comp(&[(45, 100), (73, 900), (147, 400)], 2),
            &TIMES,
            Some(1234.5),
            &ions,
            &QuantSettings::default(),
        )
    }

    #[test]
    fn component_becomes_a_compound() {
        let c = found("X");
        assert_eq!(c.spectrum.masses, vec![45.0, 73.0, 147.0]);
        assert_eq!(c.spectrum.values, vec![100, 900, 400]);
        assert_eq!(c.spectrum.rt_ms, 3000.0);
        assert_eq!(c.spectrum.ri, 1234.5);
        assert_eq!(c.spectrum.scan, 2);
        assert_eq!(c.auto_quant_ions, vec![45.0, 73.0, 147.0]);
        assert!(c.quant_ions.is_empty());
        assert_eq!(c.comb_spec, vec![(45.0, 100), (73.0, 900), (147.0, 400)]);
        assert_eq!(c.spec_count, 1);
        let none = compound_from_component(
            "Y",
            &comp(&[(45, 1)], 0),
            &TIMES,
            None,
            &[],
            &QuantSettings::default(),
        );
        assert_eq!(none.spectrum.ri, -1.0);
        assert!(none.auto_quant_ions.is_empty());
    }

    #[test]
    fn new_entry_needs_a_free_name_and_survives_a_save() {
        let mut lib = vec![LibraryCompound::new("A")];
        assert_eq!(
            add(&mut lib, found("A"), AddMode::New, ""),
            Err(AddError::NameTaken("A".into()))
        );
        assert_eq!(
            add(&mut lib, found(" "), AddMode::New, ""),
            Err(AddError::NoName)
        );
        assert_eq!(lib.len(), 1);
        assert_eq!(add(&mut lib, found("B"), AddMode::New, "from run 1"), Ok(1));
        assert_eq!(lib[1].features["COMMENT"], "from run 1");
        assert_eq!(parse_mdlib(&write_mdlib(&lib)).unwrap(), lib);
    }

    #[test]
    fn replace_takes_the_name_of_the_entry() {
        let mut old = found("Old");
        old.quant_ions = vec![73.0];
        let mut lib = vec![LibraryCompound::new("A"), old];
        assert_eq!(
            add(&mut lib, found("typed"), AddMode::Replace(1), ""),
            Ok(1)
        );
        assert_eq!(lib[1].name, "Old");
        assert!(lib[1].quant_ions.is_empty());
        assert_eq!(
            add(&mut lib, found("z"), AddMode::Replace(2), ""),
            Err(AddError::NoEntry(2))
        );
    }

    #[test]
    fn keep_features_keeps_ions_and_comment() {
        let mut old = LibraryCompound::new("Old");
        old.quant_ions = vec![73.0];
        old.auto_quant_ions = vec![45.0];
        old.features.insert("COMMENT".into(), "checked".into());
        old.features.insert("FRAG".into(), "(1/2)".into());
        old.spectrum.rt_ms = 1.0;
        let mut lib = vec![old];
        add(&mut lib, found("x"), AddMode::KeepFeatures(0), "").unwrap();
        let c = &lib[0];
        assert_eq!(
            (c.name.as_str(), c.spectrum.rt_ms, c.spectrum.ri),
            ("Old", 3000.0, 1234.5)
        );
        assert_eq!(c.spectrum.masses.len(), 3);
        assert_eq!(c.quant_ions, vec![73.0]);
        assert_eq!(c.auto_quant_ions, vec![45.0]);
        assert_eq!(c.features["COMMENT"], "checked");
        assert_eq!(c.features["FRAG"], "(1/2)");
        add(&mut lib, found("x"), AddMode::KeepFeatures(0), "new note").unwrap();
        assert_eq!(lib[0].features["COMMENT"], "checked\nnew note");
    }

    #[test]
    fn exchange_quant_ions_changes_nothing_else() {
        let mut old = LibraryCompound::new("Old");
        old.auto_quant_ions = vec![1.0];
        old.spectrum.rt_ms = 7.0;
        let mut lib = vec![old];
        add(&mut lib, found("x"), AddMode::QuantIons(0), "ignored").unwrap();
        assert_eq!(lib[0].auto_quant_ions, vec![45.0, 73.0, 147.0]);
        assert_eq!(lib[0].spectrum.rt_ms, 7.0);
        assert!(lib[0].spectrum.masses.is_empty());
        assert!(lib[0].features.is_empty());
    }

    #[test]
    fn a_added_compound_is_found_again_with_score_near_one() {
        let c = comp(&[(45, 100), (73, 900), (147, 400), (200, 250)], 2);
        let mut lib = Vec::new();
        let new = compound_from_component(
            "N",
            &c,
            &TIMES,
            Some(1000.0),
            &[],
            &QuantSettings::default(),
        );
        add(&mut lib, new, AddMode::New, "").unwrap();
        let s = SearchSettings::default();
        let id = map_compounds_to_library(
            std::slice::from_ref(&c),
            Some(&[1000.0]),
            &PreparedLibrary::new(&lib, &s),
            &s,
        );
        let hit = id.assigned[0].expect("assigned");
        assert!(hit.overall > 0.99, "{}", hit.overall);
    }

    #[test]
    fn copy_rt_plans_and_applies() {
        let mut lib = vec![found("A"), LibraryCompound::new("Other")];
        lib[0].spectrum.rt_ms = 1.0;
        lib[0].spectrum.ri = 1290.0;
        let c = comp(&[(45, 100), (73, 900), (147, 400)], 3);
        let s = SearchSettings::default();
        let plan = plan_rt_copy(std::slice::from_ref(&c), &TIMES, Some(&[1300.0]), &lib, &s);
        assert_eq!(plan.len(), 1);
        let p = &plan[0];
        assert_eq!(
            (p.target, p.name.as_str(), p.rt_ms, p.ri),
            (0, "A", 4000.0, Some(1300.0))
        );
        assert_eq!((p.old_rt_ms, p.old_ri), (1.0, 1290.0));
        apply_rt_copy(&mut lib, &plan).unwrap();
        assert_eq!(
            (lib[0].spectrum.rt_ms, lib[0].spectrum.ri),
            (4000.0, 1300.0)
        );
        assert!(!lib[0].spectrum.masses.is_empty());
        // Without a calibration the RI stays.
        let plan = plan_rt_copy(&[c], &TIMES, None, &lib, &s);
        lib[0].spectrum.ri = 5.0;
        apply_rt_copy(&mut lib, &plan).unwrap();
        assert_eq!(lib[0].spectrum.ri, 5.0);
        let bad = RtCopy {
            target: 9,
            ..plan[0].clone()
        };
        assert_eq!(apply_rt_copy(&mut lib, &[bad]), Err(AddError::NoEntry(9)));
    }

    #[test]
    fn a_component_keeps_its_ions_when_the_width_changes() {
        use crate::batch::ion_peaks_by;
        use crate::peaks::IonTraces;

        // The components were found with width 1; at width 30 their scan ranges differ, so a search
        // by range in a fresh deconvolution misses them (#312). The ions come from the stored peaks.
        let (run, _) = crate::deconv::tests::two_substances();
        let (ps, saved) = (PeakSettings::default(), DeconvSettings::default());
        let wide = DeconvSettings {
            width: Some(30.0),
            ..saved.clone()
        };
        let mut traces = IonTraces::new(&run, 0.5).unwrap();
        let found = traces.detect_peaks(&ps);
        let comps = deconvolute(&traces, &found, &run.scan_times_s, &saved);
        assert_eq!(comps.len(), 2);
        let now = deconvolute(&traces, &found, &run.scan_times_s, &wide);
        let range = |k: &Component| (k.first_scan, k.apex_scan, k.last_scan);
        assert!(!now.iter().any(|k| range(k) == range(&comps[1])));

        let ions = ion_peaks_by(&comps[1], &found, |i| traces.masses[found[i].row as usize]);
        assert!(ions.iter().all(|p| (62..74).contains(&p.mass)), "{ions:?}");
        assert!(ions.windows(2).all(|w| w[0].integral >= w[1].integral));

        let new = compound_from_component(
            "N",
            &comps[1],
            &run.scan_times_s,
            None,
            &ions,
            &QuantSettings::default(),
        );
        assert_eq!(new.auto_quant_ions.len(), 3);
        let mut lib = Vec::new();
        add(&mut lib, new, AddMode::New, "").unwrap();
        let s = SearchSettings::default();
        let id = map_compounds_to_library(&comps, None, &PreparedLibrary::new(&lib, &s), &s);
        assert!(id.assigned[0].is_none_or(|h| h.overall < 0.5));
        let hit = id.assigned[1].expect("assigned");
        assert!(hit.overall > 0.99, "{}", hit.overall);
    }

    #[test]
    fn check_names_unusable_entries() {
        let mut c = LibraryCompound::new("Test");
        c.quant_ions = vec![123566.0];
        // The reported case: a run-together ion list and no spectrum.
        assert_eq!(check(&c), [Issue::NoSpectrum]);
        c.spectrum.masses = vec![73.0, 174.0, 248.0];
        c.spectrum.values = vec![999, 0, 500];
        assert_eq!(check(&c), [Issue::QuantIonOutOfRange(123566.0)]);
        c.quant_ions = vec![74.0, 174.0, 248.0];
        assert_eq!(
            check(&c),
            [
                Issue::QuantIonNotInSpectrum(74.0),
                Issue::QuantIonNotInSpectrum(174.0)
            ]
        );
        c.quant_ions = vec![f32::NAN];
        assert_eq!(check(&c).len(), 1);
        c.quant_ions = vec![73.0, 248.0];
        assert!(check(&c).is_empty());
    }
}
