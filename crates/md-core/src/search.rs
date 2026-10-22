//! Library search: scores deconvoluted components against a compound library and assigns each library
//! compound to at most one component, the old `LibrarySearch<int,float>`
//! (docs/reverse-engineering/algorithms.md §5). Deviations are listed in docs/architecture.md.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use rayon::prelude::*;

use crate::deconv::Component;
use crate::job::Cancel;
use crate::library::LibraryCompound;
use crate::ri::{RiCalibration, pick_reference};

/// The old `LS_*` settings; defaults are the old static values (file-formats.md §5.1).
#[derive(Debug, Clone, PartialEq)]
pub struct SearchSettings {
    /// `LS_RI_DIFF`: sigma of the Gaussian RI score; differences above 3 sigma score 0.
    pub ri_diff: f64,
    /// `LS_PURE_FACTOR` and `LS_IMPURE_FACTOR`: weights of the pure and impure spectrum score.
    pub pure_factor: f64,
    pub impure_factor: f64,
    /// `LS_SCALE_FACTOR_ION` and `LS_SCALE_FACTOR_VALUE`: each peak is weighted as mass^a * intensity^b.
    pub scale_ion: f64,
    pub scale_value: f64,
    /// `LS_SCALE_FACTOR_NO_MATCH`: factor for peaks without a partner in the other spectrum.
    pub no_match: f64,
    /// `LS_USE_RI` and `LS_USE_SPEC`: which scores make up the overall score.
    pub use_ri: bool,
    pub use_spec: bool,
    /// `LS_MASS_FILTER`: inclusive mass ranges ignored on both sides.
    pub mass_filter: Vec<(f32, f32)>,
    /// `LS_CUTOFF`: minimum overall score of a hit.
    pub cutoff: f64,
    /// `SQL_MAX_NUMBER_FRAGMENTS`: the strongest query fragments the NIST pre-filter compares.
    pub max_fragments: usize,
    /// `SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS`, applied as `shared > n`: with 3 an entry must share at
    /// least 4 of them.
    pub min_equal_fragments: usize,
}

impl Default for SearchSettings {
    fn default() -> Self {
        Self {
            ri_diff: 20.0,
            pure_factor: 0.5,
            impure_factor: 0.5,
            scale_ion: 2.0,
            scale_value: 0.6,
            no_match: 1.0,
            use_ri: true,
            use_spec: true,
            mass_filter: Vec::new(),
            cutoff: 0.8,
            max_fragments: 5,
            min_equal_fragments: 3,
        }
    }
}

impl SearchSettings {
    /// Old `filterMass`: true if `mass` lies in an excluded range.
    fn filtered(&self, mass: u32) -> bool {
        let m = mass as f32;
        self.mass_filter.iter().any(|(lo, hi)| *lo <= m && m <= *hi)
    }

    /// Old `scaleCompound` for one peak.
    fn scale(&self, mass: u32, intensity: f64) -> f64 {
        (mass as f64).powf(self.scale_ion) * intensity.powf(self.scale_value)
    }
}

/// The library with every spectrum binned to nominal mass, mass-filtered and scaled once, like the old
/// pre-scaled cache (`createScaledLib`). Spectra are flat arrays: entry `i` owns `offsets[i]..offsets[i + 1]`.
#[derive(Debug, Clone)]
pub struct PreparedLibrary {
    offsets: Vec<usize>,
    masses: Vec<u32>,
    values: Vec<f64>,
    /// Sum of the squared scaled values per entry.
    norms: Vec<f64>,
    /// RI per entry; `None` where the library has none (old value -1).
    ris: Vec<Option<f64>>,
    rt_only: bool,
    /// Fragment ions per entry (old `red_spectrum.frag`), ascending: entry `i` owns
    /// `frag_offsets[i]..frag_offsets[i + 1]`.
    frag_offsets: Vec<usize>,
    frag_ions: Vec<u32>,
    fragment_filter: bool,
}

/// Old `Spectrum::getSpectrumFragments` and `Fragment` (file-formats.md §6.3) over ascending, unique
/// masses: a fragment is a falling run (each next mass is 1 higher with a lower intensity) holding at
/// least 1 % of the total signal, extended downwards by the lone peaks right below it (M-1, M-2, …).
/// Returns (ion, summed intensity), the ion being the most intense mass (the lowest on ties).
fn fragments(peaks: &[(u32, f64)]) -> Vec<(u32, f64)> {
    let total: f64 = peaks.iter().map(|p| p.1).sum();
    let mut lone = vec![false; peaks.len()];
    let mut runs = Vec::new();
    let mut start = None;
    for (i, &(m, v)) in peaks.iter().enumerate() {
        let falls = peaks.get(i + 1).is_some_and(|&(n, w)| n == m + 1 && v > w);
        match (falls, start) {
            (true, None) => start = Some(i),
            (true, Some(_)) => {}
            (false, Some(s)) => {
                let sum: f64 = peaks[s..=i].iter().map(|p| p.1).sum();
                if 100.0 * sum / total >= 1.0 {
                    runs.push((s, i));
                }
                start = None;
            }
            (false, None) => lone[i] = true,
        }
    }
    runs.into_iter()
        .map(|(mut s, e)| {
            while s > 0 && lone[s - 1] && peaks[s - 1].0 + 1 == peaks[s].0 {
                s -= 1;
            }
            let run = &peaks[s..=e];
            let ion = run
                .iter()
                .fold(run[0], |a, &p| if p.1 > a.1 { p } else { a });
            (ion.0, run.iter().map(|p| p.1).sum())
        })
        .collect()
}

impl PreparedLibrary {
    pub fn new(library: &[LibraryCompound], s: &SearchSettings) -> Self {
        Self::with_ri(library, s, |c| c.spectrum.ri, false)
    }

    /// RT-only alignment (old "Use Retention time for chromatogram alignment", §5.6): every entry's
    /// retention time in s takes the place of its RI. Use with [`RiCalibration::identity`], so the
    /// components' "RI" is their retention time in s too. Not reproducible between instruments.
    pub fn new_rt_only(library: &[LibraryCompound], s: &SearchSettings) -> Self {
        Self::with_ri(library, s, |c| c.spectrum.rt_ms / 1000.0, true)
    }

    fn with_ri(
        library: &[LibraryCompound],
        s: &SearchSettings,
        ri: fn(&LibraryCompound) -> f64,
        rt_only: bool,
    ) -> Self {
        let mut lib = Self {
            offsets: vec![0],
            masses: Vec::new(),
            values: Vec::new(),
            norms: Vec::new(),
            ris: Vec::new(),
            rt_only,
            frag_offsets: vec![0],
            frag_ions: Vec::new(),
            fragment_filter: false,
        };
        for c in library {
            lib.push_with(c, s, ri);
        }
        lib
    }

    /// Appends one entry, scaled like the others (the entry gets the next index).
    pub fn push(&mut self, c: &LibraryCompound, s: &SearchSettings) {
        self.push_with(c, s, |c| c.spectrum.ri);
    }

    fn push_with(
        &mut self,
        c: &LibraryCompound,
        s: &SearchSettings,
        ri: fn(&LibraryCompound) -> f64,
    ) {
        // shortcut: rounds to nominal mass and sums, because the old getBinnedCentroidCompound is not
        // decoded (§5.9) and components are unit-mass; high-resolution search needs the old binning.
        let mut peaks: Vec<(u32, f64)> = c
            .spectrum
            .masses
            .iter()
            .zip(&c.spectrum.values)
            .filter(|(m, v)| m.is_finite() && **m >= 0.5 && **v >= 0)
            .map(|(m, v)| (m.round() as u32, *v as f64))
            .collect();
        peaks.sort_unstable_by_key(|p| p.0);
        let merged: Vec<(u32, f64)> = peaks
            .chunk_by(|a, b| a.0 == b.0)
            .map(|c| (c[0].0, c.iter().map(|p| p.1).sum()))
            .collect();
        let mut norm = 0.0;
        for &(mass, v) in &merged {
            if s.filtered(mass) {
                continue;
            }
            let v = s.scale(mass, v);
            self.masses.push(mass);
            self.values.push(v);
            norm += v * v;
        }
        // The old index is built from the stored spectrum, before mass filter and scaling.
        // Runs are disjoint and in mass order, so the ions come out ascending.
        self.frag_ions
            .extend(fragments(&merged).into_iter().map(|f| f.0));
        self.frag_offsets.push(self.frag_ions.len());
        self.offsets.push(self.masses.len());
        self.norms.push(norm);
        self.ris.push(known_ri(ri(c)));
    }

    /// Built with [`Self::new_rt_only`].
    pub fn is_rt_only(&self) -> bool {
        self.rt_only
    }

    /// Old NIST search (`SQLLibrary::getPotentialHits`): only entries that share more than
    /// `min_equal_fragments` of the `max_fragments` strongest fragments of the component are scored
    /// ([`SearchSettings`]). Use for [`crate::library::is_nist`] libraries.
    pub fn with_fragment_filter(mut self, on: bool) -> Self {
        self.fragment_filter = on;
        self
    }

    pub fn len(&self) -> usize {
        self.norms.len()
    }

    pub fn is_empty(&self) -> bool {
        self.norms.is_empty()
    }
}

/// The old code treats a negative RI as missing.
fn known_ri(ri: f64) -> Option<f64> {
    (ri >= 0.0).then_some(ri)
}

/// One library compound scored against one component.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LibraryHit {
    /// Index into the library given to `PreparedLibrary::new`.
    pub compound: usize,
    /// `PF * pure + IF * impure`.
    pub spec_sim: f64,
    /// Gaussian RI score; `None` when the component or the library compound has no RI.
    pub ri_sim: Option<f64>,
    pub overall: f64,
}

/// Old `calculateRISimScore`.
pub fn ri_score(ri1: f64, ri2: f64, sigma: f64) -> f64 {
    let d = ri1 - ri2;
    if d.abs() > 3.0 * sigma {
        0.0
    } else {
        (-d * d / (2.0 * sigma * sigma)).exp()
    }
}

/// Old `getLibraryHits`: scores `spectrum` (nominal mass, height; ascending, unique masses like a
/// `Component`) against every library entry. Returns the hits with an overall score above zero and at
/// least `s.cutoff`, best first; ties keep library order.
pub fn library_hits(
    spectrum: &[(u32, u32)],
    ri: Option<f64>,
    lib: &PreparedLibrary,
    s: &SearchSettings,
) -> Vec<LibraryHit> {
    // Dense lookup by mass; -1 marks masses the component lacks (scaled values are >= 0).
    let top = spectrum.iter().map(|p| p.0 as usize + 1).max().unwrap_or(0);
    let mut obs = vec![-1.0; top];
    let mut obs_norm = 0.0;
    for (m, h) in spectrum {
        if !s.filtered(*m) {
            let v = s.scale(*m, *h as f64);
            obs[*m as usize] = v;
            obs_norm += v * v;
        }
    }
    let ri = ri.and_then(known_ri);
    let nm2 = s.no_match * s.no_match;
    let query: Option<Vec<u32>> = lib.fragment_filter.then(|| {
        let peaks: Vec<_> = spectrum.iter().map(|&(m, h)| (m, h as f64)).collect();
        let mut f = fragments(&peaks);
        // Old `compare_intens` order: strongest first.
        f.sort_by(|a, b| b.1.total_cmp(&a.1));
        f.iter().take(s.max_fragments).map(|f| f.0).collect()
    });
    let mut hits: Vec<LibraryHit> = (0..lib.len())
        .filter_map(|k| {
            if let Some(q) = &query {
                let ions = &lib.frag_ions[lib.frag_offsets[k]..lib.frag_offsets[k + 1]];
                let shared = q.iter().filter(|i| ions.binary_search(i).is_ok()).count();
                if shared <= s.min_equal_fragments {
                    return None;
                }
            }
            let ri_sim = ri.zip(lib.ris[k]).map(|(a, b)| ri_score(b, a, s.ri_diff));
            if s.use_ri && ri_sim.is_some_and(|r| r <= 0.0) {
                return None;
            }
            // With unique integer masses the old merge (§5.3/5.4) only needs the matched sums.
            let (mut num, mut lib_matched, mut obs_matched) = (0.0, 0.0, 0.0);
            let r = lib.offsets[k]..lib.offsets[k + 1];
            for (m, vl) in lib.masses[r.clone()].iter().zip(&lib.values[r]) {
                if let Some(vo) = obs.get(*m as usize).filter(|v| **v >= 0.0) {
                    num += vl * vo;
                    lib_matched += vl * vl;
                    obs_matched += vo * vo;
                }
            }
            let d1 = nm2 * (lib.norms[k] - lib_matched) + lib_matched;
            let cosine = |d2: f64| {
                if d1 > 0.0 && d2 > 0.0 {
                    num / (d1.sqrt() * d2.sqrt())
                } else {
                    0.0
                }
            };
            let pure = cosine(nm2 * (obs_norm - obs_matched) + obs_matched);
            let impure = cosine(obs_matched);
            let spec_sim = s.pure_factor * pure + s.impure_factor * impure;
            // An entry without RI is scored by its spectrum only.
            let ri_used = s.use_ri.then_some(ri_sim).flatten();
            let overall = match (s.use_spec, ri_used) {
                (true, Some(r)) => (spec_sim * spec_sim * r).cbrt(),
                (true, None) => spec_sim.abs(),
                (false, Some(r)) => r,
                // Only RI asked for but none known: nothing to score (old: RT fallback, not ported).
                (false, None) if s.use_ri => 0.0,
                (false, None) => 1.0,
            };
            (overall > 0.0 && overall >= s.cutoff).then_some(LibraryHit {
                compound: k,
                spec_sim,
                ri_sim,
                overall,
            })
        })
        .collect();
    hits.sort_by(|a, b| b.overall.total_cmp(&a.overall));
    hits
}

/// Result of `map_compounds_to_library`, indexed like the components.
#[derive(Debug, Clone, PartialEq)]
pub struct Identification {
    /// All hits of each component (overall >= cutoff), best first.
    pub hits: Vec<Vec<LibraryHit>>,
    /// The library compound assigned to each component, one-to-one.
    pub assigned: Vec<Option<LibraryHit>>,
}

/// Old `mapCompoundsToLibrary`: searches every component in parallel, then assigns greedily by best
/// overall score so that each library compound and each component is used at most once.
/// `ris[i]` is the RI of component `i`; pass `None` while the run is not RI-calibrated.
pub fn map_compounds_to_library(
    comps: &[Component],
    ris: Option<&[f64]>,
    lib: &PreparedLibrary,
    s: &SearchSettings,
) -> Identification {
    map_compounds_cancellable(comps, ris, lib, s, &Cancel::default())
        .expect("a fresh Cancel is never set")
}

/// [`map_compounds_to_library`] that checks `cancel` before each component's search and returns `None`
/// once it is set; otherwise the result is identical.
pub fn map_compounds_cancellable(
    comps: &[Component],
    ris: Option<&[f64]>,
    lib: &PreparedLibrary,
    s: &SearchSettings,
    cancel: &Cancel,
) -> Option<Identification> {
    if cancel.is_cancelled() {
        return None; // also for no components
    }
    // Not `collect::<Option<_>>`: that is measurably slower than a plain collect, as with `sic`.
    let hits: Vec<Vec<LibraryHit>> = comps
        .par_iter()
        .enumerate()
        .map(|(i, c)| match cancel.is_cancelled() {
            true => Vec::new(),
            false => library_hits(&c.spectrum, ris.and_then(|r| r.get(i).copied()), lib, s),
        })
        .collect();
    if cancel.is_cancelled() {
        return None;
    }
    let assigned = assign_best_first(&hits, lib.len());
    Some(Identification { hits, assigned })
}

/// The assignment step of `map_compounds_to_library` on precomputed hits (`hits[c]` best first):
/// one library compound of `lib_len` per component, best overall scores first. Each list must be sorted by
/// overall descending, then compound ascending, as `library_hits` returns it.
pub fn assign_best_first(hits: &[Vec<LibraryHit>], lib_len: usize) -> Vec<Option<LibraryHit>> {
    // Taking pairs in descending score equals the old loop, which picks the best free pair each round.
    // Ties go to the lower library index, then the lower component index (old: pointer order).
    // Every list is already in that order, so a heap over each component's current best pair replaces
    // sorting all pairs: a component whose compound is taken moves on to its next hit.
    let mut heap: BinaryHeap<Head> = (hits.iter().enumerate())
        .filter_map(|(c, h)| {
            Some(Head {
                hit: *h.first()?,
                comp: c,
                next: 1,
            })
        })
        .collect();
    let mut assigned = vec![None; hits.len()];
    let mut used = vec![false; lib_len];
    while let Some(Head { hit, comp, next }) = heap.pop() {
        if used[hit.compound] {
            if let Some(&hit) = hits[comp].get(next) {
                heap.push(Head {
                    hit,
                    comp,
                    next: next + 1,
                });
            }
        } else {
            assigned[comp] = Some(hit);
            used[hit.compound] = true;
        }
    }
    assigned
}

/// A component's best unassigned pair; the heap pops the greatest, i.e. the best score first.
struct Head {
    hit: LibraryHit,
    comp: usize,
    next: usize,
}

impl Ord for Head {
    fn cmp(&self, o: &Self) -> Ordering {
        (self.hit.overall.total_cmp(&o.hit.overall))
            .then(o.hit.compound.cmp(&self.hit.compound))
            .then(o.comp.cmp(&self.comp))
    }
}
impl PartialOrd for Head {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl PartialEq for Head {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Head {}

/// RT shift correction (§6.5) for one run: finds library compound `reference` among `comps` within
/// `3 * ri_diff` of its library RI, with an overall score of at least `cutoff`, and returns `cal` shifted
/// so that this component's apex gets the library RI. `None` when the reference has no RI or is not found.
pub fn shift_to_reference(
    cal: &RiCalibration,
    scan_times_s: &[f64],
    comps: &[Component],
    lib: &PreparedLibrary,
    reference: usize,
    cutoff: f64,
    s: &SearchSettings,
) -> Option<RiCalibration> {
    let reference_ri = lib.ris[reference]?;
    let ris = cal.component_ris(scan_times_s, comps);
    let score = |i: usize| {
        library_hits(&comps[i].spectrum, Some(ris[i]), lib, s)
            .iter()
            .find(|h| h.compound == reference)
            .map_or(0.0, |h| h.overall)
    };
    let i = pick_reference(&ris, reference_ri, 3.0 * s.ri_diff, cutoff, score)?;
    Some(cal.shifted_to(scan_times_s[comps[i].apex_scan as usize], reference_ri))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::library::LibrarySpectrum;

    fn entry(name: &str, ri: f64, peaks: &[(f32, i32)]) -> LibraryCompound {
        LibraryCompound {
            spectrum: LibrarySpectrum {
                scan: 0,
                rt_ms: 0.0,
                ri,
                bin_width: 1.0,
                masses: peaks.iter().map(|p| p.0).collect(),
                values: peaks.iter().map(|p| p.1).collect(),
            },
            name: name.into(),
            fwhm: None,
            integrals: None,
            quant_ions: Vec::new(),
            auto_quant_ions: Vec::new(),
            features: Default::default(),
            drt: 0.0,
            time_stamp: 0,
            id: 0,
            comb_spec: Vec::new(),
            ion_freq: Vec::new(),
            spec_count: 0,
            ris: Vec::new(),
            blocks: Vec::new(),
            sources: Vec::new(),
            peaks: Vec::new(),
            model: None,
            formulas: Vec::new(),
            mid_sizes: Vec::new(),
        }
    }

    fn component(spectrum: &[(u32, u32)]) -> Component {
        Component {
            first_scan: 0,
            apex_scan: 0,
            last_scan: 0,
            model: Vec::new(),
            spectrum: spectrum.to_vec(),
            peaks: Vec::new(),
            exact_spectrum: Vec::new(),
        }
    }

    const ALA: [(f32, i32); 4] = [(73.0, 1000), (116.0, 800), (147.0, 300), (190.0, 50)];
    const GLY: [(f32, i32); 4] = [(73.0, 900), (102.0, 1000), (147.0, 200), (176.0, 80)];

    fn as_comp(peaks: &[(f32, i32)], factor: u32) -> Vec<(u32, u32)> {
        peaks
            .iter()
            .map(|(m, v)| (*m as u32, *v as u32 * factor))
            .collect()
    }

    fn all(s: SearchSettings) -> SearchSettings {
        SearchSettings { cutoff: 0.0, ..s }
    }

    #[test]
    fn same_spectrum_scores_one_at_any_scale() {
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(&[entry("Ala", -1.0, &ALA)], &s);
        let h = library_hits(&as_comp(&ALA, 3), None, &lib, &s);
        assert_eq!(h.len(), 1);
        assert!((h[0].spec_sim - 1.0).abs() < 1e-12 && (h[0].overall - 1.0).abs() < 1e-12);
        assert_eq!(h[0].ri_sim, None);
    }

    #[test]
    fn foreign_spectrum_scores_zero_and_is_no_hit() {
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(&[entry("X", -1.0, &[(50.0, 10), (60.0, 10)])], &s);
        assert!(library_hits(&as_comp(&ALA, 1), None, &lib, &s).is_empty());
    }

    #[test]
    fn extra_and_shifted_peaks_lower_pure_more_than_impure() {
        // Component = Ala plus a co-eluting ion at 300, and 190 shifted by one mass.
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(&[entry("Ala", -1.0, &ALA)], &s);
        let comp = [(73, 1000), (116, 800), (147, 300), (191, 50), (300, 400)];
        let h = library_hits(&comp, None, &lib, &s)[0];
        let only = |pf, impf| {
            let t = SearchSettings {
                pure_factor: pf,
                impure_factor: impf,
                ..s.clone()
            };
            library_hits(
                &comp,
                None,
                &PreparedLibrary::new(&[entry("Ala", -1.0, &ALA)], &t),
                &t,
            )[0]
            .spec_sim
        };
        let (pure, impure) = (only(1.0, 0.0), only(0.0, 1.0));
        assert!(0.0 < pure && pure < impure && impure < 1.0);
        assert!((h.spec_sim - (pure + impure) / 2.0).abs() < 1e-12);
        // Hand-computed impure score: lib 190 has no partner, obs 191 and 300 are ignored.
        let w = |m: f64, v: f64| m.powi(2) * v.powf(0.6);
        let matched: Vec<f64> = [(73.0, 1000.0), (116.0, 800.0), (147.0, 300.0)]
            .iter()
            .map(|(m, v)| w(*m, *v))
            .collect();
        let sq: f64 = matched.iter().map(|v| v * v).sum();
        let expected = sq / ((sq + w(190.0, 50.0).powi(2)).sqrt() * sq.sqrt());
        assert!((impure - expected).abs() < 1e-12);
    }

    #[test]
    fn mass_filter_ignores_ranges_on_both_sides() {
        let s = all(SearchSettings {
            mass_filter: vec![(146.5, 147.5), (299.0, 301.0)],
            ..Default::default()
        });
        let lib = PreparedLibrary::new(&[entry("Ala", -1.0, &ALA)], &s);
        let comp = [(73, 1000), (116, 800), (147, 1), (190, 50), (300, 9999)];
        assert!((library_hits(&comp, None, &lib, &s)[0].overall - 1.0).abs() < 1e-12);
    }

    #[test]
    fn ri_score_combines_and_gates() {
        assert_eq!(ri_score(1000.0, 1061.0, 20.0), 0.0);
        assert!((ri_score(1000.0, 1020.0, 20.0) - (-0.5f64).exp()).abs() < 1e-15);
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(
            &[
                entry("near", 1020.0, &ALA),
                entry("far", 1100.0, &ALA),
                entry("no RI", -1.0, &ALA),
            ],
            &s,
        );
        let h = library_hits(&as_comp(&ALA, 1), Some(1000.0), &lib, &s);
        // "far" is beyond 3 sigma and dropped; "no RI" counts by spectrum only and ranks first.
        assert_eq!(h.iter().map(|h| h.compound).collect::<Vec<_>>(), [2, 0]);
        assert!((h[1].overall - (-0.5f64).exp().cbrt()).abs() < 1e-12);
        // Without RI use the RI only shows up as information.
        let s = SearchSettings { use_ri: false, ..s };
        let h = library_hits(&as_comp(&ALA, 1), Some(1000.0), &lib, &s);
        assert_eq!(h.len(), 3);
        assert!(h.iter().all(|h| (h.overall - 1.0).abs() < 1e-12));
        assert_eq!(h[1].ri_sim, Some(0.0));
        // RI only: entries without RI are no hits.
        let s = SearchSettings {
            use_ri: true,
            use_spec: false,
            ..s
        };
        let h = library_hits(&as_comp(&ALA, 1), Some(1000.0), &lib, &s);
        assert_eq!(h.iter().map(|h| h.compound).collect::<Vec<_>>(), [0]);
        assert!(library_hits(&as_comp(&ALA, 1), None, &lib, &s).is_empty());
    }

    #[test]
    fn library_peaks_are_binned_to_nominal_mass() {
        let s = all(SearchSettings::default());
        let split = [
            (72.9, 600),
            (73.2, 400),
            (116.0, 800),
            (147.0, 300),
            (190.0, 50),
            (0.2, 5),
            (200.0, -3),
            (300.0, 0),
        ];
        let lib = PreparedLibrary::new(&[entry("Ala", -1.0, &split)], &s);
        assert!((library_hits(&as_comp(&ALA, 1), None, &lib, &s)[0].overall - 1.0).abs() < 1e-12);
    }

    /// The former assignment: sort all (component, hit) pairs, then take the free ones in order.
    fn assign_sorted_pairs(hits: &[Vec<LibraryHit>], lib_len: usize) -> Vec<Option<LibraryHit>> {
        let mut pairs: Vec<(usize, usize)> = (hits.iter().enumerate())
            .flat_map(|(c, h)| (0..h.len()).map(move |j| (c, j)))
            .collect();
        let hit = |p: &(usize, usize)| hits[p.0][p.1];
        pairs.sort_unstable_by(|a, b| {
            let (x, y) = (hit(a), hit(b));
            (y.overall.total_cmp(&x.overall))
                .then(x.compound.cmp(&y.compound))
                .then(a.0.cmp(&b.0))
        });
        let (mut assigned, mut used) = (vec![None; hits.len()], vec![false; lib_len]);
        for p in &pairs {
            let h = hit(p);
            if assigned[p.0].is_none() && !used[h.compound] {
                assigned[p.0] = Some(h);
                used[h.compound] = true;
            }
        }
        assigned
    }

    #[test]
    fn heap_assignment_equals_sorting_all_pairs_with_many_ties() {
        // Scores from a few levels so that ties between compounds and components are everywhere.
        let mut x = 12345u64;
        let mut rnd = move |n: u64| {
            x = x
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (x >> 33) % n
        };
        for _ in 0..50 {
            let (comps, lib_len) = (1 + rnd(30) as usize, 1 + rnd(30) as usize);
            let hits: Vec<Vec<LibraryHit>> = (0..comps)
                .map(|_| {
                    let mut h = Vec::new();
                    for k in 0..lib_len {
                        if rnd(3) > 0 {
                            let overall = (1 + rnd(4)) as f64 / 4.0;
                            h.push(LibraryHit {
                                compound: k,
                                spec_sim: 1.0,
                                ri_sim: None,
                                overall,
                            });
                        }
                    }
                    h.sort_by(|a, b| b.overall.total_cmp(&a.overall));
                    h
                })
                .collect();
            assert_eq!(
                assign_best_first(&hits, lib_len),
                assign_sorted_pairs(&hits, lib_len)
            );
        }
    }

    #[test]
    fn mapping_is_one_to_one_by_best_score() {
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(&[entry("Ala", -1.0, &ALA), entry("Gly", -1.0, &GLY)], &s);
        // Two components look most like Ala; the better one gets it, the other falls back to Gly.
        let near_ala = [(73, 1000), (102, 300), (116, 800), (147, 300)];
        let comps = [
            component(&near_ala),
            component(&as_comp(&ALA, 2)),
            component(&[(50, 1)]),
        ];
        let id = map_compounds_to_library(&comps, None, &lib, &s);
        assert_eq!(id.hits[0][0].compound, 0);
        assert!(id.hits[2].is_empty());
        let got: Vec<_> = id.assigned.iter().map(|a| a.map(|h| h.compound)).collect();
        assert_eq!(got, [Some(1), Some(0), None]);
        // A cutoff above the fallback score leaves the first component unassigned.
        let cut = SearchSettings { cutoff: 0.99, ..s };
        let id = map_compounds_to_library(&comps, None, &lib, &cut);
        assert_eq!(
            id.assigned
                .iter()
                .map(|a| a.map(|h| h.compound))
                .collect::<Vec<_>>(),
            [None, Some(0), None]
        );
    }

    #[test]
    fn set_cancel_ends_the_search_without_result() {
        let lib = vec![entry("A", 100.0, &[(50.0, 999)])];
        let s = all(SearchSettings::default());
        let comps = [component(&[(50, 100)])];
        let prepared = PreparedLibrary::new(&lib, &s);
        let cancel = Cancel::default();
        assert!(map_compounds_cancellable(&comps, None, &prepared, &s, &cancel).is_some());
        cancel.cancel();
        assert!(map_compounds_cancellable(&comps, None, &prepared, &s, &cancel).is_none());
        assert!(map_compounds_cancellable(&[], None, &prepared, &s, &cancel).is_none());
        assert!(
            map_compounds_cancellable(&comps, None, &prepared, &s, &Cancel::after(2)).is_none()
        );
    }

    #[test]
    fn shift_moves_reference_onto_its_library_ri() {
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new(
            &[entry("Ala", 1100.0, &ALA), entry("Gly", 1200.0, &GLY)],
            &s,
        );
        let cal = RiCalibration::new([(100.0, 1000.0), (300.0, 1300.0)]).unwrap();
        // Scan i at i seconds; the run is 10 s late, Gly elutes at 243 s instead of 233.3 s.
        let times: Vec<f64> = (0..400).map(f64::from).collect();
        let at = |scan, peaks: &[(f32, i32)]| Component {
            apex_scan: scan,
            ..component(&as_comp(peaks, 1))
        };
        let comps = [at(187, &ALA), at(243, &GLY)];
        let shifted = shift_to_reference(&cal, &times, &comps, &lib, 1, 0.8, &s).unwrap();
        assert!((shifted.ri(243.0) - 1200.0).abs() < 1e-9);
        assert!((shifted.shift_s - (243.0 - cal.rt(1200.0))).abs() < 1e-9);
        // Another spectrum at that place, a cutoff above the score or a reference outside the window: no shift.
        assert_eq!(
            shift_to_reference(&cal, &times, &comps[..1], &lib, 1, 0.8, &s),
            None
        );
        let far = [at(330, &GLY)];
        assert_eq!(
            shift_to_reference(&cal, &times, &far, &lib, 1, 0.8, &s),
            None
        );
        let weak = [at(243, &ALA)];
        assert_eq!(
            shift_to_reference(&cal, &times, &weak, &lib, 1, 0.8, &s),
            None
        );
        // Gly at 243 s scores (1 · 0.755)^(1/3) = 0.91.
        assert_eq!(
            shift_to_reference(&cal, &times, &comps, &lib, 1, 0.95, &s),
            None
        );
        // Searched with the shifted RIs, Gly matches its library RI exactly.
        let ris = shifted.component_ris(&times, &comps);
        let id = map_compounds_to_library(&comps, Some(&ris), &lib, &s);
        assert_eq!(id.assigned[1].unwrap().ri_sim, Some(1.0));
    }

    /// The old search takes `RI >= 0` as known (`getLibraryHits`): RI 0 is a value, -1 means unknown.
    #[test]
    fn ri_zero_is_known_and_minus_one_is_not() {
        let peaks = [(57.0, 100), (71.0, 50)];
        let s = all(SearchSettings::default());
        let spectrum = as_comp(&peaks, 1);
        let ri_sim = |ri| {
            let lib = PreparedLibrary::new(&[entry("a", ri, &peaks)], &s);
            library_hits(&spectrum, Some(0.0), &lib, &s)[0].ri_sim
        };
        assert_eq!(ri_sim(0.0), Some(1.0));
        assert_eq!(ri_sim(-1.0), None);
    }

    #[test]
    fn rt_only_compares_times_in_seconds() {
        let peaks = [(57.0, 100), (71.0, 50)];
        let mut e = entry("a", 1500.0, &peaks);
        e.spectrum.rt_ms = 20_000.0;
        let s = all(SearchSettings::default());
        let lib = PreparedLibrary::new_rt_only(&[e], &s);
        let spectrum = as_comp(&peaks, 1);
        // The library RI (1500) is ignored: the component at 20 s scores 1, one at 200 s is dropped.
        assert_eq!(
            library_hits(&spectrum, Some(20.0), &lib, &s)[0].ri_sim,
            Some(1.0)
        );
        assert!(library_hits(&spectrum, Some(200.0), &lib, &s).is_empty());
    }

    #[test]
    fn fragments_match_the_old_red_spectrum() {
        let f = |p: &[(u32, f64)]| fragments(p);
        // Rows the old importNIST wrote for the synthetic records (tests/data/synthetic-nist.sdf).
        let kappa = [
            (70, 5.0),
            (71, 100.0),
            (72, 12.0),
            (77, 40.0),
            (78, 9.0),
            (91, 80.0),
        ];
        assert_eq!(f(&kappa), [(71, 117.0), (77, 49.0)]);
        let lambda = [
            (119, 10.0),
            (120, 100.0),
            (133, 100.0),
            (134, 10.0),
            (135, 50.0),
            (150, 1.0),
            (151, 1.0),
            (200, 1000.0),
        ];
        // Rising pairs and lone peaks are no fragment; 135 rises again and stays alone.
        assert_eq!(f(&lambda), [(133, 110.0)]);
        // A falling run below 1 % of the total signal is dropped.
        assert!(f(&[(100, 1.0), (101, 0.5), (200, 1000.0)]).is_empty());
    }
}
