//! The MID wizard above the solver (old `MIDWizard`, `MIDCompoundSeries`, algorithms.md §7g):
//! fragments come from the library, intensities from a batch scan, and one of four corrections
//! turns them into MIDs.

use super::{MidFit, correct_by_formula, correct_by_reference, replicate_mean};

use rayon::prelude::*;

use crate::batch::{Found, RunResult};
use crate::isotope::{IsotopeTable, parse_formula};
use crate::library::LibraryCompound;

/// Correction of the old wizard's method page (`MethodPage::getSelectedMethod`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Method {
    /// 0: mean and confidence interval of the labelled runs' normalised spectra.
    Uncorrected,
    /// 1: the labelled runs against the unlabelled ones. `m1` is the M1 correction, the tracer's
    /// natural abundance / 100, or 0 when it is off.
    Reference { m1: f64 },
    /// 2: every run against the library spectrum, corrected by the tracer's M+1 / M+0.
    Library,
    /// 3: every run against the fragment's sum formula with the tracer replaced by the label.
    Formula,
}

impl Method {
    /// Whether every run gets its own MID (old `Calculate_Individual`) instead of one per group.
    pub fn per_run(self) -> bool {
        matches!(self, Method::Library | Method::Formula)
    }
}

/// A fragment of a library substance: ions `lo..=hi`, its sum formula and MID size (0 = one share
/// per ion).
#[derive(Debug, Clone, PartialEq)]
pub struct Fragment<'a> {
    pub lo: u32,
    pub hi: u32,
    pub formula: &'a str,
    pub size: usize,
}

/// Fragments of a substance in the old order: by the (start, end) pair as entered. Fragments wider
/// than [`MAX_FRAGMENT_IONS`](super::MAX_FRAGMENT_IONS) are left out.
pub fn fragments(c: &LibraryCompound) -> Vec<Fragment<'_>> {
    let mut keys: Vec<&(u32, u32, String)> = (c.formulas.iter())
        .filter(|f| f.0.abs_diff(f.1) < super::MAX_FRAGMENT_IONS)
        .collect();
    keys.sort_by_key(|f| (f.0, f.1));
    keys.into_iter()
        .map(|(a, b, formula)| Fragment {
            lo: *a.min(b),
            hi: *a.max(b),
            formula,
            size: (c.mid_sizes.iter())
                .find(|m| (m.0, m.1) == (*a, *b))
                .map_or(0, |m| m.2 as usize),
        })
        .collect()
}

/// The library the wizard scans with (old `createMIDLibrary`): the substances with fragments,
/// with every fragment ion as quant ion, so the batch rescan fills ions that were not deconvoluted.
pub fn mid_library(lib: &[LibraryCompound]) -> Vec<LibraryCompound> {
    lib.iter()
        .filter(|c| !fragments(c).is_empty())
        .map(|c| {
            let mut c = c.clone();
            let mut ions: Vec<u32> = fragments(&c).iter().flat_map(|f| f.lo..=f.hi).collect();
            ions.sort_unstable();
            ions.dedup();
            c.quant_ions = ions.into_iter().map(|m| m as f32).collect();
            c
        })
        .collect()
}

/// Integrals of ions `lo..=hi` in one run (old `CompoundSeries::getIntegral`), 0 for an ion without peak.
pub fn intensities(found: &Found, lo: u32, hi: u32) -> Vec<f64> {
    (lo..=hi)
        .map(|m| {
            (found.ions.iter())
                .find(|p| p.mass == m)
                .map_or(0.0, |p| p.integral)
        })
        .collect()
}

/// MIDs of one fragment.
#[derive(Debug, Clone, PartialEq)]
pub struct FragmentMid {
    pub lo: u32,
    pub hi: u32,
    /// Rows of the per-run exports (old `MIDCompoundSeries`): one per ion, at most the MID size.
    pub isotopomers: usize,
    /// Per-run methods: one fit per run, `None` where the run lacks the substance or every ion is 0
    /// (the old program shows NaN there). Group methods: the single fit, `None` without usable runs.
    pub fits: Vec<Option<MidFit>>,
}

/// MIDs of one substance of the MID library.
#[derive(Debug, Clone, PartialEq)]
pub struct CompoundMid {
    /// Index into the library passed to `determine`.
    pub compound: usize,
    pub fragments: Vec<FragmentMid>,
}

/// Result of `determine`: substances in library order, and the substances left out with the reason.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Mids {
    pub compounds: Vec<CompoundMid>,
    pub errors: Vec<String>,
}

/// MIDs of every substance found in at least one run. `labelled` and `unlabelled` are indices into
/// `runs` and only matter for the group methods (`Uncorrected`, `Reference`); the per-run methods use
/// every run. A substance whose sum formula or tracer cannot be used is left out with an error, as
/// the old wizard drops it with a warning.
pub fn determine(
    runs: &[RunResult],
    lib: &[LibraryCompound],
    labelled: &[usize],
    unlabelled: &[usize],
    method: Method,
    tracer: &str,
    table: &IsotopeTable,
) -> Mids {
    // Substances are independent: one task each, results back in library order.
    let each: Vec<Result<Option<CompoundMid>, String>> = (lib.par_iter().enumerate())
        .zip(hits(runs, lib.len()))
        .map(|((i, c), found)| {
            if found.iter().all(Option::is_none) {
                return Ok(None);
            }
            let fragments = (fragments(c).iter())
                .map(|f| {
                    fragment_mid(c, &found, f, labelled, unlabelled, method, tracer, table)
                        .map(|fits| FragmentMid {
                            lo: f.lo,
                            hi: f.hi,
                            isotopomers: match (method.per_run(), f.size) {
                                (true, 1..) => f.size.min((f.hi - f.lo) as usize + 1),
                                _ => (f.hi - f.lo) as usize + 1,
                            },
                            fits,
                        })
                        .map_err(|e| format!("{} ({}-{}): {e}", c.name, f.lo, f.hi))
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok((!fragments.is_empty()).then_some(CompoundMid {
                compound: i,
                fragments,
            }))
        })
        .collect();
    let mut out = Mids::default();
    for r in each {
        match r {
            Ok(c) => out.compounds.extend(c),
            Err(e) => out.errors.push(e),
        }
    }
    out
}

/// Per substance `0..n`, its first hit in each run: one pass over the runs instead of a search per
/// substance.
pub(crate) fn hits(runs: &[RunResult], n: usize) -> Vec<Vec<Option<&Found>>> {
    let mut out = vec![vec![None; runs.len()]; n];
    for (r, run) in runs.iter().enumerate() {
        for f in run.found.iter().rev().filter(|f| f.compound < n) {
            out[f.compound][r] = Some(f);
        }
    }
    out
}

/// Rows of the runs `idx` that have the substance and some intensity (old `Calculate` skips the
/// all-zero ones).
fn group_rows(found: &[Option<&Found>], idx: &[usize], f: &Fragment) -> Vec<Vec<f64>> {
    (idx.iter())
        .filter_map(|&r| found.get(r).copied().flatten())
        .map(|x| intensities(x, f.lo, f.hi))
        .filter(|row| row.iter().any(|v| *v != 0.0))
        .collect()
}

#[allow(clippy::too_many_arguments)]
fn fragment_mid(
    c: &LibraryCompound,
    found: &[Option<&Found>],
    f: &Fragment,
    labelled: &[usize],
    unlabelled: &[usize],
    method: Method,
    tracer: &str,
    table: &IsotopeTable,
) -> Result<Vec<Option<MidFit>>, String> {
    let rows = || group_rows(found, labelled, f);
    match method {
        Method::Uncorrected => Ok(vec![replicate_mean(&rows())]),
        Method::Reference { m1 } => Ok(vec![correct_by_reference(
            &rows(),
            &group_rows(found, unlabelled, f),
            m1,
        )]),
        Method::Library => {
            // The tracer's M+1 / M+0 corrects the shifted reference (old `determineMIDs`).
            let tracer = parse_formula(tracer).and_then(|t| match t.is_empty() {
                true => Err("no tracer".into()),
                false => table.pattern(&t, 2),
            })?;
            let correction = tracer.abundance[1] / tracer.abundance[0];
            let reference: Vec<f64> = (f.lo..=f.hi).map(|m| library_intensity(c, m)).collect();
            per_run(found, f, |row| {
                Ok(correct_by_reference(
                    &[row.to_vec()],
                    std::slice::from_ref(&reference),
                    correction,
                ))
            })
        }
        Method::Formula => {
            if f.formula.is_empty() {
                return Err("no sum formula".into());
            }
            let formula = parse_formula(f.formula)?;
            let tracer = Some(tracer).filter(|t| !t.is_empty());
            per_run(found, f, |row| {
                correct_by_formula(row, &formula, table, tracer, f.size)
            })
        }
    }
}

/// Library intensity of ion `m` (old `Compound::getNormalizedIntensity`), 0 if the spectrum lacks it.
fn library_intensity(c: &LibraryCompound, m: u32) -> f64 {
    let s = &c.spectrum;
    (s.masses.iter())
        .position(|x| *x == m as f32)
        .map_or(0.0, |i| s.values[i] as f64)
}

fn per_run(
    found: &[Option<&Found>],
    f: &Fragment,
    solve: impl Fn(&[f64]) -> Result<Option<MidFit>, String>,
) -> Result<Vec<Option<MidFit>>, String> {
    (found.iter())
        .map(|x| match x.map(|x| intensities(x, f.lo, f.hi)) {
            Some(row) if row.iter().any(|v| *v != 0.0) => solve(&row),
            _ => Ok(None),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::IonPeak;

    fn run(ions: &[(u32, f64)]) -> RunResult {
        RunResult {
            path: "r".into(),
            found: vec![Found {
                compound: 0,
                score: 1.0,
                rt_s: 0.0,
                ions: (ions.iter())
                    .map(|&(mass, integral)| IonPeak {
                        mass,
                        integral,
                        quality: 0.0,
                        signal_to_noise: 0.0,
                    })
                    .collect(),
                rescanned: false,
                rescued: Vec::new(),
                alternatives: Vec::new(),
            }],
            unknown: Vec::new(),
            calibration: None,
        }
    }

    #[test]
    fn uncorrected_group_skips_empty_runs_and_exports() {
        let mut c = LibraryCompound::new("Ala");
        c.formulas.push((11, 10, "C2".into()));
        let runs = [
            run(&[(10, 3.0), (11, 1.0)]),
            run(&[(10, 1.0), (11, 1.0)]),
            run(&[(10, 0.0)]),
        ];
        let table = IsotopeTable::standard();
        let m = determine(
            &runs,
            &[c.clone()],
            &[0, 1, 2],
            &[],
            Method::Uncorrected,
            "C",
            &table,
        );
        assert!(m.errors.is_empty());
        let f = &m.compounds[0].fragments[0];
        assert_eq!((f.lo, f.hi), (10, 11));
        let fit = f.fits[0].as_ref().unwrap();
        // Mean of (0.75, 0.25) and (0.5, 0.5); the all-zero third run does not count.
        assert_eq!(fit.mids, [0.625, 0.375]);
        // 1.96 * sd / sqrt(n) with sd = 0.125 * sqrt(2), n = 2.
        let ci = 1.96 * 0.125;
        assert!((fit.cis[0].unwrap() - ci).abs() < 1e-12);
        assert_eq!(
            m.group_text(&[c]),
            format!(
                "Ala_10 (M0)\t0.625\t{}\t-1\nAla_11 (M1)\t0.375\t{}\t-1\n\n",
                crate::mid::qt_g(ci),
                crate::mid::qt_g(fit.cis[1].unwrap())
            )
        );
    }

    #[test]
    fn reference_group_explains_a_half_labelled_run() {
        let mut c = LibraryCompound::new("Ala");
        c.formulas.push((10, 12, "C2".into()));
        // Reference 0.9 / 0.1; half of the labelled run carries one label: 0.45 / 0.5 / 0.05.
        let runs = [
            run(&[(10, 0.45), (11, 0.5), (12, 0.05)]),
            run(&[(10, 9.0), (11, 1.0)]),
            run(&[(10, 0.0)]),
        ];
        let table = IsotopeTable::standard();
        let m = determine(
            &runs,
            &[c],
            &[0],
            &[1, 2],
            Method::Reference { m1: 0.0 },
            "C",
            &table,
        );
        assert!(m.errors.is_empty());
        let fit = m.compounds[0].fragments[0].fits[0].as_ref().unwrap();
        for (got, want) in fit.mids.iter().zip([0.5, 0.5, 0.0]) {
            assert!((got - want).abs() < 1e-9, "{:?}", fit.mids);
        }
    }

    #[test]
    fn too_wide_fragments_are_left_out() {
        let mut c = LibraryCompound::new("Ala");
        c.formulas.push((0, u32::MAX, "C2".into()));
        assert!(fragments(&c).is_empty());
        assert!(mid_library(&[c]).is_empty());
    }
}
