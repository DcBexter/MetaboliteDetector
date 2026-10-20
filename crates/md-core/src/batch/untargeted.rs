//! Untargeted mode: the components no library entry claimed form a temporary library over all runs
//! (old `BatchQuantificator::createLibrary` / `LibrarySearch::addCompoundSpectra`, see
//! docs/reverse-engineering/batch-quantification.md).
//!
//! Each run is matched as a whole against the library as it stood before the run
//! (`mapCompoundsToLibrary`, one component per library entry, best score first). A matched entry
//! takes the component's spectrum into its mean; the rest become new entries.
//! shortcut: the old code builds one library per replicate group and merges those; the port pools all
//! runs into one library and leaves the groups to `table`. Upgrade if group-wise libraries matter.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, Ordering::Relaxed};

use super::{Found, RunResult, ScanSettings, Unknown};
use crate::library::LibraryCompound;
use crate::search::{PreparedLibrary, assign_best_first, library_hits};

/// Pools the `unknown` components of all runs, in run order. Returns the temporary library; each
/// component becomes a `Found` of the run with index `first_index + k` for its group `k`
/// (`first_index` is the length of the reference library, so `table` finds the group in `extra`).
/// `None` once `cancel` is set; the scoring grows with the square of the number of unknowns.
pub fn add_unknowns(
    runs: &mut [RunResult],
    first_index: usize,
    s: &ScanSettings,
    cancel: &AtomicBool,
) -> Option<Vec<LibraryCompound>> {
    let mut lib: Vec<LibraryCompound> = Vec::new();
    let mut prepared = PreparedLibrary::new(&[], &s.search);
    for run in runs {
        if cancel.load(Relaxed) {
            return None;
        }
        let unknown = std::mem::take(&mut run.unknown);
        let hits: Vec<_> = (unknown.iter())
            .map(|u| library_hits(&u.spectrum, u.ri, &prepared, &s.search))
            .collect();
        let assigned = assign_best_first(&hits, lib.len());
        for (u, hit) in unknown.into_iter().zip(assigned) {
            let (k, score) = match hit {
                Some(h) => {
                    merge(&mut lib[h.compound], &u);
                    (h.compound, h.overall)
                }
                None => {
                    lib.push(new_group(lib.len() + 1, &u));
                    (lib.len() - 1, 1.0)
                }
            };
            run.found.push(Found {
                compound: first_index + k,
                score,
                rt_s: u.rt_s,
                ions: u.ions,
                rescanned: false,
            });
        }
        prepared = PreparedLibrary::new(&lib, &s.search);
    }
    Some(lib)
}

fn new_group(n: usize, u: &Unknown) -> LibraryCompound {
    let mut c = LibraryCompound::new(&format!("Unknown {n}"));
    c.spectrum.rt_ms = u.rt_s * 1000.0;
    c.spectrum.ri = u.ri.unwrap_or(-1.0);
    c.spectrum.masses = u.spectrum.iter().map(|p| p.0 as f32).collect();
    c.spectrum.values = u.spectrum.iter().map(|p| p.1 as i32).collect();
    c.ion_freq = u.spectrum.iter().map(|p| (p.0 as f32, 1)).collect();
    c.spec_count = 1;
    c
}

/// Old `LibraryCompound::addCompoundSpectrum` for one more spectrum: per ion the mean over the
/// spectra that contain it (`calculateNewValue`, cut to an integer), the retention time and RI
/// weighted by the number of spectra.
fn merge(c: &mut LibraryCompound, u: &Unknown) {
    let mut ions: BTreeMap<u32, (i32, u32)> = (c.spectrum.masses.iter().zip(&c.spectrum.values))
        .zip(&c.ion_freq)
        .map(|((m, v), (_, n))| (*m as u32, (*v, *n)))
        .collect();
    for &(m, h) in &u.spectrum {
        let e = ions.entry(m).or_insert((0, 0));
        let (n, v) = (e.1 as f64, e.0 as f64);
        *e = (((n * v + h as f64) / (n + 1.0)) as i32, e.1 + 1);
    }
    c.spectrum.masses = ions.keys().map(|m| *m as f32).collect();
    c.spectrum.values = ions.values().map(|e| e.0).collect();
    c.ion_freq = ions.iter().map(|(m, e)| (*m as f32, e.1)).collect();
    let n = c.spec_count as f64;
    c.spectrum.rt_ms = (n * c.spectrum.rt_ms + u.rt_s * 1000.0) / (n + 1.0);
    // shortcut: the old handling of a missing RI (-1) is not decoded; the known one wins.
    match u.ri {
        Some(ri) if c.spectrum.ri >= 0.0 => c.spectrum.ri = (n * c.spectrum.ri + ri) / (n + 1.0),
        Some(ri) => c.spectrum.ri = ri,
        None => {}
    }
    c.spec_count += 1;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn unknown(rt_s: f64, spectrum: &[(u32, u32)]) -> Unknown {
        Unknown {
            rt_s,
            ri: None,
            spectrum: spectrum.to_vec(),
            ions: Vec::new(),
        }
    }

    fn run(unknown: Vec<Unknown>) -> RunResult {
        RunResult {
            path: PathBuf::from("r.cdf"),
            unknown,
            ..Default::default()
        }
    }

    #[test]
    fn a_group_is_the_mean_of_its_spectra_and_a_run_uses_it_once() {
        let a = [(50, 100), (60, 200), (70, 300)];
        let b = [(50, 300), (60, 200), (70, 100), (80, 40)];
        let mut runs = vec![
            run(vec![unknown(10.0, &a)]),
            // Both components resemble the group, the second one better: it joins, the first one
            // starts a group of its own (one component per entry, best score first).
            run(vec![
                unknown(20.0, &[(50, 300), (70, 100)]),
                unknown(20.0, &b),
            ]),
        ];
        let lib = add_unknowns(
            &mut runs,
            5,
            &ScanSettings::default(),
            &AtomicBool::new(false),
        )
        .unwrap();
        assert_eq!(lib.len(), 2);
        let groups: Vec<usize> = runs[1].found.iter().map(|f| f.compound).collect();
        assert_eq!(groups, [6, 5]);
        // Mean per ion over the spectra that contain it: 80 only came with one.
        let g = &lib[0];
        assert_eq!(g.spectrum.masses, [50.0, 60.0, 70.0, 80.0]);
        assert_eq!(g.spectrum.values, [200, 200, 200, 40]);
        assert_eq!((g.spec_count, g.spectrum.rt_ms), (2, 15_000.0));
        assert_eq!(
            g.ion_freq.iter().map(|f| f.1).collect::<Vec<_>>(),
            [2, 2, 2, 1]
        );
    }
}
