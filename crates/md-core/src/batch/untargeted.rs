//! Untargeted mode: the components no library entry claimed form a temporary library over all runs.
//! shortcut: the old `LibrarySearch::addCompoundSpectra` is not decoded; a group keeps the spectrum of
//! its first member and joins a later component by the library search score (spectrum, plus RI when
//! calibrated) at the library's cutoff, at most one component per run. Upgrade when the old merging is read.

use super::{Found, RunResult, ScanSettings};
use crate::library::LibraryCompound;
use crate::search::{PreparedLibrary, library_hits};

/// Pools the `unknown` components of all runs, in run order. Returns the temporary library; each
/// component becomes a `Found` of the run with index `first_index + k` for its group `k`
/// (`first_index` is the length of the reference library, so `table` finds the group in `extra`).
pub fn add_unknowns(
    runs: &mut [RunResult],
    first_index: usize,
    s: &ScanSettings,
) -> Vec<LibraryCompound> {
    let mut lib: Vec<LibraryCompound> = Vec::new();
    let mut prepared = PreparedLibrary::new(&[], &s.search);
    for run in runs {
        let mut taken = vec![false; lib.len()];
        for u in std::mem::take(&mut run.unknown) {
            let hit = library_hits(&u.spectrum, u.ri, &prepared, &s.search)
                .into_iter()
                .find(|h| !taken[h.compound]);
            let (k, score) = match hit {
                Some(h) => (h.compound, h.overall),
                None => {
                    let mut c = LibraryCompound::new(&format!("Unknown {}", lib.len() + 1));
                    c.spectrum.rt_ms = u.rt_s * 1000.0;
                    c.spectrum.ri = u.ri.unwrap_or(-1.0);
                    c.spectrum.masses = u.spectrum.iter().map(|p| p.0 as f32).collect();
                    c.spectrum.values = u.spectrum.iter().map(|p| p.1 as i32).collect();
                    prepared.push(&c, &s.search);
                    lib.push(c);
                    taken.push(false);
                    (lib.len() - 1, 1.0)
                }
            };
            taken[k] = true;
            run.found.push(Found {
                compound: first_index + k,
                score,
                rt_s: u.rt_s,
                ions: u.ions,
                rescanned: false,
            });
        }
    }
    lib
}
