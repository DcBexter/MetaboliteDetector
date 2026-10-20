//! Rescan of missing quant ions in the ion traces (old `ResultSeries::scanIons` called from
//! `preciseSeriesCompletion`): where a substance seen in other runs lacks a quant ion, look for that
//! ion's peak near the expected retention time and integrate it. Pinned behaviour and what is not
//! decoded: docs/reverse-engineering/batch-quantification.md.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};

use super::{
    Found, IonPeak, QuantSettings, RunResult, ScanSettings, compound_ions, compound_of,
    hits_by_compound,
};
use crate::Chromatogram;
use crate::deconv::{DeconvSettings, deconvolute};
use crate::library::LibraryCompound;
use crate::peaks::{IonTraces, TooLarge};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RescanSettings {
    /// Width of the search window in seconds, added around the expected retention times (the old
    /// window is the expected times ± half of it).
    pub window_s: f64,
    /// Only substances found in at least this share of the runs are completed (the old code skips a
    /// series when too few samples have an integral).
    pub min_share: f64,
}

impl Default for RescanSettings {
    fn default() -> Self {
        Self {
            // shortcut: the old width argument is not traced to its source (batch-quantification.md).
            window_s: 6.0,
            min_share: 0.5,
        }
    }
}

/// One substance of one run with quant ions to look for, and where.
struct Todo {
    compound: usize,
    ions: Vec<u32>,
    /// Search window in the run's retention times.
    from_s: f64,
    to_s: f64,
}

/// Ion peaks found for one substance: (compound, ions, mean apex time).
type Rescued = (usize, Vec<IonPeak>, f64);

/// Looks for the quant ions each substance lacks in a run (see the module doc) and adds them to the
/// run's `found`. Substances absent from a run get a `Found` with `rescanned` set. Runs are read
/// again, at most `parallel` at a time, and only if they have gaps. Returns the number of ions added
/// and the files that could not be read (their runs stay as they were), or `None` once `cancel` is set.
/// The quant ions come from `q`, as `table` would choose them with the same settings.
#[allow(clippy::too_many_arguments)]
pub fn rescan_runs(
    runs: &mut [RunResult],
    lib: &[LibraryCompound],
    extra: &[LibraryCompound],
    s: &ScanSettings,
    q: &QuantSettings,
    r: &RescanSettings,
    parallel: usize,
    cancel: &AtomicBool,
) -> Option<(usize, Vec<String>)> {
    let work: Vec<(usize, Vec<Todo>)> = plan(runs, lib, extra, q, r).into_iter().collect();
    let paths: Vec<_> = runs.iter().map(|r| r.path.clone()).collect();
    let next = AtomicUsize::new(0);
    let mut done: Vec<(usize, Result<Vec<Rescued>, String>)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..parallel.clamp(1, work.len().max(1)))
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    while !cancel.load(Relaxed) {
                        let Some((i, todo)) = work.get(next.fetch_add(1, Relaxed)) else {
                            break;
                        };
                        let path = &paths[*i];
                        // A panic in one file becomes its error, so the other files are kept.
                        let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            crate::read_run(path)
                                .and_then(|c| rescan_run(c, todo, s).map_err(|e| e.to_string()))
                        }))
                        .unwrap_or_else(|_| Err("analysis crashed".into()))
                        .map_err(|e| format!("{}: {e}", path.display()));
                        mine.push((*i, out));
                    }
                    mine
                })
            })
            .collect();
        workers
            .into_iter()
            .flat_map(|w| w.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
            .collect()
    });
    if cancel.load(Relaxed) {
        return None;
    }
    done.sort_by_key(|d| d.0);
    let (mut added, mut errors) = (0, Vec::new());
    for (i, out) in done {
        match out {
            Ok(rescued) => {
                for (compound, mut ions, rt_s) in rescued {
                    added += ions.len();
                    let found = &mut runs[i].found;
                    let k = found
                        .iter()
                        .position(|f| f.compound == compound)
                        .unwrap_or_else(|| {
                            found.push(Found {
                                compound,
                                score: 0.0,
                                rt_s,
                                ions: Vec::new(),
                                rescanned: true,
                            });
                            found.len() - 1
                        });
                    let f = &mut found[k];
                    // A zero-integral peak of the same ion is replaced.
                    f.ions.retain(|p| ions.iter().all(|n| n.mass != p.mass));
                    f.ions.append(&mut ions);
                    f.ions.sort_by(|a, b| b.integral.total_cmp(&a.integral));
                }
            }
            Err(e) => errors.push(e),
        }
    }
    Some((added, errors))
}

/// Per run (by index) the substances with missing quant ions.
fn plan(
    runs: &[RunResult],
    lib: &[LibraryCompound],
    extra: &[LibraryCompound],
    q: &QuantSettings,
    r: &RescanSettings,
) -> BTreeMap<usize, Vec<Todo>> {
    let refs: Vec<&RunResult> = runs.iter().collect();
    let mut todo: BTreeMap<usize, Vec<Todo>> = BTreeMap::new();
    for (compound, found) in hits_by_compound(&refs) {
        let seen = found.iter().flatten().count();
        if (seen as f64) < r.min_share * runs.len() as f64 {
            continue;
        }
        let Some(ions) = compound_ions(compound_of(lib, extra, compound), &found, q) else {
            continue;
        };
        for (i, run) in runs.iter().enumerate() {
            let has = |m: &u32| {
                found[i].is_some_and(|f| f.ions.iter().any(|p| p.mass == *m && p.integral > 0.0))
            };
            let missing: Vec<u32> = ions.iter().copied().filter(|m| !has(m)).collect();
            if missing.is_empty() {
                continue;
            }
            // Expected times in this run: the others' times, through the RI when both are calibrated.
            let times = found.iter().zip(runs).filter_map(|(f, other)| {
                let f = (*f)?;
                Some(match (&run.calibration, &other.calibration) {
                    (Some(mine), Some(theirs)) => mine.rt(theirs.ri(f.rt_s)),
                    _ => f.rt_s,
                })
            });
            let (lo, hi) = times.fold((f64::MAX, f64::MIN), |(lo, hi), t| (lo.min(t), hi.max(t)));
            todo.entry(i).or_default().push(Todo {
                compound,
                ions: missing,
                from_s: lo - r.window_s / 2.0,
                to_s: hi + r.window_s / 2.0,
            });
        }
    }
    todo
}

fn rescan_run(c: Chromatogram, todo: &[Todo], s: &ScanSettings) -> Result<Vec<Rescued>, TooLarge> {
    let mut traces = IonTraces::new(&c)?;
    let times = c.scan_times_s.clone();
    drop(c);
    let peaks = traces.detect_peaks(&s.peaks);
    let at = |k: usize| times[peaks[k].apex_scan as usize];
    let mut by_row: BTreeMap<u32, Vec<usize>> = BTreeMap::new();
    for (k, p) in peaks.iter().enumerate() {
        by_row.entry(p.row).or_default().push(k);
    }
    // At least one peak per component and twice the width, as in the old rescan (§4.0).
    let local = DeconvSettings {
        min_peak_number: 1,
        width: s.deconv.width * 2.0,
        ..s.deconv.clone()
    };
    let mut out = Vec::new();
    for t in todo {
        let middle = (t.from_s + t.to_s) / 2.0;
        // The peak of each missing ion nearest to the middle of the window.
        let picked: Vec<(u32, usize)> = t
            .ions
            .iter()
            .filter_map(|m| {
                let row = traces.masses.binary_search(m).ok()? as u32;
                let k = by_row
                    .get(&row)?
                    .iter()
                    .copied()
                    .filter(|k| (t.from_s..=t.to_s).contains(&at(*k)))
                    .min_by(|a, b| (at(*a) - middle).abs().total_cmp(&(at(*b) - middle).abs()))?;
                Some((*m, k))
            })
            .collect();
        if picked.is_empty() {
            continue;
        }
        let local_peaks: Vec<_> = picked.iter().map(|(_, k)| peaks[*k]).collect();
        let comps = deconvolute(&traces, &local_peaks, &local);
        let ions: Vec<IonPeak> = picked
            .iter()
            .zip(&local_peaks)
            .enumerate()
            .map(|(j, ((mass, _), p))| {
                let fitted = comps
                    .iter()
                    .flat_map(|c| &c.peaks)
                    .filter(|f| f.peak == j)
                    .map(|f| f.integral)
                    .fold(0.0, f64::max);
                // Without a fit the old code still reports "additional peak areas": the raw sum.
                let integral = if fitted > 0.0 {
                    fitted
                } else {
                    traces.trace(p.row as usize)[p.first_scan as usize..=p.last_scan as usize]
                        .iter()
                        .map(|v| *v as f64)
                        .sum()
                };
                IonPeak {
                    mass: *mass,
                    integral,
                    quality: p.quality,
                    signal_to_noise: p.signal_to_noise,
                }
            })
            .collect();
        let rt = picked.iter().map(|(_, k)| at(*k)).sum::<f64>() / picked.len() as f64;
        out.push((t.compound, ions, rt));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ri::RiCalibration;

    #[test]
    fn window_follows_the_calibration_of_each_run() {
        let mut lib = LibraryCompound::new("a");
        lib.quant_ions = vec![50.0];
        let found = Found {
            compound: 0,
            score: 1.0,
            rt_s: 100.0,
            ions: Vec::new(),
            rescanned: false,
        };
        let shifted = RiCalibration::identity().shifted_to(120.0, 100.0);
        let run = |found: Vec<Found>, calibration| RunResult {
            found,
            calibration,
            ..Default::default()
        };
        let r = RescanSettings::default();
        let q = QuantSettings::default();
        let window = |runs: &[RunResult]| {
            let t = &plan(runs, std::slice::from_ref(&lib), &[], &q, &r)[&1][0];
            (t.from_s, t.to_s)
        };
        let id = Some(RiCalibration::identity());
        let a = run(vec![found.clone()], id.clone());
        // The same substance elutes 20 s later in a run calibrated 20 s later.
        assert_eq!(
            window(&[a.clone(), run(vec![], Some(shifted))]),
            (117.0, 123.0)
        );
        assert_eq!(window(&[a, run(vec![], id)]), (97.0, 103.0));
        // Without calibrations the retention times are taken as they are.
        assert_eq!(
            window(&[run(vec![found], None), run(vec![], None)]),
            (97.0, 103.0)
        );
    }
}
