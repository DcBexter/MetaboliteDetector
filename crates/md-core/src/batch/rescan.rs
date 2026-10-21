//! Rescan of missing quant ions in the ion traces (old `ResultSeries::scanIons` called from
//! `preciseSeriesCompletion`): where a substance seen in other runs lacks a quant ion, look for that
//! ion's peak near the expected retention time and integrate it. The window is the search's ΔRI wide
//! (old `LS_RI_DIFF`); pinned behaviour: docs/reverse-engineering/batch-quantification.md.

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
use crate::ri::RiCalibration;

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
/// again, at most `parallel` at a time, and only if they have gaps; `.cmp` inputs are skipped. Returns the number of ions
/// added and the files that could not be read (their runs stay as they were), or `None` once `cancel` is set.
/// The quant ions come from `q`, as `table` would choose them with the same settings.
pub fn rescan_runs(
    runs: &mut [RunResult],
    lib: &[LibraryCompound],
    extra: &[LibraryCompound],
    s: &ScanSettings,
    q: &QuantSettings,
    parallel: usize,
    cancel: &AtomicBool,
) -> Option<(usize, Vec<String>)> {
    // A `.cmp` input has no raw data to look into.
    let work: Vec<(usize, Vec<Todo>)> = (plan(runs, lib, extra, q, s.search.ri_diff).into_iter())
        .filter(|(i, _)| !crate::cmp::is_cmp(&runs[*i].path))
        .collect();
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
                                rescued: Vec::new(),
                                alternatives: Vec::new(),
                            });
                            found.len() - 1
                        });
                    let f = &mut found[k];
                    // A zero-integral peak of the same ion is replaced.
                    f.ions.retain(|p| ions.iter().all(|n| n.mass != p.mass));
                    if !f.rescanned {
                        f.rescued.extend(ions.iter().map(|p| p.mass));
                    }
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
    ri_diff: f64,
) -> BTreeMap<usize, Vec<Todo>> {
    let refs: Vec<&RunResult> = runs.iter().collect();
    let mut todo: BTreeMap<usize, Vec<Todo>> = BTreeMap::new();
    for (compound, found) in hits_by_compound(&refs) {
        // Old `preciseSeriesCompletion`: skipped when fewer than n/4 (rounded) samples have an integral.
        let seen = (found.iter().flatten())
            .filter(|f| f.ions.iter().any(|p| p.integral > 0.0))
            .count();
        if seen < (runs.len() as f64 / 4.0 + 0.5) as usize {
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
            let (from_s, to_s) = window(&run.calibration, &found, runs, ri_diff);
            todo.entry(i).or_default().push(Todo {
                compound,
                ions: missing,
                from_s,
                to_s,
            });
        }
    }
    todo
}

/// Old window of `scanIons`: the series' positions ± ΔRI/2. In RI when the run is calibrated and the
/// series has RIs (those of the runs with a calibration), else in retention time, ΔRI s wide.
fn window(
    mine: &Option<RiCalibration>,
    found: &[Option<&Found>],
    runs: &[RunResult],
    ri_diff: f64,
) -> (f64, f64) {
    let others = found.iter().zip(runs).filter_map(|(f, r)| Some(((*f)?, r)));
    let ris: Vec<f64> = (others.clone())
        .filter_map(|(f, r)| Some(r.calibration.as_ref()?.ri(f.rt_s)))
        .filter(|ri| *ri >= 0.0)
        .collect();
    match mine {
        Some(cal) if !ris.is_empty() => {
            let (lo, hi) = span(&ris);
            (cal.rt(lo - ri_diff / 2.0), cal.rt(hi + ri_diff / 2.0))
        }
        _ => {
            let rts: Vec<f64> = others.map(|(f, _)| f.rt_s).filter(|t| *t >= 0.0).collect();
            let (lo, hi) = span(&rts);
            (lo - ri_diff / 2.0, hi + ri_diff / 2.0)
        }
    }
}

fn span(v: &[f64]) -> (f64, f64) {
    v.iter()
        .fold((f64::MAX, f64::MIN), |(l, h), x| (l.min(*x), h.max(*x)))
}

fn rescan_run(c: Chromatogram, todo: &[Todo], s: &ScanSettings) -> Result<Vec<Rescued>, TooLarge> {
    let mut traces = IonTraces::new(&c, s.peaks.fwhm)?;
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
        width: Some(s.deconv.width_scans(&times) * 2.0),
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
                // Every trace of that nominal mass (several on exact-mass runs): the old code keeps
                // the peaks within ±0.5 of the ion.
                let rows = traces.masses.partition_point(|x| x < m)
                    ..traces.masses.partition_point(|x| x <= m);
                let k = rows
                    .filter_map(|r| by_row.get(&(r as u32)))
                    .flatten()
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
        let comps = deconvolute(&traces, &local_peaks, &times, &local);
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
    fn exact_mass_runs_search_every_trace_of_the_ion() {
        // Three traces of nominal mass 73; only the first has a peak.
        let (scans, mz) = (100, [73.047f32, 73.08, 73.09]);
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| s * mz.len()).collect(),
            mz: (0..scans).flat_map(|_| mz).collect(),
            intensity: (0..scans * mz.len())
                .map(|k| {
                    let (s, i) = ((k / mz.len()) as f32, k % mz.len());
                    let noise = ((k as u32).wrapping_mul(2_654_435_761) >> 16) % 41;
                    let peak = if i == 0 {
                        50000.0 * (-0.5 * ((s - 50.0) / 2.5).powi(2)).exp()
                    } else {
                        0.0
                    };
                    980.0 + noise as f32 + peak
                })
                .collect(),
            exact_masses: true,
        };
        let todo = Todo {
            compound: 0,
            ions: vec![73],
            from_s: 40.0,
            to_s: 60.0,
        };
        let out = rescan_run(c, &[todo], &ScanSettings::default()).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1[0].mass, 73);
        assert!((out[0].2 - 50.0).abs() <= 1.0, "{}", out[0].2);
    }

    #[test]
    fn window_follows_the_calibration_of_each_run() {
        let mut lib = LibraryCompound::new("a");
        lib.quant_ions = vec![50.0];
        let found = Found {
            compound: 0,
            score: 1.0,
            rt_s: 100.0,
            ions: vec![IonPeak {
                mass: 50,
                integral: 1.0,
                quality: 0.0,
                signal_to_noise: 1.0,
            }],
            rescanned: false,
            rescued: Vec::new(),
            alternatives: Vec::new(),
        };
        let shifted = RiCalibration::identity().shifted_to(120.0, 100.0);
        let run = |found: Vec<Found>, calibration| RunResult {
            found,
            calibration,
            ..Default::default()
        };
        let q = QuantSettings::default();
        let window = |runs: &[RunResult]| {
            let t = &plan(runs, std::slice::from_ref(&lib), &[], &q, 20.0)[&1][0];
            (t.from_s, t.to_s)
        };
        let id = Some(RiCalibration::identity());
        let a = run(vec![found.clone()], id.clone());
        // The window is ΔRI (20) wide in RI; the same substance elutes 20 s later in a run calibrated
        // 20 s later.
        assert_eq!(
            window(&[a.clone(), run(vec![], Some(shifted))]),
            (110.0, 130.0)
        );
        assert_eq!(window(&[a, run(vec![], id)]), (90.0, 110.0));
        // Without calibrations the retention times are taken as they are, ΔRI seconds wide.
        assert_eq!(
            window(&[run(vec![found], None), run(vec![], None)]),
            (90.0, 110.0)
        );
    }
}
