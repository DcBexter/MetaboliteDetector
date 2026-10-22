//! Rescan of missing quant ions in the ion traces (old `ResultSeries::scanIons` called from
//! `preciseSeriesCompletion`): where a substance seen in other runs lacks a quant ion, look for that
//! ion's peak near the expected retention time and integrate it. The window is the search's ΔRI wide
//! (old `LS_RI_DIFF`); pinned behaviour: docs/reverse-engineering/batch-quantification.md.

use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};

use crate::job::Cancel;

use super::{
    Found, IonPeak, QuantSettings, RunResult, ScanSettings, compound_ions, compound_of,
    hits_by_compound,
};
use crate::Chromatogram;
use crate::analysis::{AnalyzeError, PeakPass, RunPeaks, find_run_peaks};
use crate::deconv::{DeconvSettings, deconvolute};
use crate::library::LibraryCompound;
use crate::ri::RiCalibration;

/// One substance of one run with quant ions to look for, and where.
struct Todo {
    compound: usize,
    ions: Vec<u32>,
    /// Search window in the run's retention times.
    from_s: f64,
    to_s: f64,
    /// Mean position of the series: the old code takes the peak nearest to it. An RI when
    /// `calibration` is set, else a retention time.
    expected: f64,
    calibration: Option<RiCalibration>,
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
    cancel: &Cancel,
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
                    while !cancel.is_cancelled() {
                        let Some((i, todo)) = work.get(next.fetch_add(1, Relaxed)) else {
                            break;
                        };
                        let path = &paths[*i];
                        // A panic in one file becomes its error, so the other files are kept.
                        let out = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            super::read_scan_run(path, s).and_then(|c| {
                                rescan_run(c, todo, s, cancel).map_err(|e| e.to_string())
                            })
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
    if cancel.is_cancelled() {
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
            let (from_s, to_s, expected, calibration) =
                window(&run.calibration, &found, runs, ri_diff);
            todo.entry(i).or_default().push(Todo {
                compound,
                ions: missing,
                from_s,
                to_s,
                expected,
                calibration,
            });
        }
    }
    todo
}

/// Old window of `scanIons` in the run's retention times: the series' positions ± ΔRI/2, and their
/// mean. In RI when the run is calibrated and the series has RIs (those of the runs with a
/// calibration; the mean is then an RI and the calibration is returned), else in retention time,
/// ΔRI s wide.
fn window(
    mine: &Option<RiCalibration>,
    found: &[Option<&Found>],
    runs: &[RunResult],
    ri_diff: f64,
) -> (f64, f64, f64, Option<RiCalibration>) {
    let others = found.iter().zip(runs).filter_map(|(f, r)| Some(((*f)?, r)));
    let ris: Vec<f64> = (others.clone())
        .filter_map(|(f, r)| Some(r.calibration.as_ref()?.ri(f.rt_s)))
        .filter(|ri| *ri >= 0.0)
        .collect();
    match mine {
        Some(cal) if !ris.is_empty() => {
            let (lo, hi, mean) = span(&ris);
            let (from, to) = (cal.rt(lo - ri_diff / 2.0), cal.rt(hi + ri_diff / 2.0));
            (from, to, mean, Some(cal.clone()))
        }
        _ => {
            let rts: Vec<f64> = others.map(|(f, _)| f.rt_s).filter(|t| *t >= 0.0).collect();
            let (lo, hi, mean) = span(&rts);
            (lo - ri_diff / 2.0, hi + ri_diff / 2.0, mean, None)
        }
    }
}

/// Lowest, highest and mean value.
fn span(v: &[f64]) -> (f64, f64, f64) {
    let (lo, hi) = (v.iter()).fold((f64::MAX, f64::MIN), |(l, h), x| (l.min(*x), h.max(*x)));
    (lo, hi, v.iter().sum::<f64>() / v.len() as f64)
}

fn rescan_run(
    c: Chromatogram,
    todo: &[Todo],
    s: &ScanSettings,
    cancel: &Cancel,
) -> Result<Vec<Rescued>, AnalyzeError> {
    // Fresh, not Reordered: the old rescan (`ResultSeries::scanIons`) builds a fresh
    // `DefaultPeakDetector` and never calls `reorderPeaks` (reverse-engineering/algorithms.md:268).
    let RunPeaks {
        traces,
        times,
        peaks,
    } = find_run_peaks(c, &s.peaks, PeakPass::Fresh, cancel)?;
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
        // Old `Compound_Compare`: nearest to the series' mean position (in RI on a calibrated run),
        // the higher peak on a tie; no limit on the distance besides the window
        // (batch-quantification.md).
        let pos = |k: usize| t.calibration.as_ref().map_or(at(k), |c| c.ri(at(k)));
        let dist = |k: usize| (pos(k) - t.expected).abs();
        let height = |k: usize| traces.trace(peaks[k].row as usize)[peaks[k].apex_scan as usize];
        let (mut ions, mut rts) = (Vec::new(), Vec::new());
        for m in &t.ions {
            // Every trace of that nominal mass (several on exact-mass runs): the old code keeps the
            // peaks within ±0.5 of the ion.
            let rows =
                traces.masses.partition_point(|x| x < m)..traces.masses.partition_point(|x| x <= m);
            let found: Vec<usize> = (rows.filter_map(|r| by_row.get(&(r as u32))).flatten())
                .copied()
                .filter(|k| (t.from_s..=t.to_s).contains(&at(*k)))
                .collect();
            if found.is_empty() {
                continue;
            }
            // The old code deconvolutes the window's peaks of each ion on their own (in detection
            // order) and takes the nearest component that has the ion.
            let local_peaks: Vec<_> = found.iter().map(|k| peaks[*k]).collect();
            let comps = deconvolute(&traces, &local_peaks, &times, &local);
            let fitted = |j: usize| {
                (comps.iter().flat_map(|c| &c.peaks))
                    .filter(|f| f.peak == j)
                    .map(|f| f.integral)
                    .fold(0.0, f64::max)
            };
            let mut order: Vec<usize> = (0..found.len()).collect();
            order.sort_by(|a, b| {
                let (a, b) = (found[*a], found[*b]);
                (dist(a).total_cmp(&dist(b))).then(height(b).total_cmp(&height(a)))
            });
            let (j, integral) = (order.iter())
                .map(|j| (*j, fitted(*j)))
                .find(|(_, v)| *v > 0.0)
                .unwrap_or_else(|| {
                    // No fit at all: the old code still reports "additional peak areas", the raw sum
                    // of the nearest peak.
                    let p = &local_peaks[order[0]];
                    let raw = traces.trace(p.row as usize)
                        [p.first_scan as usize..=p.last_scan as usize]
                        .iter()
                        .map(|v| *v as f64)
                        .sum();
                    (order[0], raw)
                });
            let p = &local_peaks[j];
            ions.push(IonPeak {
                mass: *m,
                integral,
                quality: p.quality,
                signal_to_noise: p.signal_to_noise,
            });
            rts.push(at(found[j]));
        }
        if ions.is_empty() {
            continue;
        }
        out.push((t.compound, ions, rts.iter().sum::<f64>() / rts.len() as f64));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::peaks::IonTraces;
    use crate::ri::RiCalibration;

    /// A run of `scans` scans at 1 s with traces of masses `mz` and the given Gaussian peaks
    /// (trace, apex scan, height) on a noisy baseline.
    fn run(scans: usize, mz: &[f32], peaks: &[(usize, f32, f32)]) -> Chromatogram {
        Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64).collect(),
            offsets: (0..=scans).map(|s| s * mz.len()).collect(),
            mz: (0..scans).flat_map(|_| mz.iter().copied()).collect(),
            intensity: (0..scans * mz.len())
                .map(|k| {
                    let (s, i) = ((k / mz.len()) as f32, k % mz.len());
                    let noise = ((k as u32).wrapping_mul(2_654_435_761) >> 16) % 41;
                    let peak: f32 = (peaks.iter())
                        .filter(|p| p.0 == i)
                        .map(|p| p.2 * (-0.5 * ((s - p.1) / 2.5).powi(2)).exp())
                        .sum();
                    980.0 + noise as f32 + peak
                })
                .collect(),
            exact_masses: true,
            baseline: Default::default(),
        }
    }

    #[test]
    fn exact_mass_runs_search_every_trace_of_the_ion() {
        // Three traces of nominal mass 73; only the first has a peak.
        let c = run(100, &[73.047, 73.08, 73.09], &[(0, 50.0, 50000.0)]);
        let todo = Todo {
            compound: 0,
            ions: vec![73],
            from_s: 40.0,
            to_s: 60.0,
            expected: 50.0,
            calibration: None,
        };
        let out = rescan_run(c, &[todo], &ScanSettings::default(), &Cancel::default()).unwrap();
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].1[0].mass, 73);
        assert!((out[0].2 - 50.0).abs() <= 1.0, "{}", out[0].2);
    }

    #[test]
    fn takes_the_peak_nearest_to_the_mean_position() {
        // #336: series at 30, 32 and 40 s give the window 20..50 s with its middle at 35 s, but the
        // old `scanIons` sorts by the distance to the mean, 34 s: the peak at 31 s wins over 38 s.
        let todo = || Todo {
            compound: 0,
            ions: vec![73],
            from_s: 20.0,
            to_s: 50.0,
            expected: 34.0,
            calibration: None,
        };
        let c = run(80, &[73.0], &[(0, 31.0, 50000.0), (0, 38.0, 50000.0)]);
        let out = rescan_run(c, &[todo()], &ScanSettings::default(), &Cancel::default()).unwrap();
        assert!((out[0].2 - 31.0).abs() <= 1.0, "{}", out[0].2);
        // On equal distance the higher peak wins.
        let c = run(80, &[73.0], &[(0, 30.0, 20000.0), (0, 38.0, 50000.0)]);
        let out = rescan_run(c, &[todo()], &ScanSettings::default(), &Cancel::default()).unwrap();
        assert!((out[0].2 - 38.0).abs() <= 1.0, "{}", out[0].2);
    }

    #[test]
    fn calibrated_runs_rank_peaks_by_ri() {
        // 1 RI/s up to 30 s, then 10 RI/s. Expected RI 30 (30 s): the peak at 20 s (RI 20) is
        // nearer in RI than the one at 32 s (RI 50), though farther in seconds and from the window's
        // middle.
        let cal = RiCalibration::new([(0.0, 0.0), (30.0, 30.0), (60.0, 330.0)]).unwrap();
        let c = run(80, &[73.0], &[(0, 20.0, 50000.0), (0, 32.0, 50000.0)]);
        let todo = Todo {
            compound: 0,
            ions: vec![73],
            from_s: 15.0,
            to_s: 45.0,
            expected: 30.0,
            calibration: Some(cal),
        };
        let out = rescan_run(c, &[todo], &ScanSettings::default(), &Cancel::default()).unwrap();
        assert!((out[0].2 - 20.0).abs() <= 1.0, "{}", out[0].2);
    }

    #[test]
    fn peak_without_fit_falls_to_the_next_peak_with_fit() {
        // Two traces of mass 73 on an exact-mass run, one component: the weak peak at 41 s is the
        // nearest but gets no fit next to the strong one at 42 s, so the strong one is taken (#344).
        let todo = Todo {
            compound: 0,
            ions: vec![73],
            from_s: 30.0,
            to_s: 50.0,
            expected: 41.0,
            calibration: None,
        };
        let s = ScanSettings::default();
        // Alone, the weak peak has a fit and wins.
        let c = run(80, &[73.047, 73.08], &[(1, 41.0, 3000.0)]);
        let out = rescan_run(c, std::slice::from_ref(&todo), &s, &Cancel::default()).unwrap();
        assert!((out[0].2 - 41.0).abs() <= 0.5, "{}", out[0].2);
        let c = run(
            80,
            &[73.047, 73.08],
            &[(0, 42.0, 50000.0), (1, 41.0, 3000.0)],
        );
        let out = rescan_run(c, &[todo], &s, &Cancel::default()).unwrap();
        assert!((out[0].2 - 42.0).abs() <= 0.5, "{}", out[0].2);
        assert!(out[0].1[0].integral > 100_000.0, "{}", out[0].1[0].integral);
    }

    #[test]
    fn no_fit_at_all_reports_the_raw_sum_of_the_nearest_peak() {
        let c = run(80, &[73.0], &[(0, 40.0, 50000.0)]);
        // A filter wider than the run finds no component, so no peak has a fit.
        let mut s = ScanSettings::default();
        s.deconv.width = Some(25.0);
        let mut traces = IonTraces::new(&c, s.peaks.fwhm).unwrap();
        let peaks = traces.detect_peaks(&s.peaks);
        assert_eq!(peaks.len(), 1);
        let local = DeconvSettings {
            min_peak_number: 1,
            width: Some(50.0),
            ..s.deconv.clone()
        };
        assert!(deconvolute(&traces, &peaks, &c.scan_times_s, &local).is_empty());
        let p = &peaks[0];
        let raw: f64 = traces.trace(0)[p.first_scan as usize..=p.last_scan as usize]
            .iter()
            .map(|v| *v as f64)
            .sum();
        let todo = Todo {
            compound: 0,
            ions: vec![73],
            from_s: 30.0,
            to_s: 50.0,
            expected: 40.0,
            calibration: None,
        };
        let out = rescan_run(c, &[todo], &s, &Cancel::default()).unwrap();
        assert_eq!(out[0].1[0].integral, raw);
        assert!((out[0].2 - 40.0).abs() <= 0.5, "{}", out[0].2);
    }

    #[test]
    fn glycine_177_takes_the_peak_the_old_detector_finds() {
        // #327/#336: Glycine 3TMS elutes at 313 s; without RI calibration the default ΔRI 20 gives the
        // window 302.9..323.1 s. In 7235eg04 the old `DefaultPeakDetector` finds m/z 177 inside it only
        // at apex scan 2427, 321 s (crates/md-core/tests/data/old-md-3.3.tsv), the peak of another
        // substance, and `scanIons` has no other limit on the distance: it takes that peak, so do we.
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../example-data/TargetSearchData/inst/gc-ms-data/7235eg04.cdf");
        let Ok(c) = crate::read_run(&path) else {
            eprintln!("skipped: {} not found", path.display());
            return;
        };
        let apex = c.scan_times_s[2427];
        let todo = Todo {
            compound: 0,
            ions: vec![177],
            from_s: 302.9,
            to_s: 323.1,
            expected: 313.0,
            calibration: None,
        };
        let out = rescan_run(c, &[todo], &ScanSettings::default(), &Cancel::default()).unwrap();
        assert!((out[0].2 - apex).abs() < 0.11, "{} vs {apex}", out[0].2);
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
            (t.from_s, t.to_s, t.expected)
        };
        let id = Some(RiCalibration::identity());
        let a = run(vec![found.clone()], id.clone());
        // The window is ΔRI (20) wide in RI; the same substance elutes 20 s later in a run calibrated
        // 20 s later. The expected position stays an RI there (100).
        assert_eq!(
            window(&[a.clone(), run(vec![], Some(shifted))]),
            (110.0, 130.0, 100.0)
        );
        assert_eq!(window(&[a, run(vec![], id)]), (90.0, 110.0, 100.0));
        // Without calibrations the retention times are taken as they are, ΔRI seconds wide.
        assert_eq!(
            window(&[run(vec![found], None), run(vec![], None)]),
            (90.0, 110.0, 100.0)
        );
    }
}
