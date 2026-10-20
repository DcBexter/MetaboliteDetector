//! Batch quantification: finds the library substances in many runs and builds a substance × run
//! table (old `BatchQuantificator`, feature-inventory.md F10–F12). Quant ions follow algorithms.md §7.
//!
//! Two passes: `scan_runs` (expensive, once) keeps per run only the identified substances and their ion
//! integrals; `table` (cheap) turns them into the table, so filters and replicate groups can change
//! without rescanning. Targeted only: every substance comes from the reference library.
//! shortcut: no non-targeted mode and no SIC rescan of missing quant ions (old `ResultSeries::scanIons`,
//! §4.0); a substance not deconvoluted in a run stays empty there. Follow-up issue #58.

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};

use regex::Regex;

use crate::Chromatogram;
use crate::deconv::{Component, DeconvSettings, deconvolute};
use crate::library::LibraryCompound;
use crate::peaks::{IonTraces, Peak, PeakSettings, TooLarge};
use crate::ri::RiCalibration;
use crate::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};

/// Settings of the per-run pass: the single-run defaults, one calibration for every run.
#[derive(Debug, Clone, Default)]
pub struct ScanSettings {
    pub peaks: PeakSettings,
    pub deconv: DeconvSettings,
    pub search: SearchSettings,
    /// RT → RI for every run; without it the search compares spectra only.
    pub calibration: Option<RiCalibration>,
}

/// One ion peak of a found substance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IonPeak {
    pub mass: u32,
    /// Integral from the deconvolution (`FittedPeak::integral`).
    pub integral: f64,
    /// Peak quality, lower is better (§3.7).
    pub quality: f32,
}

/// A library substance assigned to a component of one run.
#[derive(Debug, Clone, PartialEq)]
pub struct Found {
    /// Index into the library.
    pub compound: usize,
    pub score: f64,
    pub rt_s: f64,
    /// Every ion peak of the component, largest integral first. A peak can sit twice (§4.5); §7.2
    /// sees both, as in the old program.
    pub ions: Vec<IonPeak>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct RunResult {
    pub path: PathBuf,
    pub found: Vec<Found>,
}

/// The ion peaks of a component, largest integral first (§7.2).
pub fn ion_peaks(k: &Component, peaks: &[Peak], traces: &IonTraces) -> Vec<IonPeak> {
    let mut ions: Vec<IonPeak> = k
        .peaks
        .iter()
        .map(|f| {
            let p = &peaks[f.peak];
            IonPeak {
                mass: traces.masses[p.row as usize],
                integral: f.integral,
                quality: p.quality,
            }
        })
        .collect();
    // The old code sorts ascending and walks backwards, so ties come last-in-list first (§7.2).
    ions.sort_by(|a, b| a.integral.total_cmp(&b.integral));
    ions.reverse();
    ions
}

/// Deconvolutes one run and returns the library substances assigned to its components.
pub fn scan_chromatogram(
    c: Chromatogram,
    lib: &PreparedLibrary,
    s: &ScanSettings,
) -> Result<Vec<Found>, TooLarge> {
    let mut traces = IonTraces::new(&c)?;
    // The ion-major copy and the times are all the rest needs; dropping the spectra now halves the
    // peak memory per run.
    let times = c.scan_times_s.clone();
    drop(c);
    let peaks = traces.detect_peaks(&s.peaks);
    let comps = deconvolute(&traces, &peaks, &s.deconv);
    let ris = s
        .calibration
        .as_ref()
        .map(|cal| cal.component_ris(&times, &comps));
    let id = map_compounds_to_library(&comps, ris.as_deref(), lib, &s.search);
    Ok(comps
        .iter()
        .zip(id.assigned)
        .filter_map(|(k, hit)| {
            let hit = hit?;
            let ions = ion_peaks(k, &peaks, &traces);
            Some(Found {
                compound: hit.compound,
                score: hit.overall,
                rt_s: times[k.apex_scan as usize],
                ions,
            })
        })
        .collect())
}

/// Scans every file with at most `parallel` runs in memory at once; each run is analysed in parallel
/// (rayon) as well. `done` counts finished files. Once `cancel` is set, no further file starts and
/// the result is `None`. Unreadable files are listed as errors and do not stop the batch.
pub fn scan_runs(
    paths: &[PathBuf],
    lib: &PreparedLibrary,
    s: &ScanSettings,
    parallel: usize,
    done: &AtomicUsize,
    cancel: &AtomicBool,
) -> Option<(Vec<RunResult>, Vec<String>)> {
    let next = AtomicUsize::new(0);
    let mut out: Vec<(usize, Result<RunResult, String>)> = std::thread::scope(|scope| {
        let workers: Vec<_> = (0..parallel.clamp(1, paths.len().max(1)))
            .map(|_| {
                scope.spawn(|| {
                    let mut mine = Vec::new();
                    while !cancel.load(Relaxed) {
                        let i = next.fetch_add(1, Relaxed);
                        let Some(path) = paths.get(i) else { break };
                        // A panic in one file becomes its error, so the other files are kept.
                        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                            crate::read_run(path).and_then(|c| {
                                scan_chromatogram(c, lib, s).map_err(|e| e.to_string())
                            })
                        }))
                        .unwrap_or_else(|_| Err("analysis crashed".into()))
                        .map(|found| RunResult {
                            path: path.clone(),
                            found,
                        })
                        .map_err(|e| format!("{}: {e}", path.display()));
                        mine.push((i, r));
                        done.fetch_add(1, Relaxed);
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
    out.sort_by_key(|r| r.0);
    let (mut runs, mut errors) = (Vec::new(), Vec::new());
    for (_, r) in out {
        match r {
            Ok(run) => runs.push(run),
            Err(e) => errors.push(e),
        }
    }
    Some((runs, errors))
}

/// The old `QUANT_*` settings; defaults are the old static values (file-formats.md §5).
#[derive(Debug, Clone, PartialEq)]
pub struct QuantSettings {
    /// `QUANT_ION_NUMBER`
    pub ion_number: usize,
    /// `QUANT_MIN_DIST_IONS`
    pub min_dist: f32,
    /// `QUANT_QUALITY_CUTOFF`: only peaks of lower (better) quality are picked first.
    pub quality_cutoff: f32,
    /// `QUANT_EXCL_IONS`: inclusive mass ranges never picked first.
    pub excluded: Vec<(f32, f32)>,
}

impl Default for QuantSettings {
    fn default() -> Self {
        Self {
            ion_number: 3,
            min_dist: 0.5,
            quality_cutoff: 1.0,
            excluded: Vec::new(),
        }
    }
}

/// Old `autoQuantificationIonAssignment` (§7.2) for ion peaks sorted by integral, largest first.
/// Returns the ions ascending. Keeps the old fallback quirk: without enough ions passing the filters
/// it restarts at the largest peak without filters, and an ion already picked still counts.
pub fn quant_ions(ions: &[IonPeak], q: &QuantSettings) -> Vec<u32> {
    let mut chosen: Vec<u32> = Vec::new();
    for p in ions {
        if chosen.len() >= q.ion_number {
            break;
        }
        let m = p.mass as f32;
        if p.quality < q.quality_cutoff
            && !q.excluded.iter().any(|(lo, hi)| *lo <= m && m <= *hi)
            && chosen.iter().all(|c| (*c as f32 - m).abs() >= q.min_dist)
        {
            chosen.push(p.mass);
        }
    }
    // Each fallback step counts, even when it adds nothing.
    let missing = q.ion_number.saturating_sub(chosen.len());
    for p in ions.iter().take(missing) {
        if !chosen.contains(&p.mass) {
            chosen.push(p.mass);
        }
    }
    chosen.sort_unstable();
    chosen
}

/// How `table` groups runs and which substances it keeps.
#[derive(Debug, Clone, PartialEq)]
pub struct TableSettings {
    pub quant: QuantSettings,
    /// Replicate group of a run: the first capture group (else the whole match) of this pattern in
    /// the file name. Runs without a match form a group of their own; an empty pattern puts all runs
    /// in one group.
    pub replicates: String,
    /// Keep a substance only if some group quantifies it in at least this fraction of its runs.
    pub min_reproducibility: f64,
    /// A value counts only if at least this many of the substance's quant ions were found (at most
    /// all of them).
    pub min_ions: usize,
}

impl Default for TableSettings {
    fn default() -> Self {
        Self {
            quant: QuantSettings::default(),
            replicates: String::new(),
            min_reproducibility: 0.0,
            min_ions: 1,
        }
    }
}

/// Substance × run table; columns ordered by replicate group, then input order.
#[derive(Debug, Clone, PartialEq)]
pub struct Table {
    /// File name of each column.
    pub runs: Vec<String>,
    pub groups: Vec<String>,
    /// Ordered by retention time.
    pub rows: Vec<Row>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Row {
    /// Index into the library.
    pub compound: usize,
    pub name: String,
    /// Mean retention time over the runs that found the substance.
    pub rt_s: f64,
    /// Quant ions: the library's own (manual, then automatic), else §7.2 on the best-scoring run.
    pub ions: Vec<u32>,
    /// Summed integral of the quant ions found, per column; `None` if not found or below `min_ions`.
    pub values: Vec<Option<f64>>,
}

impl Table {
    /// RFC 4180 CSV with `.` decimals and `,` separators, readable by R (`read.csv`); in Excel use Data → From Text/CSV (UTF-8, comma):
    /// `Substance,RT_s,Quant_ions,<one column per run>`. Quant ions are space-separated; a value
    /// that was not found is an empty cell, which R reads as `NA`.
    pub fn write_csv(&self, mut w: impl std::io::Write) -> std::io::Result<()> {
        use crate::csv::quote as field;
        let head = ["Substance", "RT_s", "Quant_ions"].map(String::from);
        let head: Vec<String> = head
            .into_iter()
            .chain(self.runs.iter().map(|r| field(r)))
            .collect();
        writeln!(w, "{}", head.join(","))?;
        for r in &self.rows {
            let ions: Vec<String> = r.ions.iter().map(u32::to_string).collect();
            let cells = [field(&r.name), r.rt_s.to_string(), ions.join(" ")]
                .into_iter()
                .chain(
                    r.values
                        .iter()
                        .map(|v| v.map(|v| v.to_string()).unwrap_or_default()),
                );
            writeln!(w, "{}", cells.collect::<Vec<_>>().join(","))?;
        }
        w.flush()
    }
}

pub fn table(
    runs: &[RunResult],
    lib: &[LibraryCompound],
    s: &TableSettings,
) -> Result<Table, regex::Error> {
    let pattern = Regex::new(&s.replicates)?;
    let name = |r: &RunResult| {
        r.path
            .file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .into_owned()
    };
    let group = |n: &str| match pattern.captures(n) {
        _ if s.replicates.is_empty() => String::new(),
        Some(c) => c.get(1).or(c.get(0)).map_or("", |m| m.as_str()).into(),
        None => n.into(),
    };
    let mut cols: Vec<(String, String, &RunResult)> = runs
        .iter()
        .map(|r| {
            let n = name(r);
            (group(&n), n, r)
        })
        .collect();
    cols.sort_by(|a, b| a.0.cmp(&b.0));
    let mut sizes: BTreeMap<&str, usize> = BTreeMap::new();
    cols.iter()
        .for_each(|c| *sizes.entry(&c.0).or_default() += 1);

    // Library compound → its hit in each column.
    let mut hits: BTreeMap<usize, Vec<Option<&Found>>> = BTreeMap::new();
    for (col, (_, _, r)) in cols.iter().enumerate() {
        for f in &r.found {
            hits.entry(f.compound)
                .or_insert_with(|| vec![None; cols.len()])[col] = Some(f);
        }
    }
    let mut rows: Vec<Row> = hits
        .into_iter()
        .filter_map(|(compound, found)| {
            let c = &lib[compound];
            let given = if c.quant_ions.is_empty() {
                &c.auto_quant_ions
            } else {
                &c.quant_ions
            };
            let mut ions: Vec<u32> = given.iter().map(|m| m.round() as u32).collect();
            if ions.is_empty() {
                let best = found
                    .iter()
                    .flatten()
                    .max_by(|a, b| a.score.total_cmp(&b.score))?;
                ions = quant_ions(&best.ions, &s.quant);
            }
            ions.dedup();
            let need = s.min_ions.clamp(1, ions.len().max(1));
            let values: Vec<Option<f64>> = found
                .iter()
                .map(|f| {
                    let f = (*f)?;
                    // The largest peak of each quant ion; a duplicate peak (§4.5) counts once.
                    let hit: Vec<f64> = ions
                        .iter()
                        .filter_map(|m| f.ions.iter().find(|p| p.mass == *m))
                        .map(|p| p.integral)
                        .collect();
                    (hit.len() >= need).then(|| hit.iter().sum())
                })
                .collect();
            let mut per_group: BTreeMap<&str, usize> = BTreeMap::new();
            for (v, col) in values.iter().zip(&cols) {
                if v.is_some() {
                    *per_group.entry(&col.0).or_default() += 1;
                }
            }
            let reproducible = per_group
                .iter()
                .any(|(g, n)| *n as f64 >= s.min_reproducibility * sizes[g] as f64);
            let rts: Vec<f64> = found.iter().flatten().map(|f| f.rt_s).collect();
            reproducible.then(|| Row {
                compound,
                name: c.name.clone(),
                rt_s: rts.iter().sum::<f64>() / rts.len() as f64,
                ions,
                values,
            })
        })
        .collect();
    rows.sort_by(|a, b| a.rt_s.total_cmp(&b.rt_s));
    let (groups, runs) = cols.into_iter().map(|(g, n, _)| (g, n)).unzip();
    Ok(Table { runs, groups, rows })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::library::parse_msl;

    fn ion(mass: u32, integral: f64, quality: f32) -> IonPeak {
        IonPeak {
            mass,
            integral,
            quality,
        }
    }

    #[test]
    fn csv_quotes_fields_and_leaves_missing_values_empty() {
        let t = Table {
            runs: vec!["a,1.cdf".into(), "b.cdf".into()],
            groups: vec![String::new(); 2],
            rows: vec![Row {
                compound: 0,
                name: "Glucose \"1\"".into(),
                rt_s: 12.5,
                ions: vec![73, 147],
                values: vec![Some(1500.0), None],
            }],
        };
        let mut out = Vec::new();
        t.write_csv(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "Substance,RT_s,Quant_ions,\"a,1.cdf\",b.cdf
\"Glucose \"\"1\"\"\",12.5,73 147,1500,
"
        );
    }

    #[test]
    fn csv_reports_a_failing_flush() {
        struct Full;
        impl std::io::Write for Full {
            fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
                Ok(b.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Err(std::io::Error::other("disk full"))
            }
        }
        let t = Table {
            runs: vec![],
            groups: vec![],
            rows: vec![],
        };
        assert!(t.write_csv(Full).is_err());
    }

    #[test]
    fn quant_ions_follow_old_rules() {
        let q = QuantSettings {
            excluded: vec![(73.0, 73.0)],
            ..Default::default()
        };
        // Sorted by integral: 73 is excluded, 147 too poor, so 100, 200 and 300 are picked.
        let ions = [
            ion(73, 90.0, 0.1),
            ion(147, 80.0, 5.0),
            ion(300, 70.0, 0.1),
            ion(100, 60.0, 0.1),
            ion(200, 50.0, 0.1),
        ];
        assert_eq!(quant_ions(&ions, &q), vec![100, 200, 300]);
        // Only 300 passes: the fallback restarts at the top and counts 300 again, so 147 is never
        // reached and fewer than three ions come back.
        let ions = [ion(73, 70.0, 0.1), ion(300, 60.0, 0.1), ion(147, 50.0, 5.0)];
        assert_eq!(quant_ions(&ions, &q), vec![73, 300]);
        let two = QuantSettings { ion_number: 2, ..q };
        assert_eq!(
            quant_ions(&[ion(300, 70.0, 0.1), ion(73, 60.0, 0.1)], &two),
            vec![300]
        );
        // A duplicated peak (§4.5) uses up a fallback step, as in the old program.
        let dup = [
            ion(100, 9.0, 5.0),
            ion(100, 9.0, 5.0),
            ion(200, 8.0, 5.0),
            ion(300, 7.0, 5.0),
        ];
        assert_eq!(quant_ions(&dup, &QuantSettings::default()), vec![100, 200]);
    }

    fn found(compound: usize, rt_s: f64, ions: &[(u32, f64)]) -> Found {
        Found {
            compound,
            score: 0.9,
            rt_s,
            ions: ions.iter().map(|(m, i)| ion(*m, *i, 0.0)).collect(),
        }
    }

    fn run(name: &str, found: Vec<Found>) -> RunResult {
        RunResult {
            path: PathBuf::from("data").join(name),
            found,
        }
    }

    #[test]
    fn table_groups_filters_and_sums() {
        let lib = parse_msl(
            "NAME: A\nQUANTIFICATION: 100 200\nNUM PEAKS: 1\n(100 1)\n\nNAME: B\nNUM PEAKS: 1\n(50 1)\n",
        )
        .unwrap();
        let runs = [
            run(
                "ctrl_1.cdf",
                vec![found(
                    0,
                    10.0,
                    &[(300, 4.0), (200, 2.0), (100, 1.0), (100, 1.0)],
                )],
            ),
            run("mut_1.cdf", vec![found(1, 20.0, &[(50, 5.0)])]),
            run(
                "ctrl_2.cdf",
                vec![found(0, 12.0, &[(100, 1.5)]), found(1, 22.0, &[(50, 7.0)])],
            ),
        ];
        let s = TableSettings {
            replicates: "^([a-z]+)_".into(),
            min_ions: 2,
            ..Default::default()
        };
        let t = table(&runs, &lib, &s).unwrap();
        assert_eq!(t.runs, ["ctrl_1.cdf", "ctrl_2.cdf", "mut_1.cdf"]);
        assert_eq!(t.groups, ["ctrl", "ctrl", "mut"]);
        // A uses the library's ions; ctrl_2 has only one of them, below min_ions.
        assert_eq!(t.rows[0].ions, [100, 200]);
        assert_eq!(t.rows[0].values, [Some(3.0), None, None]);
        assert_eq!(t.rows[0].rt_s, 11.0);
        // B has no library ions: §7.2 on its peaks gives the only one, so min_ions is capped at 1.
        assert_eq!(t.rows[1].ions, [50]);
        assert_eq!(t.rows[1].values, [None, Some(7.0), Some(5.0)]);
        // Reproducibility: A reaches 1/2 in ctrl, B 1/1 in mut.
        let strict = TableSettings {
            min_reproducibility: 0.75,
            ..s.clone()
        };
        let t = table(&runs, &lib, &strict).unwrap();
        assert_eq!(
            t.rows.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
            ["B"]
        );
        // Without a pattern all runs form one group, in input order.
        let one = table(&runs, &lib, &TableSettings::default()).unwrap();
        assert_eq!(one.runs, ["ctrl_1.cdf", "mut_1.cdf", "ctrl_2.cdf"]);
        assert!(
            table(
                &runs,
                &lib,
                &TableSettings {
                    replicates: "(".into(),
                    ..s
                }
            )
            .is_err()
        );
    }

    /// Two substances of 12 ions each, apart in time; the library holds their spectra.
    #[test]
    fn scan_finds_and_integrates_library_substances() {
        let (scans, ions) = (200, 24);
        let c = Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64 * 0.5).collect(),
            offsets: (0..=scans).map(|s| s * ions).collect(),
            mz: (0..scans * ions)
                .map(|k| 60.0 + (k % ions) as f32)
                .collect(),
            intensity: (0..scans * ions)
                .map(|k| {
                    let (s, i) = ((k / ions) as f32, k % ions);
                    let apex = if i < 12 { 60.0 } else { 140.0 };
                    // Baseline 1000 with noise of ±20 (like the peak detector tests).
                    let noise = ((k as u32).wrapping_mul(2_654_435_761) >> 16) % 41;
                    let peak =
                        5000.0 * (i % 12 + 1) as f32 * (-0.5 * ((s - apex) / 2.5).powi(2)).exp();
                    980.0 + noise as f32 + peak
                })
                .collect(),
        };
        let spectrum = |first: usize| -> String {
            (0..12)
                .map(|j| format!("({} {}) ", 60 + first + j, j + 1))
                .collect()
        };
        let msl = format!(
            "NAME: late\nNUM PEAKS: 12\n{}\n\nNAME: early\nNUM PEAKS: 12\n{}\n",
            spectrum(12),
            spectrum(0)
        );
        let lib = parse_msl(&msl).unwrap();
        let s = ScanSettings::default();
        let found = scan_chromatogram(c, &PreparedLibrary::new(&lib, &s.search), &s).unwrap();
        assert_eq!(found.len(), 2, "{found:?}");
        let early = found.iter().find(|f| f.compound == 1).unwrap();
        assert!((early.rt_s - 30.0).abs() <= 0.5, "{}", early.rt_s);
        assert_eq!(early.ions.len(), 12);
        assert_eq!(early.ions[0].mass, 71, "largest integral first");
        assert!(
            early
                .ions
                .windows(2)
                .all(|w| w[0].integral >= w[1].integral)
        );
    }

    #[test]
    fn scan_runs_reports_errors_and_cancels() {
        let lib = PreparedLibrary::new(&[], &SearchSettings::default());
        let paths = [PathBuf::from("no such run.cdf")];
        let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
        let s = ScanSettings::default();
        let (runs, errors) = scan_runs(&paths, &lib, &s, 2, &done, &cancel).unwrap();
        assert!(runs.is_empty() && errors.len() == 1 && errors[0].contains("no such run"));
        assert_eq!(done.load(Relaxed), 1);
        cancel.store(true, Relaxed);
        assert!(scan_runs(&paths, &lib, &s, 2, &done, &cancel).is_none());
    }
}
