//! Batch quantification: finds the library substances in many runs and builds a substance × run
//! table (old `BatchQuantificator`, feature-inventory.md F10–F12). Quant ions follow algorithms.md §7.
//!
//! Passes: `scan_runs` (expensive, once) keeps per run the identified substances and their ion
//! integrals; the optional untargeted pass (`untargeted`) pools the unidentified components of all runs
//! into a temporary library and `rescan` looks for quant ions a run lacks in the ion traces; `table`
//! (cheap) turns the result into the table, so filters and replicate groups can change without rescanning.
//! Details of the old program: docs/reverse-engineering/batch-quantification.md.

mod detail;
mod edit;
mod rescan;
mod untargeted;

pub use detail::{cell_ions, hit_list, ion_integrals, ions_path, settings_path};
pub use rescan::rescan_runs;
pub use untargeted::add_unknowns;

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};

use regex::Regex;

use crate::Chromatogram;
use crate::deconv::{Component, DeconvSettings, deconvolute};
use crate::library::LibraryCompound;
use crate::peaks::{IonTraces, Peak, PeakSettings, TooLarge};
use crate::ri::RiCalibration;
use crate::search::{
    PreparedLibrary, SearchSettings, map_compounds_to_library, shift_to_reference,
};

/// Settings of the per-run pass: the single-run defaults, one calibration for every run.
#[derive(Debug, Clone, Default)]
pub struct ScanSettings {
    pub peaks: PeakSettings,
    pub deconv: DeconvSettings,
    pub search: SearchSettings,
    /// RT → RI for every run; without it the search compares spectra only.
    pub calibration: Option<RiCalibration>,
    /// Shifts `calibration` per run onto a reference substance (RT shift correction, §6.5).
    pub shift: Option<ShiftReference>,
    /// Keep the components without a library hit (`Scanned::unknown`) for `add_unknowns`.
    pub untargeted: bool,
}

/// RT shift reference: library index and the minimum overall score it must reach in a run.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ShiftReference {
    pub compound: usize,
    pub cutoff: f64,
}

impl ShiftReference {
    /// The library entry called `name`; only a substance with an RI can be the reference.
    pub fn find(lib: &[LibraryCompound], name: &str, cutoff: f64) -> Option<Self> {
        let compound = lib
            .iter()
            .position(|c| c.name == name && c.spectrum.ri >= 0.0)?;
        Some(Self { compound, cutoff })
    }
}

/// One ion peak of a found substance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IonPeak {
    pub mass: u32,
    /// Integral from the deconvolution (`FittedPeak::integral`).
    pub integral: f64,
    /// Peak quality, lower is better (§3.7).
    pub quality: f32,
    /// Signal-to-noise ratio of the ion peak (`Peak::signal_to_noise`).
    pub signal_to_noise: f32,
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
    /// Made by `rescan_runs` from ion traces only: the substance was not deconvoluted in this run.
    pub rescanned: bool,
    /// Masses of the ions `rescan_runs` added to an identified component.
    pub rescued: Vec<u32>,
    /// The other library hits of the component (library index, overall score), best first; the
    /// assigned hit is `compound`/`score`. Empty after a rescan.
    pub alternatives: Vec<(usize, f64)>,
}

/// A component without library hit, kept by an untargeted scan until `add_unknowns` pools it.
#[derive(Debug, Clone, PartialEq)]
pub struct Unknown {
    pub rt_s: f64,
    /// RI at the apex under the run's calibration.
    pub ri: Option<f64>,
    /// (nominal mass, height), ascending mass.
    pub spectrum: Vec<(u32, u32)>,
    pub ions: Vec<IonPeak>,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct RunResult {
    pub path: PathBuf,
    pub found: Vec<Found>,
    /// Unidentified components; empty after `add_unknowns`.
    pub unknown: Vec<Unknown>,
    /// The calibration of this run, shifted onto the reference if one was found (`None`: not calibrated).
    pub calibration: Option<RiCalibration>,
}

/// What one run yields.
#[derive(Debug, Clone, PartialEq)]
pub struct Scanned {
    pub found: Vec<Found>,
    pub unknown: Vec<Unknown>,
    pub calibration: Option<RiCalibration>,
}

/// The ion peaks of a component, largest integral first (§7.2).
pub fn ion_peaks(k: &Component, peaks: &[Peak], traces: &IonTraces) -> Vec<IonPeak> {
    ion_peaks_by(k, peaks, |i| traces.masses[peaks[i].row as usize])
}

/// [`ion_peaks`] with the ion mass of peak `i` given by `mass_of` (peaks read from a `.cmp` have no
/// trace row).
pub fn ion_peaks_by(k: &Component, peaks: &[Peak], mass_of: impl Fn(usize) -> u32) -> Vec<IonPeak> {
    let mut ions: Vec<IonPeak> = k
        .peaks
        .iter()
        .map(|f| {
            let p = &peaks[f.peak];
            IonPeak {
                mass: mass_of(f.peak),
                integral: f.integral,
                quality: p.quality,
                signal_to_noise: p.signal_to_noise,
            }
        })
        .collect();
    // The old code sorts ascending and walks backwards, so ties come last-in-list first (§7.2).
    ions.sort_by(|a, b| a.integral.total_cmp(&b.integral));
    ions.reverse();
    ions
}

/// Deconvolutes one run and returns the library substances assigned to its components (and, if
/// `s.untargeted`, the other components).
pub fn scan_chromatogram(
    c: Chromatogram,
    lib: &PreparedLibrary,
    s: &ScanSettings,
) -> Result<Scanned, TooLarge> {
    let mut traces = IonTraces::new(&c)?;
    // The ion-major copy and the times are all the rest needs; dropping the spectra now halves the
    // peak memory per run.
    let times = c.scan_times_s.clone();
    drop(c);
    let peaks = traces.detect_peaks(&s.peaks);
    let comps = deconvolute(&traces, &peaks, &s.deconv);
    Ok(scan_components(&times, &comps, &peaks, lib, s, |i| {
        traces.masses[peaks[i].row as usize]
    }))
}

/// The run behind a `.cmp` input: its saved components instead of a deconvolution. The scan times
/// are right at the apex scans only, which is all the later passes read.
fn scan_saved(path: &Path, lib: &PreparedLibrary, s: &ScanSettings) -> Result<Scanned, String> {
    let ((comps, masses, peaks), times) =
        crate::cmp::standalone_components(path).map_err(|e| e.to_string())?;
    Ok(scan_components(&times, &comps, &peaks, lib, s, |i| {
        masses[i]
    }))
}

/// Calibration, library search and ion peaks of the components of one run.
fn scan_components(
    times: &[f64],
    comps: &[Component],
    peaks: &[Peak],
    lib: &PreparedLibrary,
    s: &ScanSettings,
    mass_of: impl Fn(usize) -> u32,
) -> Scanned {
    let calibration = s.calibration.as_ref().map(|cal| {
        let shifted = s.shift.and_then(|r| {
            shift_to_reference(cal, times, comps, lib, r.compound, r.cutoff, &s.search)
        });
        shifted.unwrap_or_else(|| cal.clone())
    });
    let ris = calibration
        .as_ref()
        .map(|cal| cal.component_ris(times, comps));
    let id = map_compounds_to_library(comps, ris.as_deref(), lib, &s.search);
    let (mut found, mut unknown) = (Vec::new(), Vec::new());
    for (i, (k, hit)) in comps.iter().zip(&id.assigned).enumerate() {
        if hit.is_none() && !s.untargeted {
            continue;
        }
        let ions = ion_peaks_by(k, peaks, &mass_of);
        let rt_s = times[k.apex_scan as usize];
        match hit {
            Some(h) => found.push(Found {
                compound: h.compound,
                score: h.overall,
                rt_s,
                ions,
                rescanned: false,
                rescued: Vec::new(),
                alternatives: (id.hits[i].iter())
                    .filter(|o| o.compound != h.compound)
                    .map(|o| (o.compound, o.overall))
                    .collect(),
            }),
            None => unknown.push(Unknown {
                rt_s,
                ri: ris.as_ref().map(|r| r[i]),
                spectrum: k.spectrum.clone(),
                ions,
            }),
        }
    }
    Scanned {
        found,
        unknown,
        calibration,
    }
}

/// Scans every file with at most `parallel` runs in memory at once; each run is analysed in parallel
/// (rayon) as well. A `.cmp` file stands for a run already analysed: its saved components are searched
/// as they are. `done` counts finished files. Once `cancel` is set, no further file starts and
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
                            match crate::cmp::is_cmp(path) {
                                true => scan_saved(path, lib, s),
                                false => crate::read_run(path).and_then(|c| {
                                    scan_chromatogram(c, lib, s).map_err(|e| e.to_string())
                                }),
                            }
                        }))
                        .unwrap_or_else(|_| Err("analysis crashed".into()))
                        .map(|r| RunResult {
                            path: path.clone(),
                            found: r.found,
                            unknown: r.unknown,
                            calibration: r.calibration,
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

/// What `complete_scan` added to the scanned runs.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct Completed {
    /// The pooled unknown substances (untargeted); `table` takes them as `extra`.
    pub extra: Vec<LibraryCompound>,
    /// Ions the rescan added.
    pub rescanned: usize,
    /// Run files the rescan could not read.
    pub errors: Vec<String>,
}

/// The passes after `scan_runs`, shared by the app and the command line: pools the unknown
/// components if `s.untargeted`, then rescans missing quant ions (chosen by `q`) if `rescan` is set.
/// `None` once `cancel` is set.
#[allow(clippy::too_many_arguments)]
pub fn complete_scan(
    runs: &mut [RunResult],
    lib: &[LibraryCompound],
    s: &ScanSettings,
    q: &QuantSettings,
    rescan: bool,
    parallel: usize,
    cancel: &AtomicBool,
) -> Option<Completed> {
    let extra = if s.untargeted {
        add_unknowns(runs, lib.len(), s, cancel)?
    } else {
        Vec::new()
    };
    let (rescanned, errors) = if rescan {
        rescan_runs(runs, lib, &extra, s, q, parallel, cancel)?
    } else {
        (0, Vec::new())
    };
    Some(Completed {
        extra,
        rescanned,
        errors,
    })
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
    /// An ion counts only if its peak has at least this signal-to-noise ratio (old "Req. S/N").
    pub min_sn: f32,
    /// A value counts only if the mean peak quality of its ions (the discrepancy index) is at most
    /// this; the default accepts everything.
    pub max_discrepancy: f32,
}

impl Default for TableSettings {
    fn default() -> Self {
        Self {
            quant: QuantSettings::default(),
            replicates: String::new(),
            min_reproducibility: 0.0,
            min_ions: 1,
            min_sn: 0.0,
            max_discrepancy: f32::INFINITY,
        }
    }
}

/// Substance × run table; columns ordered by replicate group, then input order.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Table {
    /// File name of each column.
    pub runs: Vec<String>,
    /// Index of each column's run in the list `table()` was called with: file names alone do not
    /// tell `dir1/a.cdf` from `dir2/a.cdf`.
    pub sources: Vec<usize>,
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
    /// Retention time of the hit in each column, where there is a value.
    pub rts_s: Vec<Option<f64>>,
}

/// UTF-8 byte order mark: Excel needs it to read umlauts in a CSV (R: `fileEncoding = "UTF-8-BOM"`).
const BOM: &str = "\u{feff}";

/// A text cell: RFC 4180 quoting, and a leading `'` where a spreadsheet would read a formula
/// (`=`, `+`, `-`, `@`, tab or CR first); `scripts/r/load_results.R` takes that `'` off again.
fn text_cell(s: &str) -> String {
    let formula = s.starts_with(['=', '+', '-', '@', '\t', '\r']);
    crate::csv::quote(&if formula {
        format!("'{s}")
    } else {
        s.to_owned()
    })
}

/// Where the replicate groups of a results CSV go: `results.csv` → `results.groups.csv`.
pub fn groups_path(csv: &Path) -> PathBuf {
    csv.with_extension("groups.csv")
}

impl Table {
    /// UTF-8 CSV with a byte order mark, `.` decimals and `,` separators (RFC 4180), readable by R
    /// (`read.csv`) and Excel: `Substance,RT_s,Quant_ions,<one column per run>`. Quant ions are
    /// space-separated; a value that was not found is an empty cell, which R reads as `NA`. Names
    /// that a spreadsheet would read as a formula get a leading `'`.
    pub fn write_csv(&self, mut w: impl std::io::Write) -> std::io::Result<()> {
        use text_cell as field;
        w.write_all(BOM.as_bytes())?;
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

    /// Whether any run belongs to a replicate group.
    pub fn has_groups(&self) -> bool {
        self.groups.iter().any(|g| !g.is_empty())
    }

    /// Writes the groups next to the results CSV `csv` (`groups_path`) and returns that file; if no
    /// run has a group, an old groups file of the same name is removed (a script would read it) and
    /// there is none to return.
    pub fn write_groups_file(&self, csv: &Path) -> std::io::Result<Option<PathBuf>> {
        let path = groups_path(csv);
        if !self.has_groups() {
            return match std::fs::remove_file(&path) {
                Err(e) if e.kind() != std::io::ErrorKind::NotFound => Err(e),
                _ => Ok(None),
            };
        }
        self.write_groups_csv(std::io::BufWriter::new(std::fs::File::create(&path)?))?;
        Ok(Some(path))
    }

    /// The replicate group of each run, as `Run,Group` (same encoding as `write_csv`; a run without
    /// group has an empty cell), so a script can join it to the columns of the results CSV.
    pub fn write_groups_csv(&self, mut w: impl std::io::Write) -> std::io::Result<()> {
        w.write_all(BOM.as_bytes())?;
        writeln!(w, "Run,Group")?;
        for (run, group) in self.runs.iter().zip(&self.groups) {
            writeln!(w, "{},{}", text_cell(run), text_cell(group))?;
        }
        w.flush()
    }
}

/// Library compound `i`; indices past the library address `extra` (the pooled unknowns).
fn compound_of<'a>(
    lib: &'a [LibraryCompound],
    extra: &'a [LibraryCompound],
    i: usize,
) -> &'a LibraryCompound {
    lib.get(i).unwrap_or_else(|| &extra[i - lib.len()])
}

/// Library compound → its hit in each run.
fn hits_by_compound<'a>(runs: &[&'a RunResult]) -> BTreeMap<usize, Vec<Option<&'a Found>>> {
    let mut hits: BTreeMap<usize, Vec<Option<&Found>>> = BTreeMap::new();
    for (col, r) in runs.iter().enumerate() {
        for f in &r.found {
            hits.entry(f.compound)
                .or_insert_with(|| vec![None; runs.len()])[col] = Some(f);
        }
    }
    hits
}

/// Quant ions of a substance: the library's own (manual, then automatic), else §7.2 on the
/// best-scoring run among `found`. `None` if there is no run at all.
fn compound_ions(
    c: &LibraryCompound,
    found: &[Option<&Found>],
    q: &QuantSettings,
) -> Option<Vec<u32>> {
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
            .filter(|f| !f.rescanned)
            .max_by(|a, b| a.score.total_cmp(&b.score))?;
        ions = quant_ions(&best.ions, q);
    }
    ions.dedup();
    Some(ions)
}

/// The cell of one run: the summed integral of the quant ions found, if enough of them pass the filters.
fn cell(f: &Found, ions: &[u32], need: usize, s: &TableSettings) -> Option<f64> {
    let hit: Vec<&IonPeak> = detail::quant_peaks(f, ions, s.min_sn)
        .into_iter()
        .flatten()
        .collect();
    let discrepancy = hit.iter().map(|p| p.quality).sum::<f32>() / hit.len().max(1) as f32;
    (hit.len() >= need && discrepancy <= s.max_discrepancy)
        .then(|| hit.iter().map(|p| p.integral).sum())
}

/// `extra` are the unknowns `add_unknowns` pooled; their indices follow the library's.
pub fn table(
    runs: &[RunResult],
    lib: &[LibraryCompound],
    extra: &[LibraryCompound],
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
    let mut cols: Vec<(String, String, usize, &RunResult)> = runs
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let n = name(r);
            (group(&n), n, i, r)
        })
        .collect();
    cols.sort_by(|a, b| a.0.cmp(&b.0));
    let mut sizes: BTreeMap<&str, usize> = BTreeMap::new();
    cols.iter()
        .for_each(|c| *sizes.entry(&c.0).or_default() += 1);

    let ordered: Vec<&RunResult> = cols.iter().map(|c| c.3).collect();
    let mut rows: Vec<Row> = hits_by_compound(&ordered)
        .into_iter()
        .filter_map(|(compound, found)| {
            let c = compound_of(lib, extra, compound);
            let ions = compound_ions(c, &found, &s.quant)?;
            let need = s.min_ions.clamp(1, ions.len().max(1));
            let values: Vec<Option<f64>> = found
                .iter()
                .map(|f| cell(f.as_ref()?, &ions, need, s))
                .collect();
            let rts_s = (values.iter().zip(&found))
                .map(|(v, f)| v.and(f.map(|f| f.rt_s)))
                .collect();
            let mut per_group: BTreeMap<&str, usize> = BTreeMap::new();
            // The old reproducibility filter runs before the rescan: what the rescan added does not count.
            let counts = |v: &Option<f64>, f: &Option<&Found>| match f {
                Some(f) if !f.rescanned && f.rescued.is_empty() => v.is_some(),
                Some(f) if !f.rescanned => {
                    let mut plain = (*f).clone();
                    plain.ions.retain(|p| !f.rescued.contains(&p.mass));
                    cell(&plain, &ions, need, s).is_some()
                }
                _ => false,
            };
            for ((v, col), f) in values.iter().zip(&cols).zip(&found) {
                if counts(v, f) {
                    *per_group.entry(&col.0).or_default() += 1;
                }
            }
            let reproducible = per_group
                .iter()
                .any(|(g, n)| *n as f64 >= s.min_reproducibility * sizes[g] as f64);
            // Substances only rescanned in some runs have their retention time from the others.
            let rts: Vec<f64> = found
                .iter()
                .flatten()
                .filter(|f| !f.rescanned)
                .map(|f| f.rt_s)
                .collect();
            reproducible.then(|| Row {
                compound,
                name: c.name.clone(),
                rt_s: rts.iter().sum::<f64>() / rts.len() as f64,
                ions,
                values,
                rts_s,
            })
        })
        .collect();
    rows.sort_by(|a, b| a.rt_s.total_cmp(&b.rt_s));
    let sources = cols.iter().map(|c| c.2).collect();
    let (groups, runs) = cols.into_iter().map(|(g, n, ..)| (g, n)).unzip();
    Ok(Table {
        runs,
        sources,
        groups,
        rows,
    })
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
            signal_to_noise: 10.0,
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
                rts_s: vec![Some(750.0), None],
            }],
            ..Default::default()
        };
        let mut out = Vec::new();
        t.write_csv(&mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "\u{feff}Substance,RT_s,Quant_ions,\"a,1.cdf\",b.cdf
\"Glucose \"\"1\"\"\",12.5,73 147,1500,
"
        );
    }

    #[test]
    fn csv_protects_formulas_and_lists_the_groups() {
        let t = Table {
            runs: vec!["=cmd.cdf".into(), "ctrl_1.cdf".into(), "x.cdf".into()],
            groups: vec!["+a".into(), "ctrl".into(), String::new()],
            rows: vec![Row {
                compound: 0,
                name: "-Pinene".into(),
                rt_s: 1.0,
                ions: vec![93],
                values: vec![Some(-2.0), None, None],
                rts_s: vec![None; 3],
            }],
            ..Default::default()
        };
        let (mut csv, mut groups) = (Vec::new(), Vec::new());
        t.write_csv(&mut csv).unwrap();
        t.write_groups_csv(&mut groups).unwrap();
        // Text is guarded; a negative number is still a number.
        assert_eq!(
            String::from_utf8(csv).unwrap(),
            "\u{feff}Substance,RT_s,Quant_ions,'=cmd.cdf,ctrl_1.cdf,x.cdf\n'-Pinene,1,93,-2,,\n"
        );
        assert_eq!(
            String::from_utf8(groups).unwrap(),
            "\u{feff}Run,Group\n'=cmd.cdf,'+a\nctrl_1.cdf,ctrl\nx.cdf,\n"
        );
        assert!(t.has_groups());
        assert_eq!(
            groups_path(Path::new("dir/results.csv")),
            Path::new("dir/results.groups.csv")
        );
    }

    #[test]
    fn an_old_groups_file_is_removed_when_there_are_no_groups() {
        let dir = std::env::temp_dir().join(format!("md-groups-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let csv = dir.join("results.csv");
        let mut t = Table {
            runs: vec!["a.cdf".into()],
            groups: vec!["ctrl".into()],
            rows: vec![],
            ..Default::default()
        };
        let written = t.write_groups_file(&csv).unwrap().unwrap();
        assert!(written.exists());
        t.groups = vec![String::new()];
        assert_eq!(t.write_groups_file(&csv).unwrap(), None);
        assert!(!written.exists());
        // Nothing to remove is not an error.
        assert_eq!(t.write_groups_file(&csv).unwrap(), None);
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
            ..Default::default()
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
            rescanned: false,
            rescued: Vec::new(),
            alternatives: Vec::new(),
        }
    }

    fn run(name: &str, found: Vec<Found>) -> RunResult {
        RunResult {
            path: PathBuf::from("data").join(name),
            found,
            ..Default::default()
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
        let t = table(&runs, &lib, &[], &s).unwrap();
        assert_eq!(t.runs, ["ctrl_1.cdf", "ctrl_2.cdf", "mut_1.cdf"]);
        assert_eq!(t.groups, ["ctrl", "ctrl", "mut"]);
        // A uses the library's ions; ctrl_2 has only one of them, below min_ions.
        assert_eq!(t.rows[0].ions, [100, 200]);
        assert_eq!(t.rows[0].values, [Some(3.0), None, None]);
        assert_eq!(t.rows[0].rt_s, 11.0);
        // B has no library ions: §7.2 on its peaks gives the only one, so min_ions is capped at 1.
        assert_eq!(t.rows[1].ions, [50]);
        assert_eq!(t.rows[1].values, [None, Some(7.0), Some(5.0)]);
        // The time of a hit is kept only where the cell has a value.
        assert_eq!(t.rows[0].rts_s, [Some(10.0), None, None]);
        assert_eq!(t.rows[1].rts_s, [None, Some(22.0), Some(20.0)]);
        // Reproducibility: A reaches 1/2 in ctrl, B 1/1 in mut.
        let strict = TableSettings {
            min_reproducibility: 0.75,
            ..s.clone()
        };
        let t = table(&runs, &lib, &[], &strict).unwrap();
        assert_eq!(
            t.rows.iter().map(|r| r.name.as_str()).collect::<Vec<_>>(),
            ["B"]
        );
        // Without a pattern all runs form one group, in input order.
        let one = table(&runs, &lib, &[], &TableSettings::default()).unwrap();
        assert_eq!(one.runs, ["ctrl_1.cdf", "mut_1.cdf", "ctrl_2.cdf"]);
        assert!(
            table(
                &runs,
                &lib,
                &[],
                &TableSettings {
                    replicates: "(".into(),
                    ..s
                }
            )
            .is_err()
        );
    }

    /// Two substances of 12 ions each, apart in time (early at 30 s, late at 70 s).
    fn two_substances() -> Chromatogram {
        let (scans, ions) = (200, 24);
        Chromatogram {
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
        }
    }

    /// The library of `two_substances`: "late" (index 0) and "early" (index 1).
    fn two_substances_library() -> Vec<LibraryCompound> {
        let spectrum = |first: usize| -> String {
            (0..12)
                .map(|j| format!("({} {}) ", 60 + first + j, j + 1))
                .collect()
        };
        let msl = format!(
            "NAME: late
NUM PEAKS: 12
{}

NAME: early
NUM PEAKS: 12
{}
",
            spectrum(12),
            spectrum(0)
        );
        parse_msl(&msl).unwrap()
    }

    #[test]
    fn scan_finds_and_integrates_library_substances() {
        let lib = two_substances_library();
        let s = ScanSettings::default();
        let found = scan_chromatogram(two_substances(), &PreparedLibrary::new(&lib, &s.search), &s)
            .unwrap()
            .found;
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
        assert!(early.ions.iter().all(|p| p.signal_to_noise > 0.0));
    }

    #[test]
    fn table_filters_by_signal_to_noise_and_discrepancy() {
        let lib = parse_msl(
            "NAME: A
QUANTIFICATION: 100 200
NUM PEAKS: 1
(100 1)
",
        )
        .unwrap();
        let mut good = found(0, 10.0, &[(100, 4.0), (200, 2.0)]);
        good.ions[1].signal_to_noise = 2.0;
        good.ions[0].quality = 1.0;
        good.ions[1].quality = 3.0;
        let runs = [run("a.cdf", vec![good])];
        // A substance without any value left drops out of the table.
        let value = |s: TableSettings| {
            (table(&runs, &lib, &[], &s).unwrap().rows.first()).and_then(|r| r.values[0])
        };
        assert_eq!(value(TableSettings::default()), Some(6.0));
        // S/N 5 drops ion 200; one of two ions is below min_ions 2.
        let sn = TableSettings {
            min_sn: 5.0,
            ..Default::default()
        };
        assert_eq!(value(sn.clone()), Some(4.0));
        assert_eq!(value(TableSettings { min_ions: 2, ..sn }), None);
        // Discrepancy = mean quality of the ions: (1 + 3) / 2.
        let d = |max| TableSettings {
            max_discrepancy: max,
            ..Default::default()
        };
        assert_eq!(value(d(2.0)), Some(6.0));
        assert_eq!(value(d(1.9)), None);
    }

    #[test]
    fn untargeted_runs_pool_their_unknown_components() {
        let s = ScanSettings {
            untargeted: true,
            ..Default::default()
        };
        let none = PreparedLibrary::new(&[], &s.search);
        let one = scan_chromatogram(two_substances(), &none, &s).unwrap();
        assert!(one.found.is_empty());
        assert_eq!(one.unknown.len(), 2);
        let mut runs: Vec<RunResult> = (0..2)
            .map(|i| RunResult {
                path: PathBuf::from(format!("r{i}.cdf")),
                unknown: one.unknown.clone(),
                ..Default::default()
            })
            .collect();
        // A cancelled pooling leaves the components of the runs it did not reach.
        assert!(add_unknowns(&mut runs.clone(), 0, &s, &AtomicBool::new(true)).is_none());
        // Two groups, each with one component of each run.
        let extra = add_unknowns(&mut runs, 0, &s, &AtomicBool::new(false)).unwrap();
        assert_eq!(extra.len(), 2);
        assert!(
            runs.iter()
                .all(|r| r.unknown.is_empty() && r.found.len() == 2)
        );
        let t = table(&runs, &[], &extra, &TableSettings::default()).unwrap();
        assert_eq!(t.rows.len(), 2);
        assert!(t.rows[0].rt_s < t.rows[1].rt_s);
        assert!(t.rows[0].name.starts_with("Unknown"));
        assert!(t.rows.iter().all(|r| r.values.iter().all(Option::is_some)));
    }

    /// Two runs of the same file in a temporary folder: the first identified, the second only scanned.
    fn rescan_fixture(tag: &str) -> (PathBuf, Vec<RunResult>, Vec<LibraryCompound>) {
        let dir = std::env::temp_dir().join(format!("md-rescan-{tag}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (a, b) = (dir.join("a.cdf"), dir.join("b.cdf"));
        crate::andi::write(&a, &two_substances()).unwrap();
        crate::andi::write(&b, &two_substances()).unwrap();
        let lib = two_substances_library();
        let s = ScanSettings::default();
        let found = scan_chromatogram(two_substances(), &PreparedLibrary::new(&lib, &s.search), &s)
            .unwrap()
            .found;
        let run = |path, found| RunResult {
            path,
            found,
            ..Default::default()
        };
        (dir, vec![run(a, found), run(b, Vec::new())], lib)
    }

    #[test]
    fn rescan_fills_the_ions_a_run_lacks() {
        let (dir, mut runs, lib) = rescan_fixture("fill");
        let (s, q) = (ScanSettings::default(), QuantSettings::default());
        let cancel = AtomicBool::new(false);
        let (added, errors) = rescan_runs(&mut runs, &lib, &[], &s, &q, 2, &cancel).unwrap();
        assert!(errors.is_empty(), "{errors:?}");
        assert_eq!(runs[1].found.len(), 2, "both substances are completed");
        let t = table(&runs, &lib, &[], &TableSettings::default()).unwrap();
        assert!(added >= 2 * 3, "{added}");
        for row in &t.rows {
            let (v0, v1) = (row.values[0].unwrap(), row.values[1].unwrap());
            assert!((v1 / v0 - 1.0).abs() < 0.25, "{} {v0} {v1}", row.name);
        }
        // Rescanned entries do not move the retention time of the substance.
        assert!(runs[1].found.iter().all(|f| f.rescanned && f.score == 0.0));
        // ... and do not count for reproducibility, which the old program checks before the rescan.
        let all = TableSettings {
            min_reproducibility: 1.0,
            ..Default::default()
        };
        assert!(table(&runs, &lib, &[], &all).unwrap().rows.is_empty());
        cancel.store(true, Relaxed);
        assert!(rescan_runs(&mut runs, &lib, &[], &s, &q, 1, &cancel).is_none());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn rescan_needs_a_quarter_of_the_runs_and_does_not_count_for_reproducibility() {
        let (dir, runs, lib) = rescan_fixture("limit");
        let (s, q) = (ScanSettings::default(), QuantSettings::default());
        let cancel = AtomicBool::new(false);
        // A substance seen in fewer than a quarter of the runs (rounded: 2 of 6) is left alone.
        let mut again = vec![runs[0].clone()];
        again.extend((0..5).map(|_| runs[1].clone()));
        let (added, _) = rescan_runs(&mut again, &lib, &[], &s, &q, 1, &cancel).unwrap();
        assert_eq!(added, 0);
        // At the limit (2 of 6) it is rescanned.
        again[1] = runs[0].clone();
        let (added, _) = rescan_runs(&mut again, &lib, &[], &s, &q, 1, &cancel).unwrap();
        assert!(added > 0);
        // Ions added to an identified component do not count for reproducibility either.
        let mut partial = vec![runs[0].clone(), runs[0].clone()];
        partial[1].path = runs[1].path.clone();
        // The quant ions come from the better run, which is the first.
        for f in &mut partial[1].found {
            f.ions.clear();
            f.score = 0.5;
        }
        let all = TableSettings {
            min_reproducibility: 1.0,
            ..Default::default()
        };
        let rows = |runs: &[RunResult], s: &TableSettings| table(runs, &lib, &[], s).unwrap().rows;
        assert!(rows(&partial, &all).is_empty());
        rescan_runs(&mut partial, &lib, &[], &s, &q, 1, &cancel).unwrap();
        assert!(partial[1].found.iter().all(|f| !f.ions.is_empty()));
        assert!(rows(&partial, &all).is_empty());
        assert!(!rows(&partial, &TableSettings::default()).is_empty());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_saved_cmp_scans_like_the_run_it_came_from() {
        let dir = std::env::temp_dir().join(format!("md-scan-cmp-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (run, lib) = (two_substances(), two_substances_library());
        let (s, cancel, done) = (
            ScanSettings::default(),
            AtomicBool::new(false),
            AtomicUsize::new(0),
        );
        let prepared = PreparedLibrary::new(&lib, &s.search);
        let direct = scan_chromatogram(run.clone(), &prepared, &s).unwrap();
        assert_eq!(direct.found.len(), 2);

        let mut traces = IonTraces::new(&run).unwrap();
        let peaks = traces.detect_peaks(&s.peaks);
        let comps = deconvolute(&traces, &peaks, &s.deconv);
        let masses: Vec<u32> = peaks
            .iter()
            .map(|p| traces.masses[p.row as usize])
            .collect();
        let compounds =
            crate::cmp::compounds_of_components(&comps, &masses, &peaks, &run.scan_times_s, &[]);
        let path = dir.join("run.cmp");
        crate::cmp::write_cmp(&path, &compounds).unwrap();

        let (runs, errors) = scan_runs(
            std::slice::from_ref(&path),
            &prepared,
            &s,
            1,
            &done,
            &cancel,
        )
        .unwrap();
        assert!(errors.is_empty(), "{errors:?}");
        let saved = &runs[0].found;
        assert_eq!(saved.len(), direct.found.len());
        for (a, b) in saved.iter().zip(&direct.found) {
            assert_eq!(
                (a.compound, a.score, &a.ions),
                (b.compound, b.score, &b.ions)
            );
            assert!((a.rt_s - b.rt_s).abs() < 1e-9 && !a.rescanned);
        }
        // The rescan has no raw data to look into and leaves a `.cmp` run alone.
        let mut runs = runs;
        runs[0].found.clear();
        let (added, errors) = rescan_runs(
            &mut runs,
            &lib,
            &[],
            &s,
            &QuantSettings::default(),
            1,
            &cancel,
        )
        .unwrap();
        assert_eq!((added, errors.len()), (0, 0));
        std::fs::remove_dir_all(&dir).unwrap();
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
