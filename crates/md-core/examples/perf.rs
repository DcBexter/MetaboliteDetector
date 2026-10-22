//! Performance smoke test on a synthetic large run.
//! `cargo run --release -p md-core --example perf -- [scans] [ions_per_scan]`
//! Default: 60 min at 20 scans/s with 500 ions per scan (36 000 scans, 18 M points, 300 000 peaks,
//! 12 000 substances).

use std::path::PathBuf;
use std::time::Instant;

use md_core::andi;
use md_core::batch::{
    Found, IonPeak, RunResult, ScanSettings, Study, Unknown, add_unknowns, run_study,
};
use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::isotope::IsotopeTable;
use md_core::job::{Cancel, Progress};
use md_core::library::{LibraryCompound, parse_msl};
use md_core::mid::{Method, determine, frag_export};
use md_core::peaks::{Detector, IonTraces, PeakSettings, reorder_peaks};
use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};
use md_core::{Chromatogram, decimate};

fn main() {
    let mut args = std::env::args()
        .skip(1)
        .map(|a| a.parse::<usize>().expect("number"));
    let scans = args.next().unwrap_or(36_000);
    let ions = args.next().unwrap_or(500);

    let t = Instant::now();
    let c = Chromatogram {
        scan_times_s: (0..scans).map(|s| s as f64 / 20.0).collect(),
        offsets: (0..=scans).map(|s| s * ions).collect(),
        mz: (0..scans * ions)
            .map(|k| 50.0 + (k % ions) as f32)
            .collect(),
        // Noise plus a Gaussian peak (sigma 2 scans) every 60 scans in each ion trace; every 20th
        // ion shares its apex, so each 60 scans hold 20 co-eluting substances of 25 ions.
        intensity: (0..scans * ions)
            .map(|k| {
                let x = ((k / ions + 7 * (k % ions % 20)) % 60) as f32 - 30.0;
                (k % 997) as f32 * 0.05 + 20_000.0 * (-x * x / 8.0).exp()
            })
            .collect(),
        exact_masses: false,
        baseline: Default::default(),
    };
    let points = c.mz.len();
    println!(
        "{scans} scans x {ions} ions = {points} points, {} MiB, built in {:?}",
        points * 8 / (1 << 20),
        t.elapsed()
    );

    let t = Instant::now();
    let tic = c.tic();
    println!("TIC: {:?}", t.elapsed());

    let t = Instant::now();
    let sic = c.sic(73.0, 0.5, &Cancel::default()).unwrap();
    println!("SIC m/z 73: {:?}", t.elapsed());
    // The first call also pays for touching the 137 MiB; the median of repeats shows the per-call cost.
    let mut repeats: Vec<_> = (0..15)
        .map(|_| {
            let t = Instant::now();
            std::hint::black_box(c.sic(73.0, 0.5, &Cancel::default()));
            t.elapsed()
        })
        .collect();
    repeats.sort();
    println!("SIC m/z 73, median of 15: {:?}", repeats[7]);

    let t = Instant::now();
    let d = decimate(
        &c.scan_times_s,
        &tic,
        f64::NEG_INFINITY,
        f64::INFINITY,
        2000,
    );
    println!("decimate TIC to {} points: {:?}", d.len(), t.elapsed());

    // What the app draws per frame for an overlay of 15 runs showing TIC and SIC.
    let t = Instant::now();
    for _ in 0..15 {
        for ys in [&tic, &sic] {
            std::hint::black_box(decimate(
                &c.scan_times_s,
                ys,
                f64::NEG_INFINITY,
                f64::INFINITY,
                2000,
            ));
        }
    }
    println!("overlay frame, 15 runs x TIC+SIC: {:?}", t.elapsed());

    let t = Instant::now();
    let mut traces = IonTraces::new(&c, 0.5).unwrap();
    println!("ion traces ({}): {:?}", traces.len(), t.elapsed());

    let t = Instant::now();
    let peaks = traces.detect_peaks(&PeakSettings::default());
    println!("peaks in all traces ({}): {:?}", peaks.len(), t.elapsed());

    let t = Instant::now();
    let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
    println!(
        "deconvolution ({} components): {:?}",
        comps.len(),
        t.elapsed()
    );
    assert_eq!(sic.len(), scans);

    // Matched filtering on a copy, so the default results above stay as they are.
    let mut matched = traces.clone();
    let t = Instant::now();
    let m = matched.detect_peaks(&PeakSettings {
        detector: Detector::MATCHED,
        ..PeakSettings::default()
    });
    let mc = deconvolute(&matched, &m, &c.scan_times_s, &DeconvSettings::default());
    println!(
        "matched peaks ({}) + deconvolution ({} components): {:?}",
        m.len(),
        mc.len(),
        t.elapsed()
    );

    // The same run at high resolution: two exact masses 0.04 apart per nominal mass.
    let hr = Chromatogram {
        mz: c
            .mz
            .iter()
            .map(|m| {
                let i = *m - 50.0;
                50.0 + (i / 2.0).floor() + 0.04 * (i % 2.0)
            })
            .collect(),
        ..c.clone()
    };
    let t = Instant::now();
    let mut hr_traces = IonTraces::exact(&hr).unwrap();
    let ht = t.elapsed();
    let mut hp = hr_traces.detect_peaks(&PeakSettings::default());
    reorder_peaks(
        &mut hp,
        &hr_traces,
        10.0,
        0.05,
        PeakSettings::default().number_bins,
    );
    let hc = deconvolute(
        &hr_traces,
        &hp,
        &hr.scan_times_s,
        &DeconvSettings::default(),
    );
    println!(
        "high resolution: traces ({}) {ht:?}, peaks ({}) + reorder + deconvolution ({} components): {:?}",
        hr_traces.len(),
        hp.len(),
        hc.len(),
        t.elapsed()
    );

    // 1 000 of these components against 10 000 library entries: their own spectra plus 9 000 entries of
    // 50 ions each.
    let comps = &comps[..comps.len().min(1000)];
    let msl: String = (0..10_000)
        .map(|k| {
            let peaks: Vec<(usize, usize)> = match comps.get(k) {
                Some(c) => c
                    .spectrum
                    .iter()
                    .map(|p| (p.0 as usize, p.1 as usize))
                    .collect(),
                None => (0..50)
                    .map(|j| (50 + (k * 7 + j * 13) % 500, 1 + (k * j) % 999))
                    .collect(),
            };
            let pairs: String = peaks.iter().map(|(m, v)| format!("({m} {v}) ")).collect();
            format!("NAME: c{k}\nNUM PEAKS: {}\n{pairs}\n\n", peaks.len())
        })
        .collect();
    let lib: Vec<LibraryCompound> = parse_msl(&msl).expect("synthetic MSL");
    for cutoff in [SearchSettings::default().cutoff, 0.0] {
        let s = SearchSettings {
            cutoff,
            ..Default::default()
        };
        let t = Instant::now();
        let prepared = PreparedLibrary::new(&lib, &s);
        let id = map_compounds_to_library(comps, None, &prepared, &s, &Cancel::default()).unwrap();
        println!(
            "library search, {} components x {} entries, cutoff {cutoff} ({} hits, {} assigned): {:?}",
            comps.len(),
            prepared.len(),
            id.hits.iter().map(Vec::len).sum::<usize>(),
            id.assigned.iter().flatten().count(),
            t.elapsed()
        );
    }

    // The batch window's study on files: the run written as .cdf, then three copies, so 1 and 4 runs
    // read the same bytes. Four runs are analysed four at a time.
    let s = ScanSettings::default();
    let prepared = PreparedLibrary::new(&lib, &s.analysis.search);
    let dir = std::env::temp_dir().join(format!("md-perf-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("temp dir");
    let runs: Vec<PathBuf> = (0..4).map(|r| dir.join(format!("run{r}.cdf"))).collect();
    andi::write(&runs[0], &c).expect("write run");
    for r in &runs[1..] {
        std::fs::copy(&runs[0], r).expect("copy run");
    }
    let study_row = |label: &str, runs: &[PathBuf], scan: &ScanSettings| {
        let study = Study {
            runs,
            library: &lib,
            prepared: &prepared,
            scan,
            rescan: true,
            parallel: 4,
        };
        let t = Instant::now();
        let found = run_study(&study, &Progress::default(), &Cancel::default()).unwrap();
        println!(
            "run_study, {label} ({} substances found, {} unknown): {:?}",
            found.runs.iter().map(|r| r.found.len()).sum::<usize>(),
            found.runs.iter().map(|r| r.unknown.len()).sum::<usize>(),
            t.elapsed()
        );
    };
    study_row("1 run targeted", &runs[..1], &s);
    // Untargeted: the same run, also keeping the components without library hit.
    let u = ScanSettings {
        untargeted: true,
        ..ScanSettings::default()
    };
    study_row("1 run untargeted", &runs[..1], &u);
    study_row("4 runs targeted", &runs, &s);
    std::fs::remove_dir_all(&dir).expect("remove temp dir");

    // MID of a labelling experiment: 1 000 substances with two fragments of 8 ions (sum formula
    // C6H12O2 each), found in 200 runs, then the FRAG export of the same table.
    let mut mids_lib: Vec<LibraryCompound> = lib[..1000].to_vec();
    for c in &mut mids_lib {
        for lo in [100, 200] {
            c.formulas.push((lo, lo + 7, "C6H12O2".into()));
            c.mid_sizes.push((lo, lo + 7, 6));
        }
        c.features.insert("FRAG".into(), "(100/8)(200/8)".into());
    }
    let mid_runs: Vec<RunResult> = (0..200)
        .map(|r| RunResult {
            input: r as usize,
            path: format!("run{r}.bin").into(),
            found: (0..mids_lib.len())
                .map(|k| Found {
                    compound: k,
                    score: 1.0,
                    rt_s: k as f64,
                    ions: (100..108)
                        .chain(200..208)
                        .map(|m| IonPeak {
                            mass: m,
                            integral: 1000.0 / (1 + (m + r + k as u32) % 8) as f64,
                            quality: 0.0,
                            signal_to_noise: 0.0,
                        })
                        .collect(),
                    rescanned: false,
                    rescued: Vec::new(),
                    alternatives: Vec::new(),
                })
                .collect(),
            unknown: Vec::new(),
            calibration: None,
        })
        .collect();
    let table = IsotopeTable::standard();
    let t = Instant::now();
    let m = determine(&mid_runs, &mids_lib, &[], &[], Method::Formula, "C", &table);
    println!(
        "MID by sum formula, {} substances x {} runs ({} errors): {:?}",
        m.compounds.len(),
        mid_runs.len(),
        m.errors.len(),
        t.elapsed()
    );
    let refs: Vec<(usize, &LibraryCompound)> = mids_lib.iter().enumerate().collect();
    let t = Instant::now();
    let e = frag_export(&mid_runs, &refs);
    println!(
        "FRAG export ({} KiB): {:?}",
        e.csv.len() / 1024,
        t.elapsed()
    );
    pool_unknowns();
    exact_mass_run(scans, ions);
}

/// Untargeted pooling (`add_unknowns`) with many unknown components: 20 runs share a pool of
/// n / 4 distinct spectra, so the temporary library grows to several thousand entries.
fn pool_unknowns() {
    let mut seed = 12345_u64;
    let mut rnd = move |n: u64| {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (seed >> 33) % n
    };
    for n in [3_000_usize, 6_000, 12_000] {
        let base: Vec<Vec<(u32, u32)>> = (0..n / 4)
            .map(|_| {
                let mut s: Vec<(u32, u32)> = (0..12)
                    .map(|_| (50 + rnd(350) as u32, 100 + rnd(9_900) as u32))
                    .collect();
                s.sort_unstable_by_key(|p| p.0);
                s.dedup_by_key(|p| p.0);
                s
            })
            .collect();
        let mut runs: Vec<RunResult> = (0..20)
            .map(|r| RunResult {
                input: r,
                path: format!("run{r}.cdf").into(),
                unknown: (0..n / 20)
                    .map(|k| Unknown {
                        rt_s: k as f64,
                        ri: None,
                        spectrum: base[rnd(base.len() as u64) as usize].clone(),
                        ions: Vec::new(),
                    })
                    .collect(),
                ..Default::default()
            })
            .collect();
        let s = ScanSettings::default();
        let t = Instant::now();
        let lib = add_unknowns(&mut runs, 0, &s, &Cancel::default()).unwrap();
        println!(
            "pool {n} unknowns in 20 runs ({} groups): {:?}",
            lib.len(),
            t.elapsed()
        );
    }
}

/// Exact-mass run (FastFlight-like): `ions` points per scan spread over 40 exact m/z per nominal
/// mass, each active for 1/40 of the run, so the index has `40 * ions` sparse traces (2.5 % filled).
fn exact_mass_run(scans: usize, ions: usize) {
    let block = scans.div_ceil(40);
    let c = Chromatogram {
        scan_times_s: (0..scans).map(|s| s as f64 / 20.0).collect(),
        offsets: (0..=scans).map(|s| s * ions).collect(),
        mz: (0..scans * ions)
            .map(|k| {
                let (s, j) = (k / ions, k % ions);
                50.0 + j as f32 + 0.001 * ((s / block + j) % 40) as f32
            })
            .collect(),
        intensity: (0..scans * ions)
            .map(|k| {
                let x = ((k / ions + 7 * (k % ions % 20)) % 60) as f32 - 30.0;
                (k % 997) as f32 * 0.05 + 20_000.0 * (-x * x / 8.0).exp()
            })
            .collect(),
        exact_masses: true,
        baseline: Default::default(),
    };
    let t = Instant::now();
    let mut traces = IonTraces::new(&c, 0.5).unwrap();
    let tt = t.elapsed();
    drop(c);
    let t = Instant::now();
    let peaks = traces.detect_all_peaks(&PeakSettings::default());
    let tp = t.elapsed();
    let times: Vec<f64> = (0..scans).map(|s| s as f64 / 20.0).collect();
    let t = Instant::now();
    let comps = deconvolute(&traces, &peaks, &times, &DeconvSettings::default());
    println!(
        "exact-mass run: traces ({} x {scans}) {tt:?}, peaks + reorder ({}) {tp:?}, deconvolution ({} components) {:?}",
        traces.len(),
        peaks.len(),
        comps.len(),
        t.elapsed()
    );
}
