//! Performance smoke test on a synthetic large run.
//! `cargo run --release -p md-core --example perf -- [scans] [ions_per_scan]`
//! Default: 60 min at 20 scans/s with 500 ions per scan (36 000 scans, 18 M points, 300 000 peaks,
//! 12 000 substances).

use std::time::Instant;

use md_core::batch::{Found, IonPeak, RunResult, ScanSettings, scan_chromatogram};
use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::isotope::IsotopeTable;
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
    let sic = c.sic(73.0, 0.5);
    println!("SIC m/z 73: {:?}", t.elapsed());

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
    let mut traces = IonTraces::new(&c).unwrap();
    println!("ion traces ({}): {:?}", traces.len(), t.elapsed());

    let t = Instant::now();
    let peaks = traces.detect_peaks(&PeakSettings::default());
    println!("peaks in all traces ({}): {:?}", peaks.len(), t.elapsed());

    let t = Instant::now();
    let comps = deconvolute(&traces, &peaks, &DeconvSettings::default());
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
    let mc = deconvolute(&matched, &m, &DeconvSettings::default());
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
    let hc = deconvolute(&hr_traces, &hp, &DeconvSettings::default());
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
        let id = map_compounds_to_library(comps, None, &prepared, &s);
        println!(
            "library search, {} components x {} entries, cutoff {cutoff} ({} hits, {} assigned): {:?}",
            comps.len(),
            prepared.len(),
            id.hits.iter().map(Vec::len).sum::<usize>(),
            id.assigned.iter().flatten().count(),
            t.elapsed()
        );
    }

    // Batch: one run end to end (ion traces, peaks, deconvolution, search of all components).
    let s = ScanSettings::default();
    let prepared = PreparedLibrary::new(&lib, &s.search);
    let t = Instant::now();
    let found = scan_chromatogram(c.clone(), &prepared, &s).unwrap();
    println!(
        "batch, one run end to end ({} substances found): {:?}",
        found.found.len(),
        t.elapsed()
    );
    // Untargeted: the same run, also keeping the components without library hit.
    let u = ScanSettings {
        untargeted: true,
        ..s
    };
    let t = Instant::now();
    let both = scan_chromatogram(c.clone(), &prepared, &u).unwrap();
    println!(
        "batch untargeted, one run ({} substances, {} unknown): {:?}",
        both.found.len(),
        both.unknown.len(),
        t.elapsed()
    );

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
}
