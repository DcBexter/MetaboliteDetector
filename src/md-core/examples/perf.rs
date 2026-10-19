//! Performance smoke test on a synthetic large run.
//! `cargo run --release -p md-core --example perf -- [scans] [ions_per_scan]`
//! Default: 60 min at 20 scans/s with 500 ions per scan (36 000 scans, 18 M points, 300 000 peaks,
//! 12 000 substances).

use std::time::Instant;

use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::library::{LibraryCompound, parse_msl};
use md_core::peaks::{IonTraces, PeakSettings};
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
    let traces = IonTraces::new(&c).unwrap();
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
}
