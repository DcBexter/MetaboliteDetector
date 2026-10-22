//! Untargeted groups of synthetic exact-mass runs: the group numbers and the assignment of
//! `add_unknowns` follow the deconvolution order of the components, not their apex order. The fixture
//! was written at d58c2f3 (before M9), when `analyze_run` did not sort. The apex order is a view rule
//! of the component list now. On the 15 example runs and on 400 random nominal-mass runs the two orders
//! agree (0 inversions), so #422 changed nothing on real data; these runs are the ones where they differ
//! (about 2 % of random exact-mass runs have an inversion).

use std::fmt::Write as _;
use std::path::PathBuf;

use md_core::Chromatogram;
use md_core::andi;
use md_core::batch::{ScanSettings, Study, run_study};
use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::job::{Cancel, Progress};
use md_core::peaks::{IonTraces, PeakSettings};
use md_core::search::PreparedLibrary;

/// Seeds whose runs have an apex inversion in the deconvolution order (see the module doc).
const SEEDS: [u32; 3] = [310, 323, 434];

/// Six exact masses, a few substances with apexes 0 to 12 scans apart, each drawn from a seed.
fn synthetic_run(seed0: u32) -> Chromatogram {
    const MASSES: [f32; 6] = [43.02, 57.03, 71.047, 73.047, 73.089, 147.05];
    let scans = 240usize;
    let mut seed = seed0;
    let mut rnd = move |n: u32| {
        seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
        (seed >> 8) % n
    };
    let nsub = 3 + rnd(3) as usize;
    let base = 90.0 + rnd(60) as f64;
    // Substance: apex scan, peak width in scans, (mass index, height) per ion.
    type Substance = (f64, f64, Vec<(usize, f64)>);
    let mut subs: Vec<Substance> = Vec::new();
    for _ in 0..nsub {
        let apex = base + rnd(12) as f64;
        let sigma = [1.0, 1.5, 2.0, 4.0][rnd(4) as usize];
        let nion = 2 + rnd(2) as usize;
        let ions = (0..nion)
            .map(|_| {
                (
                    rnd(MASSES.len() as u32) as usize,
                    2000.0 + rnd(20000) as f64,
                )
            })
            .collect();
        subs.push((apex, sigma, ions));
    }
    let mut intensity = Vec::with_capacity(scans * MASSES.len());
    for s in 0..scans {
        for mi in 0..MASSES.len() {
            seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
            let mut v = 300.0 + ((seed >> 16) % 60) as f64;
            for (apex, sigma, ions) in &subs {
                for (m, h) in ions {
                    if *m == mi {
                        let x = (s as f64 - apex) / sigma;
                        v += h * (-0.5 * x * x).exp();
                    }
                }
            }
            intensity.push(v as f32);
        }
    }
    Chromatogram {
        scan_times_s: (0..scans).map(|s| s as f64 * 0.5).collect(),
        offsets: (0..=scans).map(|s| s * MASSES.len()).collect(),
        mz: (0..scans).flat_map(|_| MASSES).collect(),
        intensity,
        exact_masses: true,
        baseline: Default::default(),
    }
}

/// Whether the deconvolution order of `c` (default settings, one peak per component allowed) is not
/// the apex order.
fn has_inversion(c: &Chromatogram) -> bool {
    let d = DeconvSettings {
        min_peak_number: 1,
        ..Default::default()
    };
    let mut traces = IonTraces::new(c, 0.5).unwrap();
    let peaks = traces.detect_all_peaks(&PeakSettings::default());
    let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &d);
    let apex: Vec<u32> = comps.iter().map(|k| k.apex_scan).collect();
    apex.windows(2).any(|w| w[0] > w[1])
}

/// The untargeted study of the runs with no library, as text: each pooled group, then per run its
/// components as RT and group index.
fn groups_text() -> String {
    let dir = std::env::temp_dir().join(format!("md-untargeted-order-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let runs: Vec<PathBuf> = (0..SEEDS.len())
        .map(|r| dir.join(format!("run{r}.cdf")))
        .collect();
    for (path, seed) in runs.iter().zip(SEEDS) {
        andi::write(path, &synthetic_run(seed)).unwrap();
    }
    let mut s = ScanSettings {
        untargeted: true,
        ..Default::default()
    };
    s.analysis.deconv.min_peak_number = 1;
    let prepared = PreparedLibrary::new(&[], &s.analysis.search);
    let study = Study {
        runs: &runs,
        library: &[],
        prepared: &prepared,
        scan: &s,
        rescan: true,
        parallel: 1,
    };
    let result = run_study(&study, &Progress::default(), &Cancel::default()).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    assert!(result.errors.is_empty(), "{:?}", result.errors);
    let mut text = String::new();
    for (k, group) in result.extra.iter().enumerate() {
        writeln!(text, "group {k}\t{}", group.name).unwrap();
    }
    for (input, run) in result.runs.iter().enumerate() {
        writeln!(text, "run {input}").unwrap();
        for f in &run.found {
            writeln!(text, "{:.1}\t{}", f.rt_s, f.compound).unwrap();
        }
    }
    text
}

#[test]
fn untargeted_groups_match_the_fixture_from_before_the_apex_sort() {
    // The runs must reach the study in deconvolution order that differs from the apex order.
    assert!(SEEDS.iter().all(|&s| has_inversion(&synthetic_run(s))));
    let text = groups_text();
    // The study keeps that order: some run lists a later RT before an earlier one.
    let descends = text.split("run ").skip(1).any(|run| {
        let rts: Vec<f64> = run
            .lines()
            .skip(1)
            .filter_map(|l| l.split('\t').next()?.parse().ok())
            .collect();
        rts.windows(2).any(|w| w[0] > w[1])
    });
    assert!(descends, "{text}");
    assert_eq!(text, include_str!("fixtures/untargeted_groups.tsv"));
}
