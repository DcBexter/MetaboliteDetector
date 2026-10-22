//! Untargeted groups of synthetic runs: the group numbers and the assignment of `add_unknowns` follow
//! the deconvolution order of the components. The fixture was written at d58c2f3 (before M9), when
//! `analyze_run` did not sort; the apex order is a view rule of the component list now.

use std::fmt::Write as _;
use std::path::PathBuf;

use md_core::Chromatogram;
use md_core::andi;
use md_core::batch::{ScanSettings, Study, run_study};
use md_core::job::{Cancel, Progress};
use md_core::search::PreparedLibrary;

/// Eight masses in every scan. Two substances (apex scans 100 and 104) overlap in one region, and
/// the region lists the lower mass first, so the deconvolution order is not the apex order.
/// Each substance is shifted by `shift` scans.
fn synthetic_run(shift: usize) -> Chromatogram {
    const MASSES: [u32; 8] = [43, 57, 58, 71, 73, 85, 105, 147];
    let substances: [(usize, &[(u32, f64)]); 4] = [
        (100, &[(147, 6000.0), (73, 3000.0)]),
        (104, &[(43, 9000.0), (58, 4000.0)]),
        (200, &[(85, 7000.0), (57, 3000.0)]),
        (60, &[(105, 5000.0), (71, 2500.0)]),
    ];
    let scans = 300;
    let mut seed = 7_u32 + shift as u32;
    let mut intensity = Vec::with_capacity(scans * MASSES.len());
    for s in 0..scans {
        for m in MASSES {
            seed = seed.wrapping_mul(1_103_515_245).wrapping_add(12345);
            let mut v = 500.0 + ((seed >> 16) % 50) as f64;
            for (apex, ions) in &substances {
                if let Some((_, h)) = ions.iter().find(|(mass, _)| *mass == m) {
                    let x = (s as f64 - (apex + shift) as f64) / 2.0;
                    v += h * (-0.5 * x * x).exp();
                }
            }
            intensity.push(v as f32);
        }
    }
    Chromatogram {
        scan_times_s: (0..scans).map(|s| s as f64 * 0.5).collect(),
        offsets: (0..=scans).map(|s| s * MASSES.len()).collect(),
        mz: (0..scans).flat_map(|_| MASSES.map(|m| m as f32)).collect(),
        intensity,
        exact_masses: false,
        baseline: Default::default(),
    }
}

/// The untargeted study of three runs with no library, as text: each pooled group, then per run its
/// components as RT and group index.
fn groups_text() -> String {
    let dir = std::env::temp_dir().join(format!("md-untargeted-order-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let runs: Vec<PathBuf> = (0..3).map(|r| dir.join(format!("run{r}.cdf"))).collect();
    for (r, path) in runs.iter().enumerate() {
        andi::write(path, &synthetic_run(2 * r)).unwrap();
    }
    let mut s = ScanSettings {
        untargeted: true,
        ..Default::default()
    };
    // The synthetic substances have 2 to 3 ions, far below the default of 10 peaks per component.
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
    for run in &result.runs {
        writeln!(text, "run {}", run.input).unwrap();
        for f in &run.found {
            writeln!(text, "{:.1}\t{}", f.rt_s, f.compound).unwrap();
        }
    }
    text
}

#[test]
fn untargeted_groups_match_the_fixture_from_before_the_apex_sort() {
    let text = groups_text();
    assert!(text.contains("run 2"), "{text}");
    assert_eq!(text, include_str!("fixtures/untargeted_groups.tsv"));
}
