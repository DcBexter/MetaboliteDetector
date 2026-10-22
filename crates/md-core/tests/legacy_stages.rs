//! Stage test: each stage of `analyze_run` against the old MetaboliteDetector 3.3 on the 15 example
//! runs. The fixtures `tests/fixtures/legacy/RUN.{peaks,regions,components}.tsv` come from the old
//! binary (docs/reverse-engineering/tools/golden/make-reference.sh --stages, format in stages.py).
//! The regions are recorded for later stage tests; md-core does not expose them yet. Passes without
//! checking anything when example-data/ is missing.

use std::path::Path;

use md_core::analysis::{PeakPass, analyze_run, find_run_peaks};
use md_core::deconv::DeconvSettings;
use md_core::job::Cancel;
use md_core::peaks::PeakSettings;

/// Interpolated apex, in scans: md-core keeps it as f32 (ulp 2.4e-4 at 2048..4096 scans, 4.9e-4 above).
const APEX_TOLERANCE: f64 = 5e-4;
/// S/N, quality and sharpness, relative: md-core computes them in f32 (about 7 digits), the dump prints 9.
const PEAK_VALUE_TOLERANCE: f64 = 1e-6;
/// Model peak values: plain sums of integer intensities, so exact.
const MODEL_TOLERANCE: f64 = 0.0;
/// Ion height: the old code rounds a few modelled heights the other way (as in golden_old_md.rs).
const HEIGHT_TOLERANCE: u32 = 1;
/// Model share and integral, relative: one unit in the 8th significant digit (up to 1e-7) plus the
/// dump rounding. Both codes round the share `c2` to 8 digits, but md-core's normal equations and the
/// old GSL fit land on neighbouring 8th digits for 1938 of the 4366 shares; the integral follows `c2`.
/// Issue #365 (legacy_numeric) is to make these exact.
const FIT_TOLERANCE: f64 = 2e-7;

const RUNS: [&str; 15] = [
    "7235eg04", "7235eg06", "7235eg07", "7235eg08", "7235eg09", "7235eg11", "7235eg12", "7235eg15",
    "7235eg20", "7235eg21", "7235eg22", "7235eg25", "7235eg26", "7235eg30", "7235eg32",
];

fn fixture(run: &str, stage: &str) -> Vec<Vec<String>> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures/legacy")
        .join(format!("{run}.{stage}.tsv"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
    text.lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').map(String::from).collect())
        .collect()
}

fn num(s: &str) -> f64 {
    s.parse().unwrap()
}

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / a.abs().max(b.abs()).max(1e-300)
}

fn run_path(run: &str) -> Option<std::path::PathBuf> {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data")
        .join(format!("{run}.cdf"));
    if !p.exists() {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            p.display()
        );
    }
    p.exists().then_some(p)
}

#[test]
fn fixtures_are_complete() {
    for run in RUNS {
        let peaks = fixture(run, "peaks").len();
        let regions = fixture(run, "regions");
        assert!(peaks > 0 && !regions.is_empty(), "{run}");
        // Every region lists valid peak lines.
        for r in regions {
            assert!(
                r[2].split(',').all(|k| k.parse::<usize>().unwrap() < peaks),
                "{run}"
            );
        }
        assert!(
            fixture(run, "components").iter().any(|l| l[0] == "model"),
            "{run}"
        );
    }
}

/// `find_run_peaks` (main window pass) against the old `detectAllPeaks` list, peak by peak and in order.
#[test]
fn peaks_match_old_version() {
    for run in RUNS {
        let Some(path) = run_path(run) else { return };
        let c = md_core::andi::read(path).unwrap();
        let p = find_run_peaks(
            c,
            &PeakSettings::default(),
            PeakPass::Reordered,
            &Cancel::default(),
        )
        .unwrap();
        let old = fixture(run, "peaks");
        assert_eq!(old.len(), p.peaks.len(), "{run}: peak count");
        for (k, (o, n)) in old.iter().zip(&p.peaks).enumerate() {
            let at = format!("{run}: peak {k} ({o:?} vs {n:?})");
            assert_eq!(
                num(&o[0]),
                f64::from(p.traces.masses[n.row as usize]),
                "{at}"
            );
            let scans = (n.first_scan, n.apex_scan, n.last_scan);
            assert_eq!(
                (num(&o[1]), num(&o[2]), num(&o[3])),
                (scans.0.into(), scans.1.into(), scans.2.into()),
                "{at}"
            );
            assert!(
                (num(&o[4]) - n.interpolated_max_scan()).abs() <= APEX_TOLERANCE,
                "{at}: apex"
            );
            for (old, new) in [
                (&o[5], n.signal_to_noise),
                (&o[6], n.quality),
                (&o[7], n.sharpness),
            ] {
                assert!(
                    rel(num(old), f64::from(new)) <= PEAK_VALUE_TOLERANCE,
                    "{at}: {old} vs {new}"
                );
            }
        }
    }
}

/// `analyze_run` components against the old compounds, in order: model peak, spectrum and fit.
#[test]
#[allow(
    clippy::cognitive_complexity,
    reason = "a test of many cases: each assert counts, and Sonar skips tests"
)]
fn components_match_old_version() {
    let (mut fitted, mut eighth_digit) = (0, 0);
    for run in RUNS {
        let Some(path) = run_path(run) else { return };
        let c = md_core::andi::read(path).unwrap();
        let a = analyze_run(
            c,
            &PeakSettings::default(),
            &DeconvSettings::default(),
            PeakPass::Reordered,
            &Cancel::default(),
        )
        .unwrap();
        let mut old = fixture(run, "components").into_iter().peekable();
        let mut comps = a.components.iter();
        while let Some(head) = old.next() {
            assert_eq!(head[0], "component");
            let n = comps
                .next()
                .unwrap_or_else(|| panic!("{run}: missing component at {}", head[1]));
            let at = format!("{run}: component {}", head[1]);
            assert_eq!(num(&head[1]), f64::from(n.apex_scan), "{at}");
            let m = old.next().unwrap();
            assert_eq!(m[0], "model");
            assert_eq!(
                (num(&m[1]), num(&m[2]), num(&m[3])),
                (n.first_scan.into(), n.apex_scan.into(), n.last_scan.into()),
                "{at}: model"
            );
            let values: Vec<f64> = m[4].split(',').map(num).collect();
            assert_eq!(values.len(), n.model.len(), "{at}: model length");
            assert!(
                values
                    .iter()
                    .zip(&n.model)
                    .all(|(o, n)| (o - n).abs() <= MODEL_TOLERANCE),
                "{at}: model values"
            );
            let mut ions = Vec::new();
            while let Some(l) = old.next_if(|l| l[0] == "ion") {
                ions.push(l);
            }
            assert_eq!(ions.len(), n.exact_spectrum.len(), "{at}: ion count");
            for (o, &(mass, height)) in ions.iter().zip(&n.exact_spectrum) {
                assert_eq!(num(&o[1]), f64::from(mass), "{at}");
                assert!(
                    (num(&o[2]) as u32).abs_diff(height) <= HEIGHT_TOLERANCE,
                    "{at}: ion {mass} height {} vs {height}",
                    o[2]
                );
                // The fitted peak of this ion: same mass and bounds (the old peak copy).
                let f = n
                    .peaks
                    .iter()
                    .find(|f| {
                        let p = &a.peaks[f.peak];
                        f64::from(a.masses[f.peak]) == num(&o[1])
                            && (
                                f64::from(p.first_scan),
                                f64::from(p.apex_scan),
                                f64::from(p.last_scan),
                            ) == (num(&o[3]), num(&o[4]), num(&o[5]))
                    })
                    .unwrap_or_else(|| panic!("{at}: no fitted peak for {o:?}"));
                let (di, ds) = (rel(num(&o[6]), f.integral), rel(num(&o[7]), f.model_share));
                assert!(
                    di <= FIT_TOLERANCE && ds <= FIT_TOLERANCE,
                    "{at}: ion {mass} integral {} vs {}, share {} vs {}",
                    o[6],
                    f.integral,
                    o[7],
                    f.model_share
                );
                // The dump's 9 digits resolve the 8-digit share exactly: count the neighbours.
                if ds > 1e-9 {
                    eighth_digit += 1;
                }
                fitted += 1;
            }
        }
        assert!(comps.next().is_none(), "{run}: extra components");
    }
    eprintln!("{eighth_digit} of {fitted} model shares differ in the 8th digit");
    // Known count (architecture.md, "Stage fixtures"): fails when it grows; lower it when #365 fixes some.
    assert_eq!(fitted, 4366);
    assert!(
        eighth_digit <= 1938,
        "{eighth_digit} shares off, 1938 known"
    );
}
