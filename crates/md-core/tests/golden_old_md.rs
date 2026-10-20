//! Golden test: peaks and deconvolution against the old MetaboliteDetector 3.3 on the example runs,
//! with the default and the matched peak detector. The references `tests/data/old-md-3.3*.tsv` come
//! from the old binary itself
//! (docs/reverse-engineering/tools/golden/make-reference.sh); tolerances are explained in
//! docs/architecture.md ("Golden test"). Passes without checking anything when example-data/ is
//! missing. `--nocapture` prints every difference.

use std::collections::BTreeSet;
use std::path::Path;

use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::peaks::{Detector, IonTraces, PeakSettings};

/// The reference prints S/N with 6 significant digits.
const SN_TOLERANCE: f64 = 1e-5;
/// Heights: the old code rounds some modelled apex heights the other way.
const HEIGHT_TOLERANCE: u32 = 1;
/// Integrals: the reference prints 7 significant digits.
const INTEGRAL_TOLERANCE: f64 = 1e-5;

/// Components the port may miss, pair with other ions or apex, or give other values, and why.
struct Known {
    missing: usize,
    other_ions: usize,
    other_values: usize,
}

#[derive(Default)]
struct RefRun {
    name: String,
    /// (ion, first, apex, last, S/N)
    peaks: Vec<(u32, u32, u32, u32, f64)>,
    comps: Vec<RefComp>,
}

struct RefComp {
    scan: u32,
    /// (mass, height, integral)
    ions: Vec<(u32, u32, f64)>,
}

fn reference(text: &str) -> Vec<RefRun> {
    let mut runs: Vec<RefRun> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let u = |i: usize| f[i].parse::<u32>().unwrap();
        let x = |i: usize| f[i].parse::<f64>().unwrap();
        match f[0] {
            "run" => runs.push(RefRun {
                name: f[1].into(),
                ..Default::default()
            }),
            "peak" => runs
                .last_mut()
                .unwrap()
                .peaks
                .push((u(1), u(2), u(3), u(4), x(5))),
            "component" => runs.last_mut().unwrap().comps.push(RefComp {
                scan: u(1),
                ions: Vec::new(),
            }),
            "ion" => {
                let c = runs.last_mut().unwrap().comps.last_mut().unwrap();
                c.ions.push((u(1), u(2), x(3)));
            }
            other => panic!("unknown record {other}"),
        }
    }
    runs
}

#[derive(Default, Debug)]
struct Count {
    same: usize,
    missing: usize,
    extra: usize,
    other_ions: usize,
    other_values: usize,
}

#[test]
fn matches_old_version() {
    // Precise-apex ties the old GSL fit rounds up (docs/architecture.md).
    let known = Known {
        missing: 1,
        other_ions: 3,
        other_values: 1,
    };
    compare(
        include_str!("data/old-md-3.3.tsv"),
        &PeakSettings::default(),
        known,
    );
}

#[test]
fn matches_old_version_with_matched_filtering() {
    let known = Known {
        missing: 0,
        other_ions: 0,
        other_values: 0,
    };
    let peaks = PeakSettings {
        detector: Detector::MATCHED,
        ..PeakSettings::default()
    };
    compare(include_str!("data/old-md-3.3-matched.tsv"), &peaks, known);
}

fn compare(text: &str, settings: &PeakSettings, known: Known) {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    if !dir.is_dir() {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            dir.display()
        );
        return;
    }
    let runs = reference(text);
    assert_eq!(runs.len(), 15);
    let (mut peak_runs, mut total) = (0, Count::default());
    for r in &runs {
        let c = md_core::andi::read(dir.join(&r.name)).unwrap();
        let mut traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(settings);
        let comps = deconvolute(&traces, &peaks, &DeconvSettings::default());
        let mass = |row: u32| traces.masses[row as usize];

        // Peaks: the same peaks with the same bounds and S/N.
        if !r.peaks.is_empty() {
            peak_runs += 1;
            let mut new: Vec<_> = peaks
                .iter()
                .map(|p| (mass(p.row), p.first_scan, p.apex_scan, p.last_scan))
                .collect();
            let mut old: Vec<_> = r.peaks.iter().map(|p| (p.0, p.1, p.2, p.3)).collect();
            new.sort();
            old.sort();
            assert_eq!(old, new, "{}: peaks (ion, first, apex, last)", r.name);
            for &(ion, _, apex, _, sn) in &r.peaks {
                let p = peaks
                    .iter()
                    .find(|p| mass(p.row) == ion && p.apex_scan == apex)
                    .unwrap();
                // The matched detector leaves S/N at 0 in both versions.
                assert!(
                    (p.signal_to_noise as f64 - sn).abs() <= SN_TOLERANCE * sn,
                    "{}: ion {ion} scan {apex} S/N {sn} vs {}",
                    r.name,
                    p.signal_to_noise
                );
            }
        }

        // Components: pair each old one with the unused new one at (almost) the same scan that
        // shares the most ions (more than half); several components can sit on the same scan.
        let mut used = vec![false; comps.len()];
        let mut n = Count::default();
        for o in &r.comps {
            let old_ions: BTreeSet<u32> = o.ions.iter().map(|i| i.0).collect();
            let best = (0..comps.len())
                .filter(|j| !used[*j] && comps[*j].apex_scan.abs_diff(o.scan) <= 2)
                .map(|j| {
                    let ions: BTreeSet<u32> = comps[j].spectrum.iter().map(|s| s.0).collect();
                    (old_ions.intersection(&ions).count(), j)
                })
                .max_by_key(|(c, j)| (*c, std::cmp::Reverse(comps[*j].apex_scan.abs_diff(o.scan))));
            let Some((_, j)) = best.filter(|(c, _)| 2 * c > old_ions.len()) else {
                n.missing += 1;
                eprintln!(
                    "{}: old component at scan {} {old_ions:?} missing",
                    r.name, o.scan
                );
                continue;
            };
            used[j] = true;
            let k = &comps[j];
            let new_ions: BTreeSet<u32> = k.spectrum.iter().map(|s| s.0).collect();
            if new_ions != old_ions || k.apex_scan != o.scan {
                n.other_ions += 1;
                eprintln!(
                    "{}: component scan {} -> {}: only old {:?}, only new {:?}",
                    r.name,
                    o.scan,
                    k.apex_scan,
                    old_ions.difference(&new_ions).collect::<Vec<_>>(),
                    new_ions.difference(&old_ions).collect::<Vec<_>>()
                );
                continue;
            }
            // Both spectra are sorted by mass, so they pair up ion by ion.
            let off: Vec<String> = o
                .ions
                .iter()
                .zip(&k.spectrum)
                .filter_map(|(&(ion, h, integral), &(_, nh))| {
                    let f = k.peaks.iter().find(|f| mass(peaks[f.peak].row) == ion)?;
                    let di = (f.integral - integral).abs() / integral.abs().max(1.0);
                    (h.abs_diff(nh) > HEIGHT_TOLERANCE || di > INTEGRAL_TOLERANCE).then(|| {
                        format!(
                            "{ion}: height {h} vs {nh}, integral {integral} vs {:.1}",
                            f.integral
                        )
                    })
                })
                .collect();
            if off.is_empty() {
                n.same += 1;
            } else {
                n.other_values += 1;
                eprintln!("{}: component scan {}: {}", r.name, o.scan, off.join("; "));
            }
        }
        n.extra = used.iter().filter(|u| !**u).count();
        eprintln!("{}: {} old components: {n:?}", r.name, r.comps.len());
        total.same += n.same;
        total.missing += n.missing;
        total.extra += n.extra;
        total.other_ions += n.other_ions;
        total.other_values += n.other_values;
    }
    eprintln!("total: {total:?}");
    assert_eq!(peak_runs, 3);
    assert_eq!(total.extra, 0, "new components without an old partner");
    assert!(total.missing <= known.missing, "{total:?}");
    assert!(total.other_ions <= known.other_ions, "{total:?}");
    assert!(total.other_values <= known.other_values, "{total:?}");
}
