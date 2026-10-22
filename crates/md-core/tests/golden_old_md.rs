//! Golden test: peaks and deconvolution against the old MetaboliteDetector 3.3 on the example runs,
//! with the default and the matched peak detector, and on a high-resolution run made from 7235eg04
//! (as CDF, and as a cache with exact masses like FastFlight writes). The references
//! `tests/data/old-md-3.3*.tsv` come from the old binary itself
//! (docs/reverse-engineering/tools/golden/make-reference.sh); tolerances are explained in
//! docs/architecture.md ("Golden test"). Passes without checking anything when example-data/ is
//! missing. `--nocapture` prints every difference.

use std::collections::BTreeSet;
use std::path::Path;

use md_core::Chromatogram;
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
    extra: usize,
    other_ions: usize,
    other_values: usize,
}

#[derive(Default)]
struct RefRun {
    name: String,
    /// (ion, first, apex, last, S/N)
    peaks: Vec<(f32, u32, u32, u32, f64)>,
    comps: Vec<RefComp>,
}

struct RefComp {
    scan: u32,
    /// (exact mass, height, integral)
    ions: Vec<(f32, u32, f64)>,
}

fn reference(text: &str) -> Vec<RefRun> {
    let mut runs: Vec<RefRun> = Vec::new();
    for line in text.lines().filter(|l| !l.starts_with('#')) {
        let f: Vec<&str> = line.split('\t').collect();
        let u = |i: usize| f[i].parse::<u32>().unwrap();
        let x = |i: usize| f[i].parse::<f64>().unwrap();
        let m = |i: usize| f[i].parse::<f32>().unwrap();
        match f[0] {
            "run" => runs.push(RefRun {
                name: f[1].into(),
                ..Default::default()
            }),
            "peak" => runs
                .last_mut()
                .unwrap()
                .peaks
                .push((m(1), u(2), u(3), u(4), x(5))),
            "component" => runs.last_mut().unwrap().comps.push(RefComp {
                scan: u(1),
                ions: Vec::new(),
            }),
            "ion" => {
                let c = runs.last_mut().unwrap().comps.last_mut().unwrap();
                c.ions.push((m(1), u(2), x(3)));
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
    let known = Known {
        missing: 0,
        extra: 0,
        other_ions: 0,
        other_values: 0,
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
        extra: 0,
        other_ions: 0,
        other_values: 0,
    };
    let peaks = PeakSettings {
        detector: Detector::MATCHED,
        ..PeakSettings::default()
    };
    compare(include_str!("data/old-md-3.3-matched.tsv"), &peaks, known);
}

/// The high-resolution run of `matches_old_version_on_a_high_resolution_run`, made from 7235eg04:
/// every m/z moves up by 0.047, and a second fragment at 147.089 repeats the 147 trace 5 scans later.
/// The old CDF import centroids such scans to nominal masses (algorithms.md §4.1).
const HIGH_RESOLUTION_RUN: &str = "7235eg04-hr.cdf";

fn high_resolution_run(c: &Chromatogram) -> Chromatogram {
    let mut out = Chromatogram {
        scan_times_s: c.scan_times_s.clone(),
        offsets: vec![0],
        mz: Vec::new(),
        intensity: Vec::new(),
        exact_masses: false,
        baseline: Default::default(),
    };
    let at_147 = |s: usize| {
        let (m, i) = c.spectrum(s);
        m.iter().position(|x| *x == 147.0).map_or(0.0, |k| i[k])
    };
    for s in 0..c.scan_count() {
        let (m, i) = c.spectrum(s);
        let mut points: Vec<(f32, f32)> =
            m.iter().map(|x| x + 0.047).zip(i.iter().copied()).collect();
        points.push((147.089, s.checked_sub(5).map_or(0.0, at_147)));
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        out.mz.extend(points.iter().map(|p| p.0));
        out.intensity.extend(points.iter().map(|p| p.1));
        out.offsets.push(out.mz.len());
    }
    out
}

/// The same run as a disk-scan cache (`.bin` + `.idx`) with DataType 1, as the old FastFlight import
/// writes it: the old program keeps its exact masses and runs reorderPeaks (algorithms.md §4.1).
const EXACT_MASS_RUN: &str = "7235eg04-exact.bin";

/// Writes `c` as the old `GCMSDiskScanWriter` lays out a cache (file-formats.md §2), every scan with
/// DataType 1, no features, no TOF polynomial and no RI calibration.
fn write_cache(bin_path: &Path, c: &Chromatogram) {
    let (mut idx, mut bin) = (Vec::new(), Vec::new());
    let u32s = |v: &mut Vec<u8>, xs: &[u32]| xs.iter().for_each(|x| v.extend(x.to_le_bytes()));
    idx.extend(0x4b49_4104_i64.to_le_bytes());
    u32s(&mut idx, &[4, 8, 8, c.scan_count() as u32, 4, 4, 0, 0]);
    idx.extend([0; 40]); // global baseline statistics
    let none = md_core::ri::StoredCalibration {
        points_ms: Vec::new(),
        drt_ms: 0.0,
        source: b"not calibrated".to_vec(),
    };
    idx.extend(none.to_bytes());
    let mut ions = c.mz.clone();
    ions.sort_by(f32::total_cmp);
    ions.dedup();
    u32s(&mut idx, &[ions.len() as u32]);
    ions.iter().for_each(|m| idx.extend(m.to_le_bytes()));
    for s in 0..c.scan_count() {
        let (mz, v) = c.spectrum(s);
        let (pos, n) = (bin.len() as u64, mz.len());
        let masses_ptr = pos + 36 + 4 * n as u64;
        let tic: f64 = v.iter().map(|x| f64::from(*x)).sum();
        let rt_ms = c.scan_times_s[s] * 1000.0;
        bin.extend(((36 + 8 * n) as i32).to_le_bytes());
        bin.extend(masses_ptr.to_le_bytes());
        bin.extend(tic.to_le_bytes());
        u32s(&mut bin, &[s as u32]);
        bin.extend(rt_ms.to_le_bytes());
        u32s(&mut bin, &[n as u32]);
        v.iter().for_each(|x| bin.extend((*x as i32).to_le_bytes()));
        mz.iter().for_each(|m| bin.extend(m.to_le_bytes()));
        u32s(&mut idx, &[s as u32]);
        idx.extend(pos.to_le_bytes());
        idx.extend(masses_ptr.to_le_bytes());
        idx.extend(tic.to_le_bytes());
        idx.extend(rt_ms.to_le_bytes());
        u32s(&mut idx, &[n as u32, 1]);
        // BaselineStatistics (mean, sd, msd, nf, n): the old app takes the median per-scan nf as its
        // noise factor (algorithms.md §2.3); mean = sd = 1 also moved peak edges in the old app. nf = 1
        // with mean and sd 0 is what the port assumes (architecture.md, "Noise factor 1.0"). Real FastFlight caches carry
        // `analyzeBaseline(values, 100)` here instead.
        [0.0f64, 0.0, 0.0, 1.0]
            .iter()
            .for_each(|x| idx.extend(x.to_le_bytes()));
        idx.extend(100i64.to_le_bytes());
    }
    std::fs::write(bin_path, bin).unwrap();
    std::fs::write(bin_path.with_extension("idx"), idx).unwrap();
}

/// The exact-mass run as md-core reads it: written as a cache and read back.
fn exact_mass_run(c: &Chromatogram) -> Chromatogram {
    let dir = std::env::temp_dir().join(format!("md-core-golden-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let bin = dir.join(EXACT_MASS_RUN);
    write_cache(&bin, &high_resolution_run(c));
    let run = md_core::diskscan::read(&bin).unwrap();
    std::fs::remove_dir_all(dir).ok();
    assert!(run.exact_masses, "DataType 1 means exact masses");
    run
}

/// Writes the exact-mass cache for `make-reference.sh --exact-masses` to `$MD_HIGH_RESOLUTION_OUT`.
#[test]
#[ignore]
fn write_exact_mass_run() {
    let out = std::env::var("MD_HIGH_RESOLUTION_OUT").expect("MD_HIGH_RESOLUTION_OUT");
    let c = md_core::andi::read(example_dir().join("7235eg04.cdf")).unwrap();
    write_cache(
        &Path::new(&out).join(EXACT_MASS_RUN),
        &high_resolution_run(&c),
    );
}

#[test]
fn matches_old_version_on_an_exact_mass_run() {
    let known = Known {
        missing: 0,
        extra: 0,
        other_ions: 0,
        other_values: 0,
    };
    compare(
        include_str!("data/old-md-3.3-exact-masses.tsv"),
        &PeakSettings::default(),
        known,
    );
}

/// Writes the high-resolution run for make-reference.sh to `$MD_HIGH_RESOLUTION_OUT`.
#[test]
#[ignore]
fn write_high_resolution_run() {
    let out = std::env::var("MD_HIGH_RESOLUTION_OUT").expect("MD_HIGH_RESOLUTION_OUT");
    let c = md_core::andi::read(example_dir().join("7235eg04.cdf")).unwrap();
    md_core::andi::write(
        Path::new(&out).join(HIGH_RESOLUTION_RUN),
        &high_resolution_run(&c),
    )
    .unwrap();
}

#[test]
fn matches_old_version_on_a_high_resolution_run() {
    let known = Known {
        missing: 0,
        extra: 0,
        other_ions: 0,
        other_values: 0,
    };
    compare(
        include_str!("data/old-md-3.3-high-resolution.tsv"),
        &PeakSettings::default(),
        known,
    );
}

fn example_dir() -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data")
}

#[allow(
    clippy::cognitive_complexity,
    reason = "a test of many cases: each assert counts, and Sonar skips tests"
)]
fn compare(text: &str, settings: &PeakSettings, known: Known) {
    let dir = example_dir();
    if !dir.is_dir() {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            dir.display()
        );
        return;
    }
    let runs = reference(text);
    assert!(runs.len() == 15 || [HIGH_RESOLUTION_RUN, EXACT_MASS_RUN].contains(&&*runs[0].name));
    let (mut peak_runs, mut total) = (0, Count::default());
    for r in &runs {
        let c = if r.name == HIGH_RESOLUTION_RUN {
            high_resolution_run(&md_core::andi::read(dir.join("7235eg04.cdf")).unwrap())
        } else if r.name == EXACT_MASS_RUN {
            exact_mass_run(&md_core::andi::read(dir.join("7235eg04.cdf")).unwrap())
        } else {
            md_core::andi::read(dir.join(&r.name)).unwrap()
        };
        let mut traces =
            IonTraces::new(&c, md_core::settings::AnalysisSettings::default().fwhm).unwrap();
        let peaks = traces.detect_all_peaks(settings);
        let comps = deconvolute(&traces, &peaks, &c.scan_times_s, &DeconvSettings::default());
        let mass = |row: u32| traces.mz[row as usize];

        // Peaks: the same peaks with the same bounds and S/N.
        if !r.peaks.is_empty() {
            peak_runs += 1;
            let mut new: Vec<_> = peaks
                .iter()
                .map(|p| {
                    (
                        mass(p.row).to_bits(),
                        p.first_scan,
                        p.apex_scan,
                        p.last_scan,
                    )
                })
                .collect();
            let mut old: Vec<_> = r
                .peaks
                .iter()
                .map(|p| (p.0.to_bits(), p.1, p.2, p.3))
                .collect();
            new.sort();
            old.sort();
            let only = |a: &[(u32, u32, u32, u32)], b: &[(u32, u32, u32, u32)]| -> Vec<_> {
                a.iter()
                    .filter(|p| b.binary_search(p).is_err())
                    .map(|p| (f32::from_bits(p.0), p.1, p.2, p.3))
                    .collect()
            };
            assert!(
                old == new,
                "{}: peaks (ion, first, apex, last) only old {:?}, only new {:?}",
                r.name,
                only(&old, &new),
                only(&new, &old)
            );
            for &(ion, _, apex, _, sn) in &r.peaks {
                let p = peaks
                    .iter()
                    .find(|p| mass(p.row) == ion && p.apex_scan == apex)
                    .unwrap();
                // The matched detector leaves S/N at the old default -1.
                assert!(
                    (p.signal_to_noise as f64 - sn).abs() <= SN_TOLERANCE * sn.abs(),
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
            // Masses as f32 bits: exact on the exact-mass run, nominal on the others.
            let old_ions: BTreeSet<u32> = o.ions.iter().map(|i| i.0.to_bits()).collect();
            let best = (0..comps.len())
                .filter(|j| !used[*j] && comps[*j].apex_scan.abs_diff(o.scan) <= 2)
                .map(|j| {
                    let ions: BTreeSet<u32> = (comps[j].exact_spectrum.iter())
                        .map(|s| s.0.to_bits())
                        .collect();
                    (old_ions.intersection(&ions).count(), j)
                })
                .max_by_key(|(c, j)| (*c, std::cmp::Reverse(comps[*j].apex_scan.abs_diff(o.scan))));
            let Some((_, j)) = best.filter(|(c, _)| 2 * c > old_ions.len()) else {
                n.missing += 1;
                eprintln!(
                    "{}: old component at scan {} {:?} missing",
                    r.name,
                    o.scan,
                    o.ions.iter().map(|i| i.0).collect::<Vec<_>>()
                );
                continue;
            };
            used[j] = true;
            let k = &comps[j];
            let new_ions: BTreeSet<u32> = k.exact_spectrum.iter().map(|s| s.0.to_bits()).collect();
            let floats = |d: Vec<&u32>| {
                d.into_iter()
                    .map(|b| f32::from_bits(*b))
                    .collect::<Vec<_>>()
            };
            if new_ions != old_ions || k.apex_scan != o.scan {
                n.other_ions += 1;
                eprintln!(
                    "{}: component scan {} -> {}: only old {:?}, only new {:?}",
                    r.name,
                    o.scan,
                    k.apex_scan,
                    floats(old_ions.difference(&new_ions).collect()),
                    floats(new_ions.difference(&old_ions).collect())
                );
                continue;
            }
            // Each mass against the old entry, the integral of its sharpest peak.
            let off: Vec<String> = k
                .exact_spectrum
                .iter()
                .filter_map(|&(ion, nh)| {
                    let f = k.peaks.iter().find(|f| mass(peaks[f.peak].row) == ion)?;
                    let Some(&(_, h, integral)) = o.ions.iter().find(|i| i.0 == ion) else {
                        return Some(format!("{ion}: no old ion of this exact mass"));
                    };
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
    assert_eq!(peak_runs, runs.len().min(3));
    assert!(
        total.extra <= known.extra,
        "new components without an old partner: {total:?}"
    );
    assert!(total.missing <= known.missing, "{total:?}");
    assert!(total.other_ions <= known.other_ions, "{total:?}");
    assert!(total.other_values <= known.other_values, "{total:?}");
}
