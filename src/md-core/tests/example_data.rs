//! Reads every real run in `example-data/` (fetched by `example-data/fetch.sh` or `fetch.ps1`).
//! Passes without checking anything when the data has not been downloaded.

use std::path::{Path, PathBuf};
use std::time::Instant;

use md_core::deconv::{DeconvSettings, deconvolute};
use md_core::peaks::{IonTraces, PeakSettings};

fn runs() -> Vec<PathBuf> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            dir.display()
        );
        return Vec::new();
    };
    let runs: Vec<_> = entries
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("cdf")))
        .collect();
    assert!(!runs.is_empty(), "no .cdf files in {}", dir.display());
    runs
}

fn name(path: &Path) -> std::borrow::Cow<'_, str> {
    path.file_name().unwrap().to_string_lossy()
}

#[test]
fn reads_all_example_runs() {
    for path in runs() {
        let t = Instant::now();
        let c = md_core::andi::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let tic = c.tic();
        assert!(c.scan_count() > 100, "{}", path.display());
        assert!(tic.iter().any(|&v| v > 0.0), "{}", path.display());
        let mz = |s| c.spectrum(s).0.to_vec();
        assert!((0..c.scan_count()).all(|s| mz(s).is_sorted()));
        eprintln!(
            "{}: {} scans, {} points, {:.1}-{:.1} min, read+TIC {:?}",
            name(&path),
            c.scan_count(),
            c.mz.len(),
            c.scan_times_s[0] / 60.0,
            c.scan_times_s[c.scan_count() - 1] / 60.0,
            t.elapsed()
        );
    }
}

#[test]
fn detects_peaks_in_example_runs() {
    for path in runs() {
        let c = md_core::andi::read(&path).unwrap();
        let t = Instant::now();
        let traces = IonTraces::new(&c).unwrap();
        let index = t.elapsed();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let total = t.elapsed();
        let per_ion = peaks.len() as f64 / traces.len() as f64;
        eprintln!(
            "{}: {} traces, {} peaks ({per_ion:.1} per trace), index {index:?}, index+peaks {total:?}",
            name(&path),
            traces.len(),
            peaks.len(),
        );
        // About 2 per trace on these runs (all clear peaks found, see docs/images); far outside
        // this range means a broken detector, not different data.
        assert!((1.0..200.0).contains(&per_ion), "{}", path.display());
        for p in &peaks {
            assert!(
                p.first_scan < p.apex_scan && p.apex_scan < p.last_scan,
                "{p:?}"
            );
            assert!(p.signal_to_noise >= 20.0, "{p:?}");
        }
    }
}

#[test]
fn deconvolutes_example_runs() {
    for path in runs() {
        let c = md_core::andi::read(&path).unwrap();
        let t = Instant::now();
        let traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let s = DeconvSettings::default();
        let comps = deconvolute(&traces, &peaks, &s);
        let total = t.elapsed();
        let ions: usize = comps.iter().map(|k| k.spectrum.len()).sum();
        eprintln!(
            "{}: {} peaks -> {} components ({:.1} ions each), index+peaks+deconvolution {total:?}",
            name(&path),
            peaks.len(),
            comps.len(),
            ions as f64 / comps.len().max(1) as f64,
        );
        // A derivatised plant extract has dozens of substances; none or thousands means a broken step.
        assert!((10..500).contains(&comps.len()), "{}", path.display());
        for k in &comps {
            // Like the old code, duplicated peaks count towards the minimum, distinct masses may be fewer.
            assert!(
                k.peaks.len() >= s.min_peak_number && !k.spectrum.is_empty(),
                "{k:?}"
            );
            assert!(
                k.first_scan <= k.apex_scan && k.apex_scan <= k.last_scan,
                "{k:?}"
            );
        }
    }
}

/// Tab-separated rows of a TargetSearchData text file, header skipped.
const DENSITY: usize = 40;
const THRESHOLD: i32 = 100;

fn table(path: &Path) -> Vec<Vec<String>> {
    std::fs::read_to_string(path)
        .unwrap()
        .lines()
        .skip(1)
        .map(|l| l.split('\t').map(str::to_owned).collect())
        .collect()
}

#[test]
fn ri_of_markers_matches_target_search() {
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    // FAME markers: time window (s) and RI; TargetSearch finds them on m/z 87. The RI files list
    // TargetSearch's RI for each scan, so the comparison covers interpolation and extrapolation.
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let peaks = md_core::ri::marker_peaks(&c, 87.0, DENSITY, THRESHOLD);
        let points = md_core::ri::window_points(&peaks, &windows);
        assert_eq!(points.len(), windows.len(), "{}: {peaks:?}", name(path));
        let cal = md_core::ri::RiCalibration::new(points.iter().copied()).unwrap();
        let reference = table(&dir.join(format!(
            "RI_{}.txt",
            path.file_stem().unwrap().to_string_lossy()
        )));
        let worst = reference
            .iter()
            .map(|r| {
                let (rt, ri): (f64, f64) = (r[0].parse().unwrap(), r[2].parse().unwrap());
                (cal.ri(rt) - ri).abs()
            })
            .fold(0.0, f64::max);
        eprintln!(
            "{}: markers {points:?}, worst RI difference {worst:.1}",
            name(path)
        );
        // TargetSearch put the middle marker of 7235eg06 at 321.1 s, where m/z 87 reaches 36 512;
        // the marker peak at 311.1 s (as in every other run) reaches 221 839.
        if name(path) == "7235eg06.cdf" {
            assert!((points[1].0 - 311.109).abs() < 1e-6, "{points:?}");
        } else {
            assert!(worst < 0.01, "{}: RI off by {worst}", name(path));
        }
    }
}

/// RI search with the sample library `lib_aa.lbr` (alkane RI; `MD_LBR_SAMPLE` as in `lbr_sample.rs`).
/// The example runs have FAME markers but no alkanes, so their alkane calibration is bridged over two
/// amino acids: the FAME RI in TargetSearch's `library.txt` gives the time in the run, `lib_aa` the RI.
/// Amino acids between them must then be assigned with a high RI score.
#[test]
fn amino_acids_get_ri_score_with_lib_aa() {
    use md_core::ri::{RiCalibration, marker_peaks, window_points};
    use md_core::search::{PreparedLibrary, SearchSettings, map_compounds_to_library};
    let Ok(lbr) = std::env::var("MD_LBR_SAMPLE") else {
        eprintln!("skipped: MD_LBR_SAMPLE not set");
        return;
    };
    let runs = runs();
    let Some(dir) = runs.first().and_then(|p| p.parent()) else {
        return;
    };
    let lib = md_core::library::read_lbr(&lbr).unwrap();
    let lib_ri = |n: &str| lib.iter().find(|c| c.name == n).unwrap().spectrum.ri;
    let windows: Vec<(f64, f64, f64)> = table(&dir.join("rimLimits.txt"))
        .iter()
        .map(|r| {
            (
                r[0].parse().unwrap(),
                r[1].parse().unwrap(),
                r[2].parse().unwrap(),
            )
        })
        .collect();
    // (FAME RI in library.txt, name in lib_aa)
    let bridge = [(271500.0, "Valine 2TMS"), (368800.0, "Threonine 3TMS")];
    let mut best = [
        ("Leucine 2TMS", 0.0),
        ("Glycine 3TMS", 0.0),
        ("Serine 3TMS", 0.0),
    ];
    let s = SearchSettings {
        cutoff: 0.0,
        ..Default::default()
    };
    let prepared = PreparedLibrary::new(&lib, &s);
    for path in &runs {
        let c = md_core::andi::read(path).unwrap();
        let fame = marker_peaks(&c, 87.0, DENSITY, THRESHOLD);
        let fame = RiCalibration::new(window_points(&fame, &windows)).unwrap();
        let cal = RiCalibration::new(bridge.map(|(f, n)| (fame.rt(f), lib_ri(n)))).unwrap();
        let traces = IonTraces::new(&c).unwrap();
        let peaks = traces.detect_peaks(&PeakSettings::default());
        let comps = deconvolute(&traces, &peaks, &DeconvSettings::default());
        let ris = cal.component_ris(&c.scan_times_s, &comps);
        let id = map_compounds_to_library(&comps, Some(&ris), &prepared, &s);
        for h in id.assigned.iter().flatten() {
            // Only clear spectrum matches count; weak ones are other substances.
            let (n, ri) = (&lib[h.compound].name, h.ri_sim.unwrap());
            if h.spec_sim > 0.65
                && let Some(b) = best.iter_mut().find(|b| b.0 == n)
            {
                eprintln!("{}: {n} spectrum {:.3}, RI {ri:.3}", name(path), h.spec_sim);
                b.1 = f64::max(b.1, ri);
            }
        }
    }
    if !runs.is_empty() {
        // Each is found in at least one run, less than 0.46 sigma (9 RI) off its library RI (measured 0.94-0.99).
        assert!(best.iter().all(|b| b.1 > 0.9), "{best:?}");
    }
}
