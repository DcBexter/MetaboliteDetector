//! Reads every real run in `example-data/` (fetched by `example-data/fetch.sh` or `fetch.ps1`).
//! Passes without checking anything when the data has not been downloaded.

use std::path::{Path, PathBuf};
use std::time::Instant;

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
        let traces = IonTraces::new(&c);
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
