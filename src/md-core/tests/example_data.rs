//! Reads every real run in `example-data/` (fetched by `example-data/fetch.sh` or `fetch.ps1`).
//! Passes without checking anything when the data has not been downloaded.

use std::path::Path;
use std::time::Instant;

#[test]
fn reads_all_example_runs() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    let Ok(entries) = std::fs::read_dir(&dir) else {
        eprintln!("skipped: {} not found, run example-data/fetch.*", dir.display());
        return;
    };
    let mut runs = 0;
    for path in entries.map(|e| e.unwrap().path()) {
        if path.extension().is_none_or(|e| !e.eq_ignore_ascii_case("cdf")) {
            continue;
        }
        let t = Instant::now();
        let c = md_core::andi::read(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        let tic = c.tic();
        assert!(c.scan_count() > 100, "{}", path.display());
        assert!(tic.iter().any(|&v| v > 0.0), "{}", path.display());
        let mz = |s| c.spectrum(s).0.to_vec();
        assert!((0..c.scan_count()).all(|s| mz(s).is_sorted()));
        eprintln!(
            "{}: {} scans, {} points, {:.1}-{:.1} min, read+TIC {:?}",
            path.file_name().unwrap().to_string_lossy(),
            c.scan_count(),
            c.mz.len(),
            c.scan_times_s[0] / 60.0,
            c.scan_times_s[c.scan_count() - 1] / 60.0,
            t.elapsed()
        );
        runs += 1;
    }
    assert!(runs > 0, "no .cdf files in {}", dir.display());
}
