//! Runs `md-cli` on the 15 runs in `example-data/` (fetched by `example-data/fetch.*`) with the
//! TargetSearchData library and compares its CSV with the table `md-core` builds in process, which is
//! the one of the app's batch window. Passes without checking anything when the data is missing.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, AtomicUsize};

use md_core::batch::{TableSettings, scan_runs, table};
use md_core::library::parse_msl;
use md_core::search::PreparedLibrary;
use md_core::settings::GcmsSettings;

fn data_dir() -> Option<PathBuf> {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../example-data/TargetSearchData/inst/gc-ms-data");
    if !dir.is_dir() {
        eprintln!(
            "skipped: {} not found, run example-data/fetch.*",
            dir.display()
        );
        return None;
    }
    Some(dir)
}

/// MSL text of TargetSearch's `library.txt` (name, tab-separated columns, last one `mass:intensity ...`).
fn library_msl(dir: &Path) -> String {
    let text = std::fs::read_to_string(dir.join("library.txt")).unwrap();
    text.lines()
        .skip(1)
        .map(|l| {
            let cols: Vec<&str> = l.split('\t').collect();
            let peaks: Vec<String> = cols[4]
                .split(' ')
                .map(|p| p.replace(':', " "))
                .map(|p| format!("({p})"))
                .collect();
            format!(
                "NAME: {}\nNUM PEAKS: {}\n{}\n\n",
                cols[0],
                peaks.len(),
                peaks.join(" ")
            )
        })
        .collect()
}

#[test]
fn csv_equals_the_core_table_of_the_example_runs() {
    let Some(dir) = data_dir() else { return };
    let mut runs: Vec<PathBuf> = std::fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e.eq_ignore_ascii_case("cdf")))
        .collect();
    runs.sort();
    assert_eq!(runs.len(), 15);

    let msl = library_msl(&dir);
    let lib_path = Path::new(env!("CARGO_TARGET_TMPDIR")).join("example-library.msl");
    std::fs::write(&lib_path, &msl).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_md-cli"))
        .arg("--library")
        .arg(&lib_path)
        .args(&runs)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );

    let settings = GcmsSettings::default();
    let compounds = parse_msl(&msl).unwrap();
    let prepared = PreparedLibrary::new(&compounds, &settings.search());
    let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
    let (found, errors) =
        scan_runs(&runs, &prepared, &settings.scan(None), 4, &done, &cancel).unwrap();
    assert!(errors.is_empty(), "{errors:?}");
    let t = table(&found, &compounds, &TableSettings::default()).unwrap();
    let mut expected = Vec::new();
    t.write_csv(&mut expected).unwrap();

    assert_eq!(
        String::from_utf8(out.stdout).unwrap(),
        String::from_utf8(expected).unwrap()
    );
    assert_eq!(t.runs.len(), 15);
    assert!(!t.rows.is_empty(), "no substance found in any run");
}
