//! Times the batch quantification on real files.
//! `cargo run --release -p md-core --example batch -- PARALLEL LIBRARY.msl|.lbr RUN.cdf...`

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::time::Instant;

use md_core::batch::{ScanSettings, TableSettings, scan_runs, table};
use md_core::library::{import_msl, read_lbr};
use md_core::search::PreparedLibrary;

fn main() {
    let mut args = std::env::args().skip(1);
    let usage = "PARALLEL LIBRARY RUN...";
    let parallel: usize = args.next().expect(usage).parse().expect(usage);
    let lib_path = PathBuf::from(args.next().expect(usage));
    let paths: Vec<PathBuf> = args.map(PathBuf::from).collect();
    let lib = if lib_path.extension().is_some_and(|e| e == "lbr") {
        read_lbr(&lib_path)
    } else {
        import_msl(&lib_path)
    }
    .expect("library");
    let s = ScanSettings::default();
    let prepared = PreparedLibrary::new(&lib, &s.search);

    let t = Instant::now();
    let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
    let (runs, errors) = scan_runs(&paths, &prepared, &s, parallel, &done, &cancel).unwrap();
    let scanned = t.elapsed();
    let t = Instant::now();
    let tab = table(&runs, &lib, &TableSettings::default()).unwrap();
    println!(
        "{} runs ({} errors) x {} library entries, {parallel} in parallel: scan {scanned:?}, table {:?}",
        runs.len(),
        errors.len(),
        lib.len(),
        t.elapsed()
    );
    let cells = tab.rows.iter().flat_map(|r| &r.values).flatten().count();
    println!(
        "{} substances, {cells} of {} cells filled",
        tab.rows.len(),
        tab.rows.len() * tab.runs.len()
    );
    for e in errors {
        println!("error: {e}");
    }
}
