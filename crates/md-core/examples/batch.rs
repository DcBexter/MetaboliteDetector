//! Times the batch quantification on real files.
//! `cargo run --release -p md-core --example batch -- [--untargeted] [--rescan] PARALLEL LIBRARY.msl|.lbr RUN.cdf...`

use std::path::PathBuf;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::time::Instant;

use md_core::batch::{
    QuantSettings, ScanSettings, TableSettings, add_unknowns, rescan_runs, scan_runs, table,
};
use md_core::library::{import_msl, read_lbr};
use md_core::search::PreparedLibrary;

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let mut flag = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .map(|i| args.remove(i))
            .is_some()
    };
    let (untargeted, rescan) = (flag("--untargeted"), flag("--rescan"));
    let mut args = args.into_iter();
    let usage = "[--untargeted] [--rescan] PARALLEL LIBRARY RUN...";
    let parallel: usize = args.next().expect(usage).parse().expect(usage);
    let lib_path = PathBuf::from(args.next().expect(usage));
    let paths: Vec<PathBuf> = args.map(PathBuf::from).collect();
    let lib = if lib_path.extension().is_some_and(|e| e == "lbr") {
        read_lbr(&lib_path)
    } else {
        import_msl(&lib_path)
    }
    .expect("library");
    let s = ScanSettings {
        untargeted,
        ..Default::default()
    };
    let prepared = PreparedLibrary::new(&lib, &s.search);

    let t = Instant::now();
    let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
    let (mut runs, errors) = scan_runs(&paths, &prepared, &s, parallel, &done, &cancel).unwrap();
    let scanned = t.elapsed();
    let t = Instant::now();
    let extra = if untargeted {
        add_unknowns(&mut runs, lib.len(), &s, &cancel).unwrap()
    } else {
        Vec::new()
    };
    let pooled = t.elapsed();
    let t = Instant::now();
    let rescanned = rescan.then(|| {
        rescan_runs(
            &mut runs,
            &lib,
            &extra,
            &s,
            &QuantSettings::default(),
            parallel,
            &cancel,
        )
        .unwrap()
    });
    let rescan_time = t.elapsed();
    let t = Instant::now();
    let tab = table(&runs, &lib, &extra, &TableSettings::default()).unwrap();
    println!(
        "{} runs ({} errors) x {} library entries, {parallel} in parallel: scan {scanned:?}, pool {pooled:?} ({} unknown), rescan {rescan_time:?} ({:?} ions added), table {:?}",
        runs.len(),
        errors.len(),
        lib.len(),
        extra.len(),
        rescanned.as_ref().map(|r| r.0),
        t.elapsed()
    );
    let cells = tab.rows.iter().flat_map(|r| &r.values).flatten().count();
    println!(
        "{} substances, {cells} of {} cells filled",
        tab.rows.len(),
        tab.rows.len() * tab.runs.len()
    );
    for e in errors.iter().chain(rescanned.iter().flat_map(|r| &r.1)) {
        println!("error: {e}");
    }
}
