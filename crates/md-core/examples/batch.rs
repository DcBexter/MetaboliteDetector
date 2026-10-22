//! Times the batch quantification on real files.
//! `cargo run --release -p md-core --example batch -- [--untargeted] [--rescan] PARALLEL LIBRARY.msl|.lbr RUN.cdf...`

use md_core::job::{Cancel, Progress};
use std::path::PathBuf;
use std::time::Instant;

use md_core::batch::{ScanSettings, Study, TableSettings, run_study, table};
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
    let prepared = PreparedLibrary::new(&lib, &s.analysis.search);
    let study = Study {
        runs: &paths,
        library: &lib,
        prepared: &prepared,
        scan: &s,
        rescan,
        parallel,
    };

    let t = Instant::now();
    let r = run_study(&study, &Progress::default(), &Cancel::default()).unwrap();
    let took = t.elapsed();
    let t = Instant::now();
    let tab = table(&r, &lib, &TableSettings::default(), None).unwrap();
    println!(
        "{} runs ({} errors) x {} library entries, {parallel} in parallel: study {took:?} ({} unknown, {} ions added by the rescan), table {:?}",
        r.runs.len(),
        r.errors.len(),
        lib.len(),
        r.extra.len(),
        r.rescanned,
        t.elapsed()
    );
    let cells = tab.rows.iter().flat_map(|r| &r.values).flatten().count();
    println!(
        "{} substances, {cells} of {} cells filled",
        tab.rows.len(),
        tab.rows.len() * tab.runs.len()
    );
    for e in &r.errors {
        println!("error: {e}");
    }
}
