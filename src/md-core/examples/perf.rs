//! Performance smoke test on a synthetic large run.
//! `cargo run --release -p md-core --example perf -- [scans] [ions_per_scan]`
//! Default: 60 min at 20 scans/s with 500 ions per scan (36 000 scans, 18 M points).

use std::time::Instant;

use md_core::{Chromatogram, decimate};

fn main() {
    let mut args = std::env::args()
        .skip(1)
        .map(|a| a.parse::<usize>().expect("number"));
    let scans = args.next().unwrap_or(36_000);
    let ions = args.next().unwrap_or(500);

    let t = Instant::now();
    let c = Chromatogram {
        scan_times_s: (0..scans).map(|s| s as f64 / 20.0).collect(),
        offsets: (0..=scans).map(|s| s * ions).collect(),
        mz: (0..scans * ions)
            .map(|k| 50.0 + (k % ions) as f32)
            .collect(),
        intensity: (0..scans * ions).map(|k| (k % 997) as f32).collect(),
    };
    let points = c.mz.len();
    println!(
        "{scans} scans x {ions} ions = {points} points, {} MiB, built in {:?}",
        points * 8 / (1 << 20),
        t.elapsed()
    );

    let t = Instant::now();
    let tic = c.tic();
    println!("TIC: {:?}", t.elapsed());

    let t = Instant::now();
    let sic = c.sic(73.0, 0.5);
    println!("SIC m/z 73: {:?}", t.elapsed());

    let t = Instant::now();
    let d = decimate(
        &c.scan_times_s,
        &tic,
        f64::NEG_INFINITY,
        f64::INFINITY,
        2000,
    );
    println!("decimate TIC to {} points: {:?}", d.len(), t.elapsed());
    assert_eq!(sic.len(), scans);
}
