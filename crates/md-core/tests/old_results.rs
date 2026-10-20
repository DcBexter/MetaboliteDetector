//! `.cmp` and `.calib` files written by the old MetaboliteDetector 3.3 itself
//! (`MD_RESULTS=example-data/old-results docs/reverse-engineering/tools/golden/make-reference.sh`,
//! never committed). Checks every value against the old program's own dump (`raw.tsv`) and that
//! writing the files again gives the old bytes. Passes without checking anything when they are missing.
//! It also saves its own files to `old-results/md-core/` so that `make-reference.sh --reread` can
//! check that the old program reads them.

use std::path::{Path, PathBuf};

use md_core::cmp::{cmp_payload, read_cmp, write_cmp};
use md_core::ri::{RiCalibration, StoredCalibration};

fn dir() -> Option<PathBuf> {
    let d = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../example-data/old-results");
    d.join("raw.tsv").exists().then_some(d)
}

fn ours(d: &Path) -> PathBuf {
    let o = d.join("md-core");
    std::fs::create_dir_all(&o).unwrap();
    o
}

/// Components of one run in the dump: (scan, rt, [(ion, height)], [peak fields]).
type RefComp = (u32, f64, Vec<(f32, i32)>, Vec<Vec<f64>>);

fn dump(d: &Path) -> Vec<(String, Vec<RefComp>)> {
    let mut runs: Vec<(String, Vec<RefComp>)> = Vec::new();
    for line in std::fs::read_to_string(d.join("raw.tsv")).unwrap().lines() {
        let f: Vec<&str> = line.split('\t').collect();
        let comps = runs.last_mut().map(|r| &mut r.1);
        match f[0] {
            "run" => runs.push((f[1].trim_end_matches(".cdf").into(), Vec::new())),
            "component" => {
                comps
                    .unwrap()
                    .push((f[1].parse().unwrap(), f[2].parse().unwrap(), vec![], vec![]))
            }
            "ion" => (comps.unwrap().last_mut().unwrap().2)
                .push((f[1].parse().unwrap(), f[2].parse().unwrap())),
            "cpeak" => (comps.unwrap().last_mut().unwrap().3)
                .push(f[1..].iter().map(|x| x.parse().unwrap()).collect()),
            _ => {}
        }
    }
    runs
}

/// The dump prints 9 significant digits.
fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-8 * a.abs().max(b.abs()).max(1e-300)
}

#[test]
fn reads_and_rewrites_old_cmp_files() {
    let Some(d) = dir() else {
        eprintln!("skipped: example-data/old-results/ missing");
        return;
    };
    let runs = dump(&d);
    assert!(!runs.is_empty());
    for (name, comps) in &runs {
        let path = d.join(format!("{name}.cmp"));
        let cmp = read_cmp(&path).unwrap();
        assert_eq!(cmp.len(), comps.len(), "{name}");
        for (c, (scan, rt, ions, peaks)) in cmp.iter().zip(comps) {
            let s = &c.spectrum;
            assert_eq!((s.scan, s.ri), (*scan, -1.0), "{name}");
            assert!(
                close(s.rt_ms, *rt),
                "{name} scan {scan}: rt {} vs {rt}",
                s.rt_ms
            );
            let spectrum: Vec<_> = s
                .masses
                .iter()
                .copied()
                .zip(s.values.iter().copied())
                .collect();
            assert_eq!(&spectrum, ions, "{name} scan {scan}");
            assert_eq!(c.peaks.len(), peaks.len(), "{name} scan {scan}");
            for (p, r) in c.peaks.iter().zip(peaks) {
                // old_md_dump.cpp: ion first apex last interp sn quality sharpness integral model
                let got = [
                    p.ion as f64,
                    p.first_scan as f64,
                    p.max_scan as f64,
                    p.last_scan as f64,
                    p.interpolated_max_scan,
                    p.signal_to_noise,
                    p.quality,
                    p.sharpness,
                    p.integral,
                    p.model_parameter,
                ];
                for (i, (g, r)) in got.iter().zip(r).enumerate() {
                    assert!(
                        close(*g, *r),
                        "{name} scan {scan} ion {} field {i}: {g} vs {r}",
                        p.ion
                    );
                }
            }
            let k = c.component().unwrap();
            assert!(
                k.first_scan <= k.apex_scan && k.apex_scan <= k.last_scan,
                "{name}"
            );
            assert_eq!(
                k.model.len() as u32,
                k.last_scan - k.first_scan + 1,
                "{name}"
            );
        }
        // Round trip: our writer gives the old payload byte for byte.
        let mut old = Vec::new();
        use std::io::Read;
        flate2::read::GzDecoder::new(std::fs::File::open(&path).unwrap())
            .read_to_end(&mut old)
            .unwrap();
        assert!(cmp_payload(&cmp) == old, "{name}: payload differs");
        // The run's raw data, when it is here, accepts its own .cmp (scan times match).
        let cdf = (d.join("../TargetSearchData/inst/gc-ms-data")).join(format!("{name}.cdf"));
        if cdf.exists() {
            let times = md_core::read_run(&cdf).unwrap().scan_times_s;
            let k = md_core::cmp::saved_components(&d.join(format!("{name}.cdf")), &times);
            assert_eq!(k.unwrap().unwrap().0.len(), cmp.len(), "{name}");
        }
        let mine = ours(&d).join(format!("{name}.cmp"));
        write_cmp(&mine, &cmp).unwrap();
        assert_eq!(read_cmp(&mine).unwrap(), cmp);
    }
}

#[test]
fn reads_and_rewrites_old_calibration() {
    let Some(d) = dir() else {
        eprintln!("skipped: example-data/old-results/ missing");
        return;
    };
    // old_md_dump.cpp write_calibration: these pairs, source "harness".
    let pairs = [
        (300000.0, 1000.0),
        (612345.678901, 1500.5),
        (901234.5, 2000.25),
        (1.5e6, 3000.0),
    ];
    let bin = std::fs::read(d.join("harness.calib.bin")).unwrap();
    let (stored, len) = StoredCalibration::parse(&bin).unwrap();
    assert_eq!(len, bin.len());
    assert_eq!(stored.points_ms, pairs);
    assert_eq!((stored.drt_ms, &stored.source[..]), (0.0, &b"harness"[..]));
    assert_eq!(stored.to_bytes(), bin);
    std::fs::write(ours(&d).join("harness.calib.bin"), stored.to_bytes()).unwrap();

    // The text file has 6 significant digits.
    let text = RiCalibration::read_calib(d.join("harness.calib")).unwrap();
    for ((t, ri), (t_ms, ri_old)) in text.points().iter().zip(pairs) {
        assert!(
            (t * 1000.0 - t_ms).abs() <= 1e-6 * t_ms && *ri == ri_old,
            "{t} {ri}"
        );
    }
    let again = RiCalibration::parse_calib(&text.to_calib()).unwrap();
    assert_eq!(again, text);
    std::fs::write(ours(&d).join("harness.calib"), text.to_calib()).unwrap();
}
