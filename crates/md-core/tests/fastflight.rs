//! FastFlight import against a synthetic run. No real FastFlight files are available, so the run is
//! built from the reverse-engineered layout; `matches_old_importer` checks that layout against the
//! old program's own importer (docs/reverse-engineering/tools/fastflight/make-cache.sh):
//!   MD_FF_SYNTH=DIR cargo test -p md-core --test fastflight write_synthetic
//!   docs/reverse-engineering/tools/fastflight/make-cache.sh DIR
//!   MD_FF_OLD=DIR cargo test -p md-core --test fastflight matches_old_importer

use std::path::{Path, PathBuf};

use md_core::fastflight;
use md_core::tof::TofCalibration;

enum Item {
    Base(u32),
    Pos(u16),
    Val(u32),
}
use Item::{Base, Pos, Val};

struct Spectrum {
    n: i32,
    /// FFC time factor; 0 stands for the default 2.5.
    factor: f64,
    /// First frame carries the time word, the last one is final.
    frames: Vec<Vec<Item>>,
    /// TIC in the FFC; the trailer of the FFT gets the real one.
    ffc_tic: Option<u64>,
}

fn words(w: &[u32]) -> Vec<u8> {
    w.iter()
        .flat_map(|w| w.to_le_bytes()[..3].to_vec())
        .collect()
}

impl Spectrum {
    fn intensity_sum(&self) -> u64 {
        let all = self.frames.iter().flatten();
        all.map(|i| if let Val(v) = i { u64::from(*v) } else { 0 })
            .sum()
    }

    fn encode(&self) -> Vec<u8> {
        let mut w = Vec::new();
        let last = self.frames.len() - 1;
        for (i, frame) in self.frames.iter().enumerate() {
            let mut data = Vec::new();
            for item in frame {
                match item {
                    Base(b) => data.extend([0xFF_8000, *b]),
                    Pos(p) => data.push(0xFF_0000 | u32::from(*p)),
                    Val(v) => data.push(*v),
                }
            }
            let (flag, skip) = match i {
                0 => (0x4000, 3),
                i if i == last => (0x8000, 4),
                _ => (0, 2),
            };
            w.extend([0xFF_FFFF, flag | (skip + data.len() as u32)]);
            if i == 0 {
                w.push(777); // acquisition time, unused
            }
            w.extend(data);
        }
        let tic = self.intensity_sum();
        w.extend([(tic & 0xFF_FFFF) as u32, (tic >> 24) as u32]);
        words(&w)
    }
}

/// 60 000 points, as many as the old importer needs before it stores a scan uncompressed.
fn big(n: i32, seed: u32) -> Spectrum {
    let mut x = seed;
    let mut values = || {
        x = x.wrapping_mul(1_664_525).wrapping_add(1_013_904_223);
        Val(3 + (x >> 8) % 4000)
    };
    let mut frames: Vec<Vec<Item>> = (0..4)
        .map(|_| (0..15_000).map(|_| values()).collect())
        .collect();
    frames[0].splice(0..0, [Base(5000), Pos(0)]);
    Spectrum {
        n,
        factor: 2.5,
        frames,
        ffc_tic: None,
    }
}

fn spectra() -> Vec<Spectrum> {
    let s = |n, factor, frames| Spectrum {
        n,
        factor,
        frames,
        ffc_tic: None,
    };
    vec![
        // Two frames, a base offset and a jump of the position inside the first frame.
        s(
            100,
            2.5,
            vec![
                vec![
                    Base(1000),
                    Pos(5),
                    Val(30),
                    Val(31),
                    Val(35),
                    Pos(20),
                    Val(90),
                    Val(12),
                ],
                vec![Pos(40), Val(20), Val(25)],
            ],
        ),
        // Factor 0 means 2.5; a middle frame; a second base offset.
        s(
            101,
            0.0,
            vec![
                vec![Pos(7), Val(8), Val(9)],
                vec![Base(300), Pos(1), Val(4), Val(6)],
                vec![Val(5)],
            ],
        ),
        // Time before the others although stored last; TIC above 24 bit.
        s(
            50,
            4.0,
            vec![
                vec![Pos(2), Val(0xF0_0000), Val(0xF0_0000)],
                vec![Val(0xF0_0000), Val(3)],
            ],
        ),
        // Empty scan: TIC 0 is skipped by the old importer.
        Spectrum {
            n: 103,
            factor: 2.5,
            frames: vec![vec![Pos(1), Val(0)], vec![]],
            ffc_tic: Some(0),
        },
        // The same time as the first scan: the first one stays.
        s(
            100,
            2.5,
            vec![vec![Pos(9), Val(77), Val(78)], vec![Val(79)]],
        ),
        // A single point: nothing is left after the first one is dropped.
        s(104, 2.5, vec![vec![], vec![Pos(3), Val(10)]]),
        big(200, 1),
        big(201, 2),
    ]
}

fn calibration() -> TofCalibration {
    TofCalibration {
        coefficients: vec![10.0, 0.25, 1.5e-5],
    }
}

/// `MsData.FFT`, `MsData.FFC` and `calib.mzcalib` of the synthetic run.
fn write_run(dir: &Path) {
    std::fs::create_dir_all(dir).unwrap();
    let (mut fft, mut ffc) = (Vec::new(), Vec::new());
    ffc.extend((spectra().len() as i32).to_le_bytes());
    for s in spectra() {
        ffc.extend([0u8; 4]);
        ffc.extend(s.n.to_le_bytes());
        ffc.extend([0u8; 4]);
        ffc.extend(s.ffc_tic.unwrap_or(s.intensity_sum()).to_le_bytes());
        ffc.extend((fft.len() as u64).to_le_bytes());
        ffc.extend(s.factor.to_le_bytes());
        fft.extend(s.encode());
    }
    std::fs::write(dir.join("MsData.FFT"), fft).unwrap();
    std::fs::write(dir.join("MsData.FFC"), ffc).unwrap();
    let pairs = [(400.0, 110.0), (1500.0, 407.5)];
    calibration()
        .write(dir.join("calib.mzcalib"), &pairs)
        .unwrap();
}

fn temp_dir(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!("md-core-ff-{}-{name}", std::process::id()))
}

#[test]
fn reads_the_synthetic_run() {
    let dir = temp_dir("read");
    write_run(&dir);
    let cal = TofCalibration::read(dir.join("calib.mzcalib")).unwrap();
    for path in [dir.join("MsData.FFT"), dir.join("MsData.FFC"), dir.clone()] {
        let run = fastflight::read(&path, &cal).unwrap();
        // Sorted by time; the scan without signal and the repeated time are left out.
        assert_eq!(run.scan_times_s, [0.2, 0.25, 0.2525, 0.26, 0.5, 0.5025]);
        // Per scan the first point is dropped and the smallest intensity subtracted.
        let expected: [(&[f32], &[f32]); 4] = [
            (&[3.0, 4.0, 5.0], &[15_728_637.0, 15_728_637.0, 0.0]),
            (
                &[1006.0, 1007.0, 1020.0, 1021.0, 1040.0, 1041.0],
                &[19.0, 23.0, 78.0, 0.0, 8.0, 13.0],
            ),
            (&[8.0, 301.0, 302.0, 303.0], &[5.0, 0.0, 2.0, 1.0]),
            (&[], &[]),
        ];
        for (scan, (bins, values)) in expected.iter().enumerate() {
            let (mz, intensity) = run.spectrum(scan);
            let want: Vec<f32> = bins.iter().map(|b| cal.mz(*b)).collect();
            assert_eq!((mz, intensity), (&want[..], *values), "scan {scan}");
        }
        assert_eq!(run.spectrum(4).0.len(), 59_999);
        assert!(run.spectrum(5).1.iter().all(|v| (0.0..4000.0).contains(v)));
    }
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn damaged_runs_are_errors() {
    let dir = temp_dir("damaged");
    write_run(&dir);
    let cal = calibration();
    let fft = std::fs::read(dir.join("MsData.FFT")).unwrap();
    let ffc = std::fs::read(dir.join("MsData.FFC")).unwrap();
    // A shortened data file, a shortened index and a TIC that differs between them.
    std::fs::write(dir.join("MsData.FFT"), &fft[..fft.len() - 1]).unwrap();
    assert!(fastflight::read(&dir, &cal).is_err());
    std::fs::write(dir.join("MsData.FFT"), &fft).unwrap();
    std::fs::write(dir.join("MsData.FFC"), &ffc[..ffc.len() - 1]).unwrap();
    assert!(fastflight::read(&dir, &cal).is_err());
    let mut bad = ffc.clone();
    bad[4 + 12] ^= 1; // TIC of the first scan
    std::fs::write(dir.join("MsData.FFC"), bad).unwrap();
    let e = fastflight::read(&dir, &cal).unwrap_err().to_string();
    assert!(e.contains("TIC mismatch"), "{e}");
    assert!(fastflight::read(&dir.join("missing"), &cal).is_err());
    // Entries that all point at the first spectrum would multiply its points.
    let mut shared = ffc.clone();
    for entry in shared[4..].chunks_mut(36) {
        entry[20..28].fill(0);
    }
    std::fs::write(dir.join("MsData.FFC"), shared).unwrap();
    let e = fastflight::read(&dir, &cal).unwrap_err().to_string();
    assert!(e.contains("same data") || e.contains("TIC"), "{e}");
    // A negative scan number is refused.
    let mut negative = ffc.clone();
    negative[8..12].copy_from_slice(&(-5i32).to_le_bytes());
    std::fs::write(dir.join("MsData.FFC"), negative).unwrap();
    assert!(fastflight::read(&dir, &cal).is_err());
    std::fs::remove_dir_all(dir).ok();
}

#[test]
fn read_run_finds_the_calibration_next_to_the_run() {
    let dir = temp_dir("read-run");
    write_run(&dir);
    let run = md_core::read_run(&dir.join("MsData.FFT")).unwrap();
    assert_eq!(run.scan_count(), 6);
    std::fs::remove_file(dir.join("calib.mzcalib")).unwrap();
    let e = md_core::read_run(&dir.join("MsData.FFT")).unwrap_err();
    assert!(e.contains(".mzcalib"), "{e}");
    std::fs::remove_dir_all(dir).ok();
}

/// Writes the synthetic run for the old importer (see the top of this file).
#[test]
fn write_synthetic() {
    match std::env::var_os("MD_FF_SYNTH") {
        Some(dir) => write_run(Path::new(&dir)),
        None => eprintln!("skipped: MD_FF_SYNTH not set"),
    }
}

/// The old importer's cache of the same run (`DIR/old/synthetic.bin`) holds the same points. Run with
/// `importData(.., true)` it keeps the scans of more than 50 000 points and skips the others. Its
/// writer repeats points and loses about half a percent of them, so the spectra are compared as a
/// sequence: every point the old cache holds is one of ours, in the same order, and the times agree.
#[test]
fn matches_old_importer() {
    let Some(dir) = std::env::var_os("MD_FF_OLD").map(PathBuf::from) else {
        eprintln!("skipped: MD_FF_OLD not set");
        return;
    };
    let cal = TofCalibration::read(dir.join("calib.mzcalib")).unwrap();
    let ours = fastflight::read(&dir, &cal).unwrap();
    let old = md_core::diskscan::read(dir.join("old/synthetic.bin")).unwrap();
    let big: Vec<usize> = (0..ours.scan_count())
        .filter(|s| ours.spectrum(*s).0.len() > 50_000)
        .collect();
    assert_eq!(old.scan_count(), big.len());
    assert_eq!(big.len(), 2);
    for (o, s) in big.into_iter().enumerate() {
        assert_eq!(old.scan_times_s[o], ours.scan_times_s[s], "scan {s}");
        let (mz, iv) = ours.spectrum(s);
        let mut mine = mz.iter().zip(iv).peekable();
        let (omz, oiv) = old.spectrum(o);
        let mut found = 0;
        for point in omz.iter().zip(oiv) {
            // Skip forward to the same point; a miss means the old one has a point we do not.
            while mine.next_if(|p| *p != point).is_some() {}
            assert!(
                mine.peek().is_some(),
                "scan {s}: {point:?} is not in our spectrum"
            );
            found += 1;
        }
        assert!(found > 100_000, "scan {s}: {found}");
        let distinct = |v: &[f32]| {
            v.iter()
                .map(|x| x.to_bits())
                .collect::<std::collections::BTreeSet<_>>()
                .len()
        };
        assert!(distinct(omz) * 100 > mz.len() * 98, "scan {s}");
    }
}
