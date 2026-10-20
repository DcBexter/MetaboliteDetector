//! Centroiding and nominal binning against the old program's own functions (a LD_PRELOAD driver,
//! docs/reverse-engineering/tools/spectrum/run.sh). Without MD_SPEC_* the oracle tests do nothing.

use md_core::centroid::{CentroidSettings, centroid, nominal};
use std::fmt::Write as _;
use std::path::PathBuf;

/// Deterministic generator (xorshift) so the cases are the same on every run.
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// A profile spectrum: Gaussian peaks of random height and width plus integer noise.
fn spectrum(rng: &mut Rng) -> (Vec<f32>, Vec<i32>) {
    let n = [3, 5, 6, 12, 25, 26, 40, 80, 200, 600][rng.below(10) as usize];
    let step = 0.02 + 0.2 * rng.unit();
    let mut mz = Vec::with_capacity(n);
    let mut m = 30.0 + 100.0 * rng.unit();
    for _ in 0..n {
        mz.push(m as f32);
        m += step * (0.5 + rng.unit());
    }
    let peaks: Vec<(f64, f64, f64)> = (0..1 + rng.below(6))
        .map(|_| {
            (
                n as f64 * rng.unit(),
                20.0 + 5000.0 * rng.unit().powi(2),
                1.0 + 6.0 * rng.unit(),
            )
        })
        .collect();
    let noise = [0.0, 2.0, 15.0][rng.below(3) as usize];
    let v = (0..n)
        .map(|i| {
            let signal: f64 = peaks
                .iter()
                .map(|(c, h, w)| h * (-((i as f64 - c) / w).powi(2) / 2.0).exp())
                .sum();
            (signal + noise * rng.unit() + 3.0) as i32
        })
        .collect();
    (mz, v)
}

fn cases() -> String {
    let mut rng = Rng(0x9E37_79B9_7F4A_7C15);
    let mut s = String::new();
    for k in 0..1000 {
        let (mz, v) = spectrum(&mut rng);
        if k % 4 == 3 {
            let _ = writeln!(s, "B {} 1", v.len());
        } else {
            let _ = writeln!(s, "C {} 10 -5 30", v.len());
        }
        for (m, y) in mz.iter().zip(&v) {
            let _ = writeln!(s, "{m:.9} {y}");
        }
    }
    s
}

fn dir(var: &str) -> Option<PathBuf> {
    std::env::var_os(var).map(PathBuf::from)
}

#[test]
fn write_cases() {
    if let Some(d) = dir("MD_SPEC_CASES") {
        std::fs::write(d.join("cases.txt"), cases()).unwrap();
    }
}

#[test]
fn matches_old_program() {
    let Some(d) = dir("MD_SPEC_OLD") else { return };
    let old = std::fs::read_to_string(d.join("old.txt")).unwrap();
    let mut old = old.lines();
    let cases = cases();
    let mut input = cases.lines();
    let mut compared = 0;
    let mut peaks = 0;
    while let Some(head) = input.next() {
        let mut f = head.split(' ');
        let kind = f.next().unwrap();
        let n: usize = f.next().unwrap().parse().unwrap();
        let (mut mz, mut v) = (Vec::new(), Vec::new());
        for _ in 0..n {
            let l = input.next().unwrap();
            let (m, y) = l.split_once(' ').unwrap();
            mz.push(m.parse::<f64>().unwrap() as f32);
            v.push(y.parse::<i32>().unwrap());
        }
        let expected_head = old.next().unwrap();
        let count: usize = expected_head[2..].parse().unwrap();
        let rows: Vec<&str> = (0..count).map(|_| old.next().unwrap()).collect();
        let label = format!("case {compared} ({kind}, {n} points)");
        let parse = |row: &str| -> (f32, i32, f64) {
            let mut c = row.split(' ');
            let m = f32::from_bits(u32::from_str_radix(c.next().unwrap(), 16).unwrap());
            let y: i32 = c.next().unwrap().parse().unwrap();
            let w = c.next().map_or(0.0, |w| w.parse().unwrap());
            (m, y, w)
        };
        let expected: Vec<_> = rows.iter().map(|r| parse(r)).collect();
        let got: Vec<(f32, i32, f64)> = if kind == "B" {
            let (m, y) = nominal(&mz, &v, 1.0);
            m.into_iter().zip(y).map(|(m, y)| (m, y, 0.0)).collect()
        } else {
            let c = centroid(&mz, &v, &CentroidSettings::default());
            peaks += c.mz.len();
            (0..c.mz.len())
                .map(|k| (c.mz[k], c.intensity[k], c.fwhm[k]))
                .collect()
        };
        assert_eq!(got, expected, "{label}");
        compared += 1;
    }
    assert!(old.next().is_none());
    eprintln!("{compared} spectra match, {peaks} centroids");
}
