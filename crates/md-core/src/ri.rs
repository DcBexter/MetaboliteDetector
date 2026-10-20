//! Retention index (RI): piecewise linear RT ↔ RI calibration (old `gcms::RICalculator`), marker
//! search in a calibration run, RT shift correction and the `.calib` text file.
//! Algorithms: `docs/reverse-engineering/algorithms.md` §6; file: `file-formats.md` §4.1.

use std::io::Read;
use std::path::Path;

use crate::Chromatogram;
use crate::deconv::Component;
use crate::peaks::detect_peaks_matched;

#[derive(Debug, thiserror::Error)]
pub enum CalibError {
    #[error("cannot read calibration file: {0}")]
    Io(#[from] std::io::Error),
    #[error("calibration file line {line}: {reason}")]
    Line { line: usize, reason: String },
    #[error("invalid calibration: {0}")]
    Invalid(String),
}

/// RT → RI calibration: at least two points, retention time and RI both strictly ascending.
#[derive(Debug, Clone, PartialEq)]
pub struct RiCalibration {
    /// (retention time in s, RI).
    points: Vec<(f64, f64)>,
    /// Constant retention time shift `drt` in s (§6.5), subtracted before the lookup.
    pub shift_s: f64,
}

impl RiCalibration {
    /// Like the old `addCalibrationPair` (std::map insert), the first point of a repeated retention
    /// time wins. Unlike the old code, RI must rise with retention time and two points are required:
    /// otherwise RT → RI has no inverse and the old code returned -1 for every RI.
    pub fn new(points: impl IntoIterator<Item = (f64, f64)>) -> Result<Self, CalibError> {
        let mut points: Vec<_> = points.into_iter().collect();
        if points
            .iter()
            .any(|(t, ri)| !t.is_finite() || !ri.is_finite())
        {
            return Err(CalibError::Invalid("not a finite number".into()));
        }
        points.sort_by(|a, b| a.0.total_cmp(&b.0));
        points.dedup_by(|b, a| a.0 == b.0);
        if points.len() < 2 {
            return Err(CalibError::Invalid(format!(
                "{} point(s), at least 2 needed",
                points.len()
            )));
        }
        if let Some(w) = points.windows(2).find(|w| w[1].1 <= w[0].1) {
            return Err(CalibError::Invalid(format!(
                "RI {} at {} s is not above RI {} at {} s",
                w[1].1, w[1].0, w[0].1, w[0].0
            )));
        }
        // Finite numbers can still give slopes that overflow on extrapolation (huge RI steps,
        // near-equal times). Real series have about 1 to 1000 RI units per second.
        if let Some(w) = points.windows(2).find(|w| {
            let slope = (w[1].1 - w[0].1) / (w[1].0 - w[0].0);
            !(1e-9..=1e9).contains(&slope)
        }) {
            return Err(CalibError::Invalid(format!(
                "slope between {} s and {} s is outside 1e-9 to 1e9 RI per second",
                w[0].0, w[1].0
            )));
        }
        Ok(Self {
            points,
            shift_s: 0.0,
        })
    }

    /// (retention time in s, RI), ascending; without the shift.
    pub fn points(&self) -> &[(f64, f64)] {
        &self.points
    }

    /// Old `calculateRI(rt)` (§6.2).
    pub fn ri(&self, rt_s: f64) -> f64 {
        interpolate(&self.points, rt_s - self.shift_s, |p| p)
    }

    /// Old `calculateRT(ri)` (§6.3), the inverse of [`Self::ri`].
    pub fn rt(&self, ri: f64) -> f64 {
        interpolate(&self.points, ri, |(t, ri)| (ri, t)) + self.shift_s
    }

    /// Old `calculateRI(compounds)` (§6.4): RI of each component at its apex scan.
    pub fn component_ris(&self, scan_times_s: &[f64], components: &[Component]) -> Vec<f64> {
        components
            .iter()
            .map(|c| self.ri(scan_times_s[c.apex_scan as usize]))
            .collect()
    }

    /// RT shift correction (§6.5): the calibration shifted so that the reference substance, found
    /// at `found_rt_s`, gets its library RI. The old code subtracts the previous shift once more
    /// (correct only from shift 0); here the result does not depend on the previous shift.
    pub fn shifted_to(&self, found_rt_s: f64, reference_ri: f64) -> Self {
        Self {
            shift_s: found_rt_s - interpolate(&self.points, reference_ri, |(t, ri)| (ri, t)),
            ..self.clone()
        }
    }

    /// Reads a `.calib` file. Its content comes from users, so it is checked line by line.
    pub fn read_calib(path: impl AsRef<Path>) -> Result<Self, CalibError> {
        // A calibration has tens of lines; the cap stops a wrong file from filling memory.
        const MAX_BYTES: u64 = 1 << 20;
        let mut bytes = Vec::new();
        std::fs::File::open(path)?
            .take(MAX_BYTES + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_BYTES {
            return Err(CalibError::Invalid(format!(
                "file larger than {MAX_BYTES} bytes"
            )));
        }
        let text = std::str::from_utf8(&bytes)
            .map_err(|_| CalibError::Invalid("not a text file".into()))?;
        Self::parse_calib(text)
    }

    /// Parses `.calib` text: one `<rt in ms> <RI>` pair per line, separated by whitespace
    /// (the old `fromDisk` splits on `\s+`). Blank lines are skipped.
    pub fn parse_calib(text: &str) -> Result<Self, CalibError> {
        let mut points = Vec::new();
        for (i, l) in text.lines().enumerate() {
            let err = |reason: &str| CalibError::Line {
                line: i + 1,
                reason: reason.into(),
            };
            let fields: Vec<_> = l.split_whitespace().collect();
            match fields[..] {
                [] => {}
                [t, ri] => {
                    let num = |s: &str| s.parse::<f64>().ok().filter(|v| v.is_finite());
                    let (Some(t), Some(ri)) = (num(t), num(ri)) else {
                        return Err(err("not a number"));
                    };
                    points.push((t / 1000.0, ri));
                }
                _ => return Err(err("expected two numbers: retention time (ms) and RI")),
            }
        }
        Self::new(points)
    }

    /// `.calib` text as the old `toDisk` writes it (`<rt in ms>` TAB `<RI>`), but with full
    /// precision instead of 6 significant digits. Like the old file, the shift is not stored.
    pub fn to_calib(&self) -> String {
        self.points
            .iter()
            .map(|(t, ri)| format!("{}\t{ri}\n", t * 1000.0))
            .collect()
    }
}

/// Piecewise linear `y(x)` through `points` (ascending in x after `xy`), extrapolated with the first
/// or last segment. The old code dereferenced `end()` above the last point (§6.2); here that case
/// uses the last segment, as the old code does for every other point above the range.
fn interpolate(points: &[(f64, f64)], x: f64, xy: fn((f64, f64)) -> (f64, f64)) -> f64 {
    let at = |i: usize| xy(points[i]);
    let i = points.partition_point(|p| xy(*p).0 < x);
    if i < points.len() && at(i).0 == x {
        return at(i).1;
    }
    let i = i.clamp(1, points.len() - 1);
    let ((x1, y1), (x2, y2)) = (at(i - 1), at(i));
    y1 + (x - x1) * (y2 - y1) / (x2 - x1)
}

/// Peaks of the marker ion `mz` (± 0.5) in a calibration run, found like the old calibration dialog
/// with the matched peak detector (§3.3; `density` and `threshold` are the dialog's settings).
/// Returns (time in s, intensity) of each peak's highest scan, ascending in time.
pub fn marker_peaks(c: &Chromatogram, mz: f32, density: usize, threshold: i32) -> Vec<(f64, f64)> {
    let sic = c.sic(mz, 0.5);
    let trace: Vec<f32> = sic.iter().map(|v| *v as f32).collect();
    detect_peaks_matched(&trace, density, threshold, 1)
        .iter()
        .map(|p| {
            // Highest raw scan of the peak: the filtered apex can sit a few scans off.
            let apex = (p.first_scan as usize..=p.last_scan as usize)
                .max_by(|a, b| sic[*a].total_cmp(&sic[*b]))
                .unwrap();
            (c.scan_times_s[apex], sic[apex])
        })
        .collect()
}

/// Calibration points for a homologous series (alkanes C10, C11, … → RI 1000, 1100, …): the
/// `ris.len()` highest peaks, in time order, get `ris` in order. Fewer peaks give fewer points.
pub fn series_points(peaks: &[(f64, f64)], ris: &[f64]) -> Vec<(f64, f64)> {
    let mut top = peaks.to_vec();
    top.sort_by(|a, b| b.1.total_cmp(&a.1));
    top.truncate(ris.len());
    top.sort_by(|a, b| a.0.total_cmp(&b.0));
    top.iter().zip(ris).map(|(p, ri)| (p.0, *ri)).collect()
}

/// Calibration points for markers with known time windows `(from s, to s, RI)`: the highest peak in
/// each window. Windows without a peak give no point.
pub fn window_points(peaks: &[(f64, f64)], windows: &[(f64, f64, f64)]) -> Vec<(f64, f64)> {
    windows
        .iter()
        .filter_map(|&(from, to, ri)| {
            peaks
                .iter()
                .filter(|p| (from..=to).contains(&p.0))
                .max_by(|a, b| a.1.total_cmp(&b.1))
                .map(|p| (p.0, ri))
        })
        .collect()
}

/// Reference search of the RT shift correction (§6.5): among the substances whose RI lies within
/// `ri_window` of the reference RI, the one with the highest library score for the reference
/// (at least `cutoff` and above 0; the first wins a tie). `ri_window` is `3 * LS_RI_DIFF` in the old app.
pub fn pick_reference(
    ris: &[f64],
    reference_ri: f64,
    ri_window: f64,
    cutoff: f64,
    score: impl Fn(usize) -> f64,
) -> Option<usize> {
    let mut best = None;
    let mut best_score = 0.0;
    for (i, ri) in ris.iter().enumerate() {
        if (ri - reference_ri).abs() < ri_window {
            let s = score(i);
            if s >= cutoff && s > best_score {
                (best, best_score) = (Some(i), s);
            }
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::*;

    fn alkanes() -> RiCalibration {
        RiCalibration::new([(100.0, 1000.0), (200.0, 1100.0), (400.0, 1200.0)]).unwrap()
    }

    #[test]
    fn interpolates_and_extrapolates_both_ways() {
        let c = alkanes();
        assert_eq!(c.ri(100.0), 1000.0);
        assert_eq!(c.ri(150.0), 1050.0);
        assert_eq!(c.ri(300.0), 1150.0);
        assert_eq!(c.ri(400.0), 1200.0);
        // Outside the range the first or last segment continues.
        assert_eq!(c.ri(50.0), 950.0);
        assert_eq!(c.ri(600.0), 1300.0);
        for t in [-10.0, 50.0, 100.0, 123.0, 400.0, 999.0] {
            assert!((c.rt(c.ri(t)) - t).abs() < 1e-9);
        }
    }

    #[test]
    fn rejects_too_few_or_broken_points() {
        assert!(RiCalibration::new([]).is_err());
        assert!(RiCalibration::new([(1.0, 1000.0)]).is_err());
        // Repeated time: the first point wins, leaving one.
        assert!(RiCalibration::new([(1.0, 1000.0), (1.0, 1100.0)]).is_err());
        assert!(RiCalibration::new([(1.0, 1000.0), (2.0, 1000.0)]).is_err());
        assert!(RiCalibration::new([(1.0, 1000.0), (2.0, 900.0)]).is_err());
        assert!(RiCalibration::new([(1.0, 1000.0), (f64::NAN, 1100.0)]).is_err());
        // Finite points whose slope overflows.
        assert!(RiCalibration::new([(0.0, -1e308), (1.0, 1e308)]).is_err());
        assert!(RiCalibration::new([(0.0, 1000.0), (1e-300, 1100.0)]).is_err());
        // Order of input does not matter.
        let c = RiCalibration::new([(400.0, 1200.0), (100.0, 1000.0), (200.0, 1100.0)]).unwrap();
        assert_eq!(c, alkanes());
    }

    #[test]
    fn shift_moves_reference_onto_its_ri() {
        let c = alkanes().shifted_to(160.0, 1050.0);
        assert_eq!(c.shift_s, 10.0);
        assert_eq!(c.ri(160.0), 1050.0);
        assert_eq!(c.rt(1050.0), 160.0);
        // A second correction replaces the shift instead of stacking on it.
        assert_eq!(c.shifted_to(150.0, 1050.0).shift_s, 0.0);
    }

    #[test]
    fn calib_round_trip_and_broken_files() {
        let c = alkanes();
        let text = c.to_calib();
        assert_eq!(text, "100000\t1000\n200000\t1100\n400000\t1200\n");
        assert_eq!(RiCalibration::parse_calib(&text).unwrap(), c);
        // Old files: any whitespace, 6 significant digits in exponent form, blank lines.
        let old = "1e+05 1000\n\n  2e+05\t\t1100  \r\n4e+05 1200";
        assert_eq!(RiCalibration::parse_calib(old).unwrap(), c);
        for bad in [
            "",
            "1 2\n",
            "1 2 3\n4 5\n",
            "abc 1000\n2 1100\n",
            "1 1000\n2 nan\n",
            "1 1000\n2 inf\n",
        ] {
            assert!(RiCalibration::parse_calib(bad).is_err(), "{bad:?}");
        }
        match RiCalibration::parse_calib("1 1000\nx 1100\n") {
            Err(CalibError::Line { line: 2, .. }) => {}
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn read_calib_rejects_binary_and_huge_files() {
        let dir = std::env::temp_dir().join(format!("md-ri-test-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let (bin, huge) = (dir.join("bin.calib"), dir.join("huge.calib"));
        std::fs::write(&bin, [0xff, 0xfe, 0, 1]).unwrap();
        std::fs::write(&huge, "1 1000\n".repeat(200_000)).unwrap();
        assert!(RiCalibration::read_calib(&bin).is_err());
        assert!(RiCalibration::read_calib(&huge).is_err());
        assert!(RiCalibration::read_calib(dir.join("missing.calib")).is_err());
        std::fs::remove_dir_all(dir).ok();
    }

    #[test]
    fn markers_from_peaks() {
        let peaks = [
            (10.0, 5.0),
            (20.0, 100.0),
            (30.0, 3.0),
            (40.0, 80.0),
            (50.0, 90.0),
        ];
        assert_eq!(
            series_points(&peaks, &[1000.0, 1100.0, 1200.0]),
            [(20.0, 1000.0), (40.0, 1100.0), (50.0, 1200.0)]
        );
        assert_eq!(
            series_points(&peaks[..1], &[1000.0, 1100.0]),
            [(10.0, 1000.0)]
        );
        assert_eq!(
            window_points(
                &peaks,
                &[
                    (5.0, 35.0, 1000.0),
                    (60.0, 70.0, 1100.0),
                    (35.0, 55.0, 1200.0)
                ]
            ),
            [(20.0, 1000.0), (50.0, 1200.0)]
        );
    }

    #[test]
    fn marker_peaks_finds_gaussians() {
        let times: Vec<f64> = (0..400).map(|i| i as f64 * 0.5).collect();
        let gauss = |i: f64, at: f64| 1000.0 * (-0.5 * ((i - at) / 3.0f64).powi(2)).exp();
        let c = Chromatogram {
            offsets: (0..=400).collect(),
            mz: vec![57.0; 400],
            intensity: (0..400)
                .map(|i| (gauss(i as f64, 100.0) + gauss(i as f64, 300.0)) as f32)
                .collect(),
            scan_times_s: times,
        };
        let peaks = marker_peaks(&c, 57.0, 12, 10);
        assert_eq!(peaks.iter().map(|p| p.0).collect::<Vec<_>>(), [50.0, 150.0]);
    }

    #[test]
    fn reference_pick_follows_old_rules() {
        let ris = [1000.0, 1040.0, 1060.0, 1200.0];
        let scores = [0.9, 0.5, 0.8, 0.99];
        let pick = |cutoff| pick_reference(&ris, 1050.0, 30.0, cutoff, |i| scores[i]);
        assert_eq!(pick(0.0), Some(2));
        assert_eq!(pick(0.85), None);
        assert_eq!(pick_reference(&ris, 1050.0, 30.0, 0.0, |_| 0.0), None);
    }
}
