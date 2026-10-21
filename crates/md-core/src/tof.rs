//! TOF m/z calibration of FastFlight raw data: a polynomial of the TOF bin, stored by the old
//! TOF calibration window in `.mzcalib` text files. Layout: docs/reverse-engineering/file-formats.md §4.2.

use std::path::{Path, PathBuf};

use nalgebra::{DMatrix, DVector};

#[derive(Debug, thiserror::Error)]
pub enum TofError {
    #[error("cannot read {0}: {1}")]
    Io(PathBuf, std::io::Error),
    #[error("invalid TOF calibration file: {0}")]
    Invalid(String),
    #[error("Not enough values for polynomial regression.")]
    NotEnoughPairs,
    #[error("the regression has no unique solution (the calibration ions are too similar)")]
    Singular,
}

/// `m/z = Σ coefficients[i] · bin^i` (`FilteredDataProvider::calculatePolynomial`).
#[derive(Clone, Debug, PartialEq)]
pub struct TofCalibration {
    pub coefficients: Vec<f64>,
}

impl TofCalibration {
    /// The m/z of a TOF bin; the old code kept m/z as `float`.
    pub fn mz(&self, bin: f32) -> f32 {
        let bin = f64::from(bin);
        self.coefficients
            .iter()
            .enumerate()
            .map(|(i, c)| c * bin.powi(i as i32))
            .sum::<f64>() as f32
    }

    pub fn read(path: impl AsRef<Path>) -> Result<Self, TofError> {
        let path = path.as_ref();
        let text =
            std::fs::read_to_string(path).map_err(|e| TofError::Io(path.to_path_buf(), e))?;
        Self::parse(&text)
    }

    /// `MainWindow::readTOFCalibData`: whitespace-separated tokens `Polynomial order: <n>`, then
    /// `Parameter: <i> <value>` for `i` in `0..=n`. Whatever follows (the calibration table) is ignored.
    pub fn parse(text: &str) -> Result<Self, TofError> {
        let invalid = |what: &str| TofError::Invalid(what.to_string());
        let mut tokens = text.split_whitespace();
        tokens.nth(1).ok_or_else(|| invalid("missing header"))?;
        let order: i32 = tokens
            .next()
            .and_then(|t| t.parse().ok())
            .ok_or_else(|| invalid("missing polynomial order"))?;
        // Real files stay far below this; it caps the allocation for a corrupt order.
        if !(0..=64).contains(&order) {
            return Err(TofError::Invalid(format!(
                "polynomial order {order} is outside 0..=64"
            )));
        }
        let mut coefficients = Vec::with_capacity(order as usize + 1);
        for i in 0..=order {
            // The old reader skips the label and the index.
            let value = tokens
                .nth(2)
                .and_then(|t| t.parse::<f64>().ok())
                .filter(|v| v.is_finite())
                .ok_or_else(|| TofError::Invalid(format!("parameter {i} is missing or invalid")))?;
            coefficients.push(value);
        }
        Ok(Self { coefficients })
    }

    /// Text of a `.mzcalib` file as `TOFCalibrationWindow::on_act_save_triggered` writes it:
    /// 18 significant digits, then the calibration pairs (TOF bin, m/z).
    pub fn to_text(&self, pairs: &[(f32, f64)]) -> String {
        let mut s = format!(
            "Polynomial order:\t{}\n",
            self.coefficients.len() as i64 - 1
        );
        for (i, c) in self.coefficients.iter().enumerate() {
            s += &format!("Parameter:\t{i}\t{}\n", g18(*c));
        }
        s += "Calibration table\n";
        for (bin, mz) in pairs {
            s += &format!("{}\t{}\n", g18(f64::from(*bin)), g18(*mz));
        }
        s
    }

    /// The calibration pairs (TOF bin, m/z) a `.mzcalib` file lists after its parameters; lines that
    /// are not two numbers are skipped.
    pub fn parse_pairs(text: &str) -> Vec<(f32, f64)> {
        text.lines()
            .skip_while(|l| !l.starts_with("Calibration table"))
            .skip(1)
            .filter_map(|l| {
                let mut f = l.split_whitespace();
                let (bin, mz) = (
                    f.next()?.parse::<f32>().ok()?,
                    f.next()?.parse::<f64>().ok()?,
                );
                (f.next().is_none() && bin.is_finite() && mz.is_finite()).then_some((bin, mz))
            })
            .collect()
    }

    pub fn write(&self, path: impl AsRef<Path>, pairs: &[(f32, f64)]) -> Result<(), TofError> {
        let path = path.as_ref();
        std::fs::write(path, self.to_text(pairs)).map_err(|e| TofError::Io(path.to_path_buf(), e))
    }

    /// Least-squares polynomial of `degree` through (TOF bin, m/z) pairs, like
    /// `TOFCalibrationCalculator::regress` (GSL `multifit_linear`, an SVD solution).
    pub fn fit(pairs: &[(f32, f64)], degree: usize) -> Result<Self, TofError> {
        if pairs.len() < degree + 1 {
            return Err(TofError::NotEnoughPairs);
        }
        let mut x = DMatrix::from_fn(pairs.len(), degree + 1, |r, c| {
            f64::from(pairs[r].0).powi(c as i32)
        });
        // Bins of 10^4 to the power of 3 next to 1: balance the columns first, like GSL does.
        let scale: Vec<f64> = x
            .column_iter()
            .map(|c| c.norm())
            .map(|n| if n > 0.0 { n } else { 1.0 })
            .collect();
        for (mut col, s) in x.column_iter_mut().zip(&scale) {
            col /= *s;
        }
        let y = DVector::from_iterator(pairs.len(), pairs.iter().map(|p| p.1));
        let svd = x.svd(true, true);
        if svd.rank(1e-12) < degree + 1 {
            return Err(TofError::Singular);
        }
        let solution = svd.solve(&y, 1e-12).map_err(|_| TofError::Singular)?;
        let coefficients: Vec<f64> = solution.iter().zip(&scale).map(|(c, s)| c / s).collect();
        match coefficients.iter().all(|c| c.is_finite()) {
            true => Ok(Self { coefficients }),
            false => Err(TofError::Singular),
        }
    }
}

/// The TOF bin of the peak near `near` in a summed bin spectrum (bin order): the intensity-weighted
/// centre of the points at or above `fwhm` (`SF_FWHM`) times the highest one within `half_width`
/// bins. None when the bins around `near` carry no signal.
pub fn peak_bin(spectrum: &[(f32, f64)], near: f32, half_width: f32, fwhm: f64) -> Option<f32> {
    let lo = spectrum.partition_point(|p| p.0 < near - half_width);
    let hi = spectrum.partition_point(|p| p.0 <= near + half_width);
    let window = &spectrum[lo..hi];
    let top = window.iter().map(|p| p.1).fold(0.0, f64::max);
    if top <= 0.0 {
        return None;
    }
    let (weighted, total) = window
        .iter()
        .filter(|p| p.1 >= top * fwhm)
        .fold((0.0, 0.0), |(w, t), p| (w + f64::from(p.0) * p.1, t + p.1));
    Some((weighted / total) as f32)
}

/// C's `%.18g`, which `QTextStream::setRealNumberPrecision(18)` selects for doubles.
fn g18(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return if v == 0.0 { "0".into() } else { v.to_string() };
    }
    const PRECISION: i32 = 18;
    let sci = format!("{:.*e}", PRECISION as usize - 1, v);
    let (mantissa, exp) = sci.split_once('e').expect("exponent form");
    let exp: i32 = exp.parse().expect("integer exponent");
    let trim = |s: &str| -> String {
        match s.contains('.') {
            true => s.trim_end_matches('0').trim_end_matches('.').to_string(),
            false => s.to_string(),
        }
    };
    if !(-4..PRECISION).contains(&exp) {
        format!(
            "{}e{}{:02}",
            trim(mantissa),
            if exp < 0 { '-' } else { '+' },
            exp.abs()
        )
    } else {
        trim(&format!("{:.*}", (PRECISION - 1 - exp) as usize, v))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates_the_polynomial() {
        let c = TofCalibration {
            coefficients: vec![1.0, 2.0, 3.0],
        };
        assert_eq!(c.mz(2.0), 17.0);
        assert_eq!(c.mz(0.0), 1.0);
    }

    #[test]
    fn mzcalib_text_round_trips() {
        let c = TofCalibration {
            coefficients: vec![-12.345678901234567, 0.0123456789, 2.5e-9],
        };
        let text = c.to_text(&[(1000.0, 28.0), (2500.5, 73.0)]);
        assert!(text.starts_with("Polynomial order:\t2\nParameter:\t0\t-12.34567890123456"));
        assert!(text.contains("Calibration table\n1000\t28\n2500.5\t73\n"));
        assert_eq!(TofCalibration::parse(&text).unwrap(), c);
    }

    #[test]
    fn numbers_are_printed_like_printf_g() {
        assert_eq!(g18(0.1), "0.100000000000000006");
        assert_eq!(g18(1e-5), "1.00000000000000008e-05");
        assert_eq!(g18(123456.0), "123456");
        assert_eq!(g18(1e18), "1e+18");
        assert_eq!(g18(-2.5), "-2.5");
        assert_eq!(g18(0.0), "0");
    }

    #[test]
    fn pairs_come_back_from_the_file_text() {
        let c = TofCalibration {
            coefficients: vec![1.0, 2.0],
        };
        let pairs = [(1000.0, 28.0), (2500.5, 73.25)];
        assert_eq!(TofCalibration::parse_pairs(&c.to_text(&pairs)), pairs);
        assert!(TofCalibration::parse_pairs("Polynomial order: 0").is_empty());
    }

    #[test]
    fn peak_bin_is_the_centre_of_the_peak_near_the_hint() {
        let spectrum: Vec<(f32, f64)> = (990..1011)
            .map(|b| (b as f32, 100.0 - 10.0 * ((b - 1000) as f64).abs()))
            .chain([(5000.0, 1e9)])
            .collect();
        assert_eq!(peak_bin(&spectrum, 1003.0, 8.0, 0.5), Some(1000.0));
        assert_eq!(peak_bin(&spectrum, 3000.0, 8.0, 0.5), None);
        // SF_FWHM decides which points of a skewed peak count.
        let skew = [(1000.0, 100.0), (1001.0, 60.0), (1002.0, 30.0)];
        assert_eq!(peak_bin(&skew, 1001.0, 8.0, 0.5), Some(1000.375));
        assert_eq!(peak_bin(&skew, 1001.0, 8.0, 0.7), Some(1000.0));
    }

    #[test]
    fn broken_files_are_errors() {
        for text in [
            "",
            "Polynomial order:",
            "Polynomial order:\tx\n",
            "Polynomial order:\t-1\n",
            "Polynomial order:\t1\nParameter:\t0\t1\n",
            "Polynomial order:\t0\nParameter:\t0\tabc\n",
            "Polynomial order:\t0\nParameter:\t0\tnan\n",
            "Polynomial order:\t99999999\n",
        ] {
            assert!(TofCalibration::parse(text).is_err(), "{text:?}");
        }
    }

    #[test]
    fn fit_recovers_the_polynomial() {
        let truth = TofCalibration {
            coefficients: vec![5.0, 0.01, 4e-8],
        };
        let pairs: Vec<_> = [800.0, 1500.0, 2600.0, 4000.0, 6100.0]
            .map(|b| (b, f64::from(truth.mz(b))))
            .to_vec();
        let fit = TofCalibration::fit(&pairs, 2).unwrap();
        for (bin, mz) in &pairs {
            assert!((f64::from(fit.mz(*bin)) - mz).abs() < 1e-3, "bin {bin}");
        }
        let same = [(1000.0, 5.0), (1000.0, 6.0), (1000.0, 7.0)];
        assert!(matches!(
            TofCalibration::fit(&same, 2),
            Err(TofError::Singular)
        ));
        assert!(matches!(
            TofCalibration::fit(&pairs[..2], 2),
            Err(TofError::NotEnoughPairs)
        ));
    }
}
