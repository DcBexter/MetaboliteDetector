//! The old MID result windows (algorithms.md §7g.7): group statistics of the per-run MIDs
//! (`ReplicateMIDData`) and the text exports, byte for byte with Qt's number format.

use std::fmt::Write;

use super::{FragmentMid, Method, MidFit, Mids};
use crate::library::LibraryCompound;

/// A number as Qt's `QTextStream` writes it: `%g` with 6 significant digits, `nan`, `inf`, and `0` for
/// -0 (printf: `-0`).
pub fn qt_g(v: f64) -> String {
    if v.is_nan() {
        return "nan".into();
    }
    if v.is_infinite() {
        return if v < 0.0 { "-inf" } else { "inf" }.into();
    }
    if v == 0.0 {
        return "0".into();
    }
    // The exponent after rounding to 6 digits decides between the two notations, as in printf.
    let sci = format!("{v:.5e}");
    let (mantissa, exp) = sci.split_once('e').expect("exponent");
    let exp: i32 = exp.parse().expect("exponent");
    let trim = |s: String| match s.contains('.') {
        true => s.trim_end_matches('0').trim_end_matches('.').to_owned(),
        false => s,
    };
    if (-4..6).contains(&exp) {
        trim(format!("{v:.*}", (5 - exp) as usize))
    } else {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{sign}{:02}", trim(mantissa.to_owned()), exp.abs())
    }
}

/// Mean and standard error of one isotopomer over the labelled runs (old `ReplicateMIDData`, group
/// `GL`, always called with the filter on): shares that are 0, NaN or infinite do not count, as runs
/// without a fit do not. No value gives mean 0, fewer than two give a NaN standard error (GSL).
pub fn group_stat(f: &FragmentMid, k: usize) -> (f64, f64) {
    let v: Vec<f64> = (f.fits.iter().flatten())
        .filter_map(|x| x.mids.get(k).copied())
        .filter(|v| *v != 0.0 && v.is_finite())
        .collect();
    let n = v.len() as f64;
    let mean = match v.is_empty() {
        true => 0.0,
        false => v.iter().sum::<f64>() / n,
    };
    let stderr = match v.len() {
        0 | 1 => f64::NAN,
        _ => (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt() / n.sqrt(),
    };
    (mean, stderr)
}

impl Mids {
    /// Export of the group methods (old `MIDResultWidget::actExportTriggered` 0x4e9212): per
    /// fragment a row `Name_ion (Mk)`, MID, CI and R² (-1 where there is none), then an empty line.
    /// A fragment without MID is the old solver's empty result: one row `M0` with -1 everywhere.
    pub fn group_text(&self, lib: &[LibraryCompound]) -> String {
        let none = MidFit {
            mids: vec![-1.0],
            r2: None,
            cis: vec![None],
        };
        let mut s = String::new();
        for c in &self.compounds {
            let name = &lib[c.compound].name;
            for f in &c.fragments {
                for fit in f.fits.iter().map(|x| x.as_ref().unwrap_or(&none)) {
                    for (k, m) in fit.mids.iter().take(f.isotopomers).enumerate() {
                        let ci = fit.cis.get(k).copied().flatten().unwrap_or(-1.0);
                        let r2 = fit.r2.unwrap_or(-1.0);
                        _ = writeln!(
                            s,
                            "{name}_{} (M{k})\t{}\t{}\t{}",
                            f.lo as usize + k,
                            qt_g(*m),
                            qt_g(ci),
                            qt_g(r2)
                        );
                    }
                }
                s.push('\n');
            }
        }
        s
    }

    /// Export of the per-run methods (old `MIDSingleResultWidget::outputMIDDataToStream` 0x4ef548):
    /// a column per run for the MIDs, the CIs and R². A run without a fit is NaN, its CI too with the
    /// sum formula and -1 with the library spectrum, as in the old program.
    pub fn run_text(&self, lib: &[LibraryCompound], method: Method, columns: &[String]) -> String {
        let mut s = String::new();
        for (pre, post) in [("", ""), ("CI (", ")"), ("R2 (", ")")] {
            columns
                .iter()
                .for_each(|c| _ = write!(s, "\t{pre}{c}{post}"));
            s.push('\t');
        }
        s.pop();
        s.push('\n');
        let missing_ci = match method {
            Method::Library => -1.0,
            _ => f64::NAN,
        };
        let cell = |x: Option<&super::MidFit>, column: usize, k: usize| match (x, column) {
            (None, 1) => missing_ci,
            (None, _) => f64::NAN,
            (Some(x), 0) => x.mids.get(k).copied().unwrap_or(f64::NAN),
            (Some(x), 1) => x.cis.get(k).copied().flatten().unwrap_or(-1.0),
            (Some(x), _) => x.r2.unwrap_or(-1.0),
        };
        self.blocks(lib, &mut s, |s, _, f, k| {
            for column in 0..3 {
                for fit in &f.fits {
                    _ = write!(s, "\t{}", qt_g(cell(fit.as_ref(), column, k)));
                }
                s.push('\t');
            }
            s.pop();
        });
        s
    }

    /// Statistics export of the per-run methods (old `outputMIDStatisticsDataToStream` 0x4f0090):
    /// mean and standard error of the labelled runs ([`group_stat`]) per isotopomer, then substance,
    /// first fragment ion and isotopomer again.
    pub fn stats_text(&self, lib: &[LibraryCompound]) -> String {
        let mut s =
            String::from("\tGL\t\tStdErr (GL)\tMetabolite name\tFragment ion\tMass isotopomer\n");
        self.blocks(lib, &mut s, |s, name, f, k| {
            let (mean, stderr) = group_stat(f, k);
            _ = write!(
                s,
                "\t{}\t\t{}\t{name}\t{}\tM{k}",
                qt_g(mean),
                qt_g(stderr),
                f.lo
            );
        });
        s
    }

    /// The substance blocks both per-run exports share: `Name;RI;RT`, then a row `ion (Mk)` per
    /// isotopomer with an empty line before each fragment, then two line breaks.
    fn blocks(
        &self,
        lib: &[LibraryCompound],
        s: &mut String,
        row: impl Fn(&mut String, &str, &FragmentMid, usize),
    ) {
        for c in &self.compounds {
            let lc = &lib[c.compound];
            let sp = &lc.spectrum;
            _ = write!(s, "{};{};{}", lc.name, qt_g(sp.ri), qt_g(sp.rt_ms));
            for f in &c.fragments {
                for k in 0..f.isotopomers {
                    if k == 0 {
                        s.push('\n');
                    }
                    _ = write!(s, "{} (M{k})", f.lo as usize + k);
                    row(s, &lc.name, f, k);
                    s.push('\n');
                }
            }
            s.push_str("\n\n");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mid::MidFit;

    #[test]
    fn numbers_print_like_qt() {
        let cases = [
            (1042570.0, "1.04257e+06"),
            (3599.99, "3599.99"),
            (-15.0602, "-15.0602"),
            (0.0236703, "0.0236703"),
            (0.0001, "0.0001"),
            (0.00001234, "1.234e-05"),
            (999999.5, "1e+06"),
            (123456.0, "123456"),
            (-1.0, "-1"),
            (0.75, "0.75"),
            (0.0, "0"),
            (-0.0, "0"),
            (f64::NAN, "nan"),
            (f64::NEG_INFINITY, "-inf"),
        ];
        for (v, want) in cases {
            assert_eq!(qt_g(v), want, "{v}");
        }
    }

    #[test]
    fn group_statistics_skip_zero_and_missing_shares() {
        let fit = |m: f64| {
            Some(MidFit {
                mids: vec![m],
                r2: None,
                cis: vec![None],
            })
        };
        let f = |fits| FragmentMid {
            lo: 1,
            hi: 1,
            isotopomers: 1,
            fits,
        };
        let (mean, se) = group_stat(
            &f(vec![fit(1.0), None, fit(0.0), fit(3.0), fit(f64::NAN)]),
            0,
        );
        assert_eq!(mean, 2.0);
        // sd of (1, 3) is sqrt(2), divided by sqrt(2).
        assert!((se - 1.0).abs() < 1e-12);
        let (mean, se) = group_stat(&f(vec![fit(0.5), None]), 0);
        assert!(mean == 0.5 && se.is_nan());
        let (mean, se) = group_stat(&f(vec![None]), 0);
        assert!(mean == 0.0 && se.is_nan());
    }
}
