//! Mass isotopomer distributions (MID) of labelled fragments, the old `tools::MIDSolver`
//! (algorithms.md §7g). A measured fragment spectrum `M+0 … M+n` is explained as a mix of the
//! fragment with 0, 1, 2 … labelled atoms; the shares of that mix are the MID. Every correction is a
//! linear least-squares fit; R² and confidence intervals come from the same fit as in the old program.

mod export;
mod frag;
mod report;
mod series;

pub use report::{MidReport, Move};

pub use export::{group_stat, qt_g};

pub use frag::{FragExport, frag_export, parse_fragments};
pub use series::{
    CompoundMid, Fragment, FragmentMid, Method, Mids, determine, fragments, intensities,
    mid_library,
};

use crate::isotope::{Formula, IsotopeTable};

/// Widest fragment taken from library text or a `FRAG` feature: GC-MS spectra end near m/z 1000,
/// so anything wider is a typo that would freeze the scan or the export.
pub const MAX_FRAGMENT_IONS: u32 = 1000;

/// One fitted MID. `None` stands for the old program's `-1` ("not available").
#[derive(Debug, Clone, PartialEq)]
pub struct MidFit {
    /// Share of `M+0`, `M+1`, … labelled atoms.
    pub mids: Vec<f64>,
    /// `1 − RSS/TSS` of the fit; `None` for replicate statistics.
    pub r2: Option<f64>,
    /// Half width of the confidence interval of each share; `None` when the fit has no degree of
    /// freedom left (or for a single replicate).
    pub cis: Vec<Option<f64>>,
}

/// `v / Σv`.
fn normalize(v: &[f64]) -> Vec<f64> {
    let sum: f64 = v.iter().sum();
    v.iter().map(|x| x / sum).collect()
}

/// Correction by the sum formula of the fragment (old `MIDSolver(vv, formula, calc, tracer, n)`).
/// Column `k` of the design matrix is the natural pattern of the fragment with `k` tracer atoms
/// replaced by the label (one atom fewer of `tracer` per column, never below 0), shifted by `k`.
/// `size` is the number of shares (0 = one per measured ion). `None` for an empty spectrum.
pub fn correct_by_formula(
    labelled: &[f64],
    formula: &Formula,
    table: &IsotopeTable,
    tracer: Option<&str>,
    size: usize,
) -> Result<Option<MidFit>, String> {
    if labelled.is_empty() {
        return Ok(None);
    }
    let len = labelled.len();
    let cols = if size == 0 { len } else { size.min(len) };
    let mut atoms = formula.clone();
    let mut columns = Vec::with_capacity(cols);
    for k in 0..cols {
        if let Some(n) = tracer.filter(|_| k > 0).and_then(|t| atoms.get_mut(t)) {
            *n = n.saturating_sub(1);
        }
        columns.push(normalize(&table.pattern(&atoms, len)?.abundance));
    }
    let m: Vec<Vec<f64>> = (0..len)
        .map(|j| {
            (0..cols)
                .map(|k| if j >= k { columns[k][j - k] } else { 0.0 })
                .collect()
        })
        .collect();
    Ok(Some(fit(&m, &normalize(labelled))))
}

/// Correction by unlabelled reference spectra of the same fragment (old `MIDSolver(vv, vv, double)`).
/// Every pair (reference, labelled spectrum) adds one equation per ion: the labelled spectrum is
/// the reference spectrum shifted by `k` for `k` labelled atoms. `m1` corrects the shifted `M+0`/`M+1`
/// of the reference for the tracer's own natural heavy isotope (0 = off). `None` when either side is empty.
pub fn correct_by_reference(
    labelled: &[Vec<f64>],
    unlabelled: &[Vec<f64>],
    m1: f64,
) -> Option<MidFit> {
    let cols = labelled.first()?.len();
    unlabelled.first()?;
    let (mut m, mut y) = (Vec::new(), Vec::new());
    for u in unlabelled {
        // shortcut: a reference longer than the labelled spectra is cut to their length (the old
        // program's behaviour there is untested, file-formats.md does not apply).
        let r = normalize(&u[..cols.min(u.len())]);
        let at = |i: usize| r.get(i).copied().unwrap_or(0.0);
        for l in labelled {
            let target = normalize(l);
            for n in 0..cols {
                let row = (0..cols)
                    .map(|j| {
                        let cf = j as f64 * m1;
                        let shift = at(0) * cf / ((at(1) / at(0) + 1.0) - cf);
                        match n.checked_sub(j) {
                            Some(0) if n != 0 => (shift + at(0)).max(0.0),
                            Some(1) => (at(1) - shift).max(0.0),
                            Some(d) => at(d),
                            None => 0.0,
                        }
                    })
                    .collect();
                m.push(row);
                y.push(target.get(n).copied().unwrap_or(0.0));
            }
        }
    }
    Some(fit(&m, &y))
}

/// Mean MID of replicates and `1.96 · sd / √n` (old `MIDSolver(vv)`). `None` without replicates.
pub fn replicate_mean(rows: &[Vec<f64>]) -> Option<MidFit> {
    let cols = rows.first()?.len();
    let rows: Vec<Vec<f64>> = rows.iter().map(|r| normalize(r)).collect();
    let n = rows.len() as f64;
    let (mut mids, mut cis) = (Vec::with_capacity(cols), Vec::with_capacity(cols));
    for j in 0..cols {
        let v: Vec<f64> = rows
            .iter()
            .map(|r| r.get(j).copied().unwrap_or(0.0))
            .collect();
        let mean = v.iter().sum::<f64>() / n;
        mids.push(if rows.len() > 1 { mean } else { v[0] });
        cis.push((rows.len() > 1).then(|| {
            let sd = (v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
            1.96 * (sd / n.sqrt())
        }));
    }
    Some(MidFit {
        mids,
        r2: None,
        cis,
    })
}

/// Ordinary least squares `m · x ≈ y` like `gsl_multifit_linear`, with R² and the 80 % two-sided
/// (t at 0.9) confidence half widths of the old program.
fn fit(m: &[Vec<f64>], y: &[f64]) -> MidFit {
    let (rows, cols) = (y.len(), m[0].len());
    let (x, cov_diag, rss) = least_squares(m, y);
    let mean = y.iter().sum::<f64>() / rows as f64;
    let tss: f64 = y.iter().map(|v| (v - mean).powi(2)).sum();
    let cis = if rows > cols {
        let df = rows - cols;
        let t = t_quantile(0.9, df);
        let sigma2 = rss / df as f64;
        cov_diag
            .iter()
            .map(|c| Some((c * sigma2).sqrt() * t))
            .collect()
    } else {
        vec![None; cols]
    };
    MidFit {
        mids: x,
        r2: Some(1.0 - rss / tss),
        cis,
    }
}

/// Householder QR solution of `m · x ≈ y`: `x`, the diagonal of `(mᵀm)⁻¹` and the residual sum of
/// squares. shortcut: no rank cutoff like GSL's SVD; a rank-deficient `m` gives non-finite values.
// Matrix code: row and column indices read clearer than iterator chains.
#[allow(clippy::needless_range_loop)]
fn least_squares(m: &[Vec<f64>], y: &[f64]) -> (Vec<f64>, Vec<f64>, f64) {
    let (rows, cols) = (y.len(), m[0].len());
    let mut a: Vec<Vec<f64>> = m.to_vec();
    let mut b = y.to_vec();
    for k in 0..cols.min(rows) {
        let norm = (k..rows).map(|i| a[i][k].powi(2)).sum::<f64>().sqrt();
        if norm == 0.0 {
            continue;
        }
        let alpha = if a[k][k] >= 0.0 { -norm } else { norm };
        let mut v = vec![0.0; rows];
        v[k] = a[k][k] - alpha;
        for i in k + 1..rows {
            v[i] = a[i][k];
        }
        let vv: f64 = v[k..].iter().map(|x| x * x).sum();
        if vv == 0.0 {
            continue;
        }
        for j in k..cols {
            let s = 2.0 * (k..rows).map(|i| v[i] * a[i][j]).sum::<f64>() / vv;
            (k..rows).for_each(|i| a[i][j] -= s * v[i]);
        }
        let s = 2.0 * (k..rows).map(|i| v[i] * b[i]).sum::<f64>() / vv;
        (k..rows).for_each(|i| b[i] -= s * v[i]);
    }
    let r = |i: usize, j: usize| if i < rows { a[i][j] } else { 0.0 };
    let mut x = vec![0.0; cols];
    for i in (0..cols).rev() {
        let s: f64 = (i + 1..cols).map(|j| r(i, j) * x[j]).sum();
        x[i] = (b.get(i).copied().unwrap_or(0.0) - s) / r(i, i);
    }
    // R⁻¹ column by column; diag((mᵀm)⁻¹) = row sums of squares of R⁻¹.
    let mut inv = vec![vec![0.0; cols]; cols];
    for c in 0..cols {
        for i in (0..cols).rev() {
            let s: f64 = (i + 1..cols).map(|j| r(i, j) * inv[j][c]).sum();
            inv[i][c] = (f64::from(u8::from(i == c)) - s) / r(i, i);
        }
    }
    let cov_diag = inv
        .iter()
        .map(|row| row.iter().map(|v| v * v).sum())
        .collect();
    let rss = b.iter().skip(cols).map(|v| v * v).sum();
    (x, cov_diag, rss)
}

/// Student t distribution function for `df` degrees of freedom (exact series for integer `df`,
/// Abramowitz & Stegun 26.7.3/26.7.4).
fn t_cdf(t: f64, df: usize) -> f64 {
    let theta = (t / (df as f64).sqrt()).atan();
    let (s, c) = theta.sin_cos();
    let c2 = c * c;
    // A(t|df) = P(|T| < t) · sign(t).
    let a = if df % 2 == 1 {
        let (mut term, mut sum) = (c, if df > 1 { c } else { 0.0 });
        for k in 1..(df.saturating_sub(1)) / 2 {
            term *= c2 * (2 * k) as f64 / (2 * k + 1) as f64;
            sum += term;
        }
        2.0 / std::f64::consts::PI * (theta + s * sum)
    } else {
        let (mut term, mut sum) = (1.0, 1.0);
        for k in 1..df / 2 {
            term *= c2 * (2 * k - 1) as f64 / (2 * k) as f64;
            sum += term;
        }
        s * sum
    };
    0.5 + 0.5 * a
}

/// `t` with `t_cdf(t, df) = p` for `p` in (0.5, 1), by bisection to the last bit.
fn t_quantile(p: f64, df: usize) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 1e4f64);
    while hi - lo > f64::EPSILON * hi {
        let mid = 0.5 * (lo + hi);
        if t_cdf(mid, df) < p {
            lo = mid
        } else {
            hi = mid
        }
    }
    0.5 * (lo + hi)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isotope::parse_formula;

    fn close(actual: &[f64], expected: &[f64], tol: f64) {
        assert_eq!(actual.len(), expected.len(), "{actual:?} != {expected:?}");
        for (a, e) in actual.iter().zip(expected) {
            assert!((a - e).abs() <= tol, "{actual:?} != {expected:?}");
        }
    }

    fn check(fit: &MidFit, mids: &[f64], r2: Option<f64>, cis: &[f64]) {
        close(&fit.mids, mids, 1e-12);
        match (fit.r2, r2) {
            (Some(a), Some(e)) => close(&[a], &[e], 1e-12),
            (a, e) => assert_eq!(a, e),
        }
        let got: Vec<f64> = fit.cis.iter().map(|c| c.unwrap_or(-1.0)).collect();
        close(&got, cis, 1e-12);
    }

    /// Outputs of `tools::MIDSolver` of MetaboliteDetector 3.3 (called through `LD_PRELOAD`) for
    /// synthetic spectra, with the old program's isotope values.
    #[test]
    #[allow(clippy::excessive_precision)] // the old program's output, verbatim
    fn matches_old_program() {
        let table = IsotopeTable::new([
            ("H", 1, 1.007825, 99.9844),
            ("H", 2, 2.0141, 0.0156),
            ("C", 12, 12.0, 98.9184),
            ("C", 13, 13.0034, 1.0816),
            ("O", 16, 15.9949, 99.758),
            ("O", 17, 16.9991, 0.038),
            ("O", 18, 17.9992, 0.204),
        ]);
        let glucose = parse_formula("C6H12O6").unwrap();
        let fit = |l: &[f64], size| correct_by_formula(l, &glucose, &table, Some("C"), size);
        check(
            &fit(
                &[
                    0.04382309380000551,
                    0.0987630827695157,
                    0.2915233996369361,
                    0.254251869638238,
                    0.15731233107088466,
                    0.10482616059912415,
                    0.05412568349566571,
                ],
                4,
            )
            .unwrap()
            .unwrap(),
            &[
                0.047320823129663403,
                0.10217506321779779,
                0.30322823565178686,
                0.25698329844437467,
            ],
            Some(0.38858621085933787),
            &[
                0.18861810372063209,
                0.18703147316425203,
                0.18487671143216131,
                0.18264823417345924,
            ],
        );
        check(
            &fit(
                &[
                    0.048951803624578794,
                    0.09800667116419133,
                    0.2870567819939816,
                    0.2517667490911604,
                    0.15982399596783295,
                    0.10528692504043535,
                    0.055102060828925506,
                ],
                5,
            )
            .unwrap()
            .unwrap(),
            &[
                0.052800647197742689,
                0.10091623504536015,
                0.29631191806607987,
                0.24861725772962037,
                0.15538764907650482,
            ],
            Some(0.76427398416198333),
            &[
                0.16056684180352371,
                0.15921618950255548,
                0.15739209644863147,
                0.15558943869605929,
                0.15378470613718168,
            ],
        );
        let labelled = [
            vec![
                0.7762520132840594,
                0.21166490115481168,
                0.9786028078247297,
                0.49111174221368115,
                0.310014768707883,
            ],
            vec![
                0.29113770434592806,
                0.28679215309817624,
                0.33361617929582166,
                0.06253042741271105,
                0.587724297964751,
            ],
        ];
        let unlabelled = [
            vec![
                0.8249686262185113,
                0.8210097435995273,
                0.9109037469481489,
                0.4099651277163867,
                0.0757021611274613,
            ],
            vec![
                0.6710963730114814,
                0.3583320146969797,
                0.1290084770503076,
                0.36017615295558836,
                0.5234683169759594,
            ],
            vec![
                0.8115647885441393,
                0.21070352094904599,
                0.6688142667861245,
                0.5577813234574209,
                0.20266365561123653,
            ],
        ];
        check(
            &correct_by_reference(&labelled, &unlabelled, 0.05).unwrap(),
            &[
                0.62472158165702241,
                0.031057473827248716,
                0.44233221415869961,
                -0.1900164190916683,
                0.25651136545551945,
            ],
            Some(0.24091137227495485),
            &[
                0.14390967470643817,
                0.16761870474437382,
                0.16418782535841606,
                0.17123602302101743,
                0.17268999294115819,
            ],
        );
        check(
            &correct_by_reference(
                &[vec![
                    0.22436419552056724,
                    0.8765244170234823,
                    0.06812592179272249,
                ]],
                &[vec![
                    0.6986666457683174,
                    0.7234427664281333,
                    0.783942718174666,
                ]],
                0.02,
            )
            .unwrap(),
            &[0.6060094576997368, 1.7229027279632516, -2.2184020858245757],
            Some(1.0),
            &[-1.0, -1.0, -1.0],
        );
        let replicates = [
            vec![
                0.8339726048754931,
                0.11554833472238987,
                0.9619514544943147,
                0.5648473167860983,
            ],
            vec![
                0.4179681265452158,
                0.23541613809925818,
                0.7452328214431115,
                0.46509541421524303,
            ],
            vec![
                0.6386930910354672,
                0.059712435921681725,
                0.3972228347614469,
                0.7575581826204791,
            ],
            vec![
                0.4996298234121156,
                0.8981173639117144,
                0.05628849528117099,
                0.3739747428361489,
            ],
            vec![
                0.9565618520213375,
                0.8204779707774443,
                0.6755386525485985,
                0.7858776164879152,
            ],
            vec![
                0.3750673895436141,
                0.3441684416135428,
                0.3780193578561891,
                0.9783779130919619,
            ],
        ];
        check(
            &replicate_mean(&replicates).unwrap(),
            &[
                0.27584773006241814,
                0.18594592576435298,
                0.23736404909779218,
                0.30084229507543675,
            ],
            None,
            &[
                0.05129465386815206,
                0.13616138309887529,
                0.11105742652248796,
                0.088599783893362372,
            ],
        );
    }

    #[test]
    fn t_quantiles_match_tables() {
        // t_{0.9}: 3.0777 (1), 1.8856 (2), 1.4759 (5), 1.3830 (9), 1.2816 (∞).
        for (df, t) in [
            (1, 3.077_683_537_175),
            (2, 1.885_618_083_164),
            (5, 1.475_884_048_824),
            (9, 1.383_028_738_396),
        ] {
            assert!((t_quantile(0.9, df) - t).abs() < 1e-11, "{df}");
        }
        assert!((t_quantile(0.9, 10_000) - 1.281_6).abs() < 1e-3);
    }

    #[test]
    fn exact_label_is_recovered() {
        let table = IsotopeTable::standard();
        let formula = parse_formula("C6H12O6").unwrap();
        let truth = [0.1, 0.2, 0.3, 0.4];
        // Measured = Σ share_k · (pattern with k fewer C) shifted by k.
        let mut atoms = formula.clone();
        let mut measured = vec![0.0; 7];
        for (k, share) in truth.iter().enumerate() {
            if k > 0 {
                *atoms.get_mut("C").unwrap() -= 1;
            }
            let p = normalize(&table.pattern(&atoms, 7).unwrap().abundance);
            for j in k..7 {
                measured[j] += share * p[j - k];
            }
        }
        let fit = correct_by_formula(&measured, &formula, &table, Some("C"), 4)
            .unwrap()
            .unwrap();
        // The fit sees the spectrum normalised to the 7 measured ions, so the shares scale with it.
        let scale: f64 = measured.iter().sum();
        let mids: Vec<f64> = fit.mids.iter().map(|m| m * scale).collect();
        close(&mids, &truth, 1e-12);
        assert!((fit.r2.unwrap() - 1.0).abs() < 1e-12);
        assert!(
            correct_by_formula(&[], &formula, &table, Some("C"), 0)
                .unwrap()
                .is_none()
        );
        assert!(correct_by_reference(&[], &[vec![1.0]], 0.0).is_none());
        assert!(replicate_mean(&[]).is_none());
    }
}
