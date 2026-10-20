//! Biomarker statistics, box plot numbers and PCA of the batch table (F14, F15,
//! `docs/reverse-engineering/algorithms.md` §7f). Pure functions on [`Table`]; the table is the
//! one the user sees, so normalisation (`normalize`) is applied before.

use std::io::Write;

use nalgebra::DMatrix;

use crate::batch::Table;

/// How the runs of the two groups become the samples of the test.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Samples {
    /// Every run is a sample; Welch's t-test (the old "independent samples").
    #[default]
    Independent,
    /// One sample per replicate group, its mean over the runs; Welch's t-test (the old "not
    /// independent samples").
    GroupMeans,
    /// Group k of side 1 pairs with group k of side 2 (in the order given), run by run in table
    /// order; paired t-test. Needs the same number of groups and runs per pair. Not in the old
    /// program.
    Paired,
}

/// Two replicate-group lists to compare.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoGroups {
    pub group1: Vec<String>,
    pub group2: Vec<String>,
    pub samples: Samples,
    /// A substance not found in a run counts as 0 (the old table); else the run is left out.
    pub missing_as_zero: bool,
}

/// One substance, group 1 against group 2.
#[derive(Debug, Clone, PartialEq)]
pub struct Biomarker {
    /// Index into `Table::rows`.
    pub row: usize,
    /// Runs with a value, per group.
    pub hits: [usize; 2],
    pub mean: [f64; 2],
    /// Sample standard deviation (N-1); NaN below two samples.
    pub sd: [f64; 2],
    /// `mean[0] / mean[1]`; not guarded, like the old program (inf/NaN).
    pub ratio: f64,
    /// Two-sided p-value; `None` where the test is undefined (paired with unequal or fewer than
    /// two pairs, or no spread in the differences).
    pub p: Option<f64>,
}

/// Values of one substance in the runs of `groups`, per group in table order. `None` is a missing
/// value.
fn cells<'a>(
    t: &'a Table,
    row: usize,
    groups: &'a [String],
) -> impl Iterator<Item = (&'a str, Option<f64>)> + 'a {
    t.groups
        .iter()
        .zip(&t.rows[row].values)
        .filter(|(g, _)| groups.contains(g))
        .map(|(g, v)| (g.as_str(), *v))
}

/// The values of `groups` for one substance (all runs of those groups, table order).
pub fn group_values(t: &Table, row: usize, groups: &[String], missing_as_zero: bool) -> Vec<f64> {
    cells(t, row, groups)
        .filter_map(|(_, v)| v.or(missing_as_zero.then_some(0.0)))
        .collect()
}

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

/// Sample variance, N-1.
fn variance(x: &[f64]) -> f64 {
    // Equal values have exactly no spread; `mean` of them can be off by a rounding error.
    if x.windows(2).all(|w| w[0] == w[1]) {
        return 0.0;
    }
    let m = mean(x);
    x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (x.len() as f64 - 1.0)
}

/// The batch table's per-group columns of one substance (old `BatchQuantificationStatisticsTableModel`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Summary {
    /// Samples with a value above 0 ("#").
    pub hits: usize,
    /// Mean of the values that are not 0, NaN or infinite; 0 if there are none.
    pub mean: f64,
    /// Standard error over the same values divided by their mean, as a fraction
    /// ("Rel. Std. Err."); NaN below two values.
    pub rel_std_err: f64,
}

/// `x` are the group's values of one substance, a missing value as 0.
pub fn summary(x: &[f64]) -> Summary {
    let kept: Vec<f64> = x
        .iter()
        .copied()
        .filter(|v| *v != 0.0 && v.is_finite())
        .collect();
    let m = if kept.is_empty() { 0.0 } else { mean(&kept) };
    let sd = if kept.len() > 1 {
        variance(&kept).sqrt()
    } else {
        f64::NAN
    };
    Summary {
        hits: x.iter().filter(|v| **v > 0.0).count(),
        mean: m,
        rel_std_err: sd / (kept.len() as f64).sqrt() / m,
    }
}

/// Two-sided p of Student's t-test with pooled variance (the old "P (blank)"); `None` below three
/// values in total, where the degrees of freedom vanish.
pub fn student(x: &[f64], y: &[f64]) -> Option<f64> {
    let (n1, n2) = (x.len() as f64, y.len() as f64);
    let df = n1 + n2 - 2.0;
    if df < 1.0 || x.is_empty() || y.is_empty() {
        return None;
    }
    let pooled = ((n1 - 1.0) * variance_or_zero(x) + (n2 - 1.0) * variance_or_zero(y)) / df;
    let t = ((mean(x) - mean(y)) / pooled.sqrt() * (n1 * n2 / (n1 + n2)).sqrt()).abs();
    t.is_finite().then(|| t_two_sided(t, df))
}

/// A single value has no spread; the pooled test counts it as 0.
fn variance_or_zero(x: &[f64]) -> f64 {
    if x.len() > 1 { variance(x) } else { 0.0 }
}

/// The samples of the test for one substance and group list.
fn samples(t: &Table, row: usize, groups: &[String], c: &TwoGroups) -> Vec<f64> {
    match c.samples {
        Samples::GroupMeans => groups
            .iter()
            .filter_map(|g| {
                let x = group_values(t, row, std::slice::from_ref(g), c.missing_as_zero);
                (!x.is_empty()).then(|| mean(&x))
            })
            .collect(),
        _ => group_values(t, row, groups, c.missing_as_zero),
    }
}

/// Group 1 against group 2 for every substance of the table.
pub fn compare(t: &Table, c: &TwoGroups) -> Vec<Biomarker> {
    (0..t.rows.len())
        .map(|row| {
            let hits = [&c.group1, &c.group2]
                .map(|g| cells(t, row, g).filter(|(_, v)| v.is_some()).count());
            let (x, y) = (samples(t, row, &c.group1, c), samples(t, row, &c.group2, c));
            let sd = [&x, &y].map(|s| {
                if s.len() > 1 {
                    variance(s).sqrt()
                } else {
                    f64::NAN
                }
            });
            let (m1, m2) = (mean(&x), mean(&y));
            let p = match c.samples {
                Samples::Paired => paired(t, row, c),
                _ => Some(welch(&x, &y)),
            };
            Biomarker {
                row,
                hits,
                mean: [m1, m2],
                sd,
                ratio: m1 / m2,
                p,
            }
        })
        .collect()
}

/// Whether a paired test is possible: group i of side 1 is paired with group i of side 2, so both
/// sides need the same number of groups and each pair the same number of runs.
pub fn pairs_match(t: &Table, c: &TwoGroups) -> bool {
    c.group1.len() == c.group2.len()
        && c.group1.iter().zip(&c.group2).all(|(a, b)| {
            let runs = |g: &String| t.groups.iter().filter(|x| *x == g).count();
            runs(a) == runs(b)
        })
}

/// Paired t-test: the i-th run of group k on side 1 against the i-th run of group k on side 2
/// (the order the user gave the groups). Pairs where a value is missing are left out.
fn paired(t: &Table, row: usize, c: &TwoGroups) -> Option<f64> {
    if !pairs_match(t, c) {
        return None;
    }
    let pick = |g: &String| -> Vec<Option<f64>> {
        cells(t, row, std::slice::from_ref(g))
            .map(|(_, v)| v.or(c.missing_as_zero.then_some(0.0)))
            .collect()
    };
    let d: Vec<f64> = c
        .group1
        .iter()
        .zip(&c.group2)
        .flat_map(|(a, b)| pick(a).into_iter().zip(pick(b)))
        .filter_map(|(a, b)| Some(a? - b?))
        .collect();
    let n = d.len() as f64;
    let se = (variance(&d) / n).sqrt();
    (d.len() > 1 && se > 0.0).then(|| t_two_sided(mean(&d) / se, n - 1.0))
}

/// Welch's t-test, two-sided p (alglib `unequalvariancettest`, `both`): 1 for fewer than two values
/// or no spread in a sample.
pub fn welch(x: &[f64], y: &[f64]) -> f64 {
    if x.len() < 2 || y.len() < 2 {
        return 1.0;
    }
    let (nx, ny) = (x.len() as f64, y.len() as f64);
    let (a, b) = (variance(x) / nx, variance(y) / ny);
    if a == 0.0 || b == 0.0 {
        return 1.0;
    }
    let df = (a + b).powi(2) / (a * a / (nx - 1.0) + b * b / (ny - 1.0));
    t_two_sided((mean(x) - mean(y)) / (a + b).sqrt(), df)
}

/// P(|T| >= |t|) for Student's t with `df` degrees of freedom: `I_{df/(df+t²)}(df/2, 1/2)`.
fn t_two_sided(t: f64, df: f64) -> f64 {
    regularized_beta(df / 2.0, 0.5, df / (df + t * t))
}

/// Regularised incomplete beta function `I_x(a, b)` (continued fraction of Lentz's method,
/// Numerical Recipes §6.4). Accurate to about 1e-14 for the arguments used here.
fn regularized_beta(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    if x >= 1.0 {
        return 1.0;
    }
    let front =
        (ln_gamma(a + b) - ln_gamma(a) - ln_gamma(b) + a * x.ln() + b * (1.0 - x).ln()).exp();
    // The fraction converges quickly below the mean of the distribution; use the symmetry above.
    if x < (a + 1.0) / (a + b + 2.0) {
        front * beta_fraction(a, b, x) / a
    } else {
        1.0 - front * beta_fraction(b, a, 1.0 - x) / b
    }
}

fn beta_fraction(a: f64, b: f64, x: f64) -> f64 {
    const TINY: f64 = 1e-300;
    let (mut c, mut d) = (1.0, 1.0 - (a + b) * x / (a + 1.0));
    if d.abs() < TINY {
        d = TINY;
    }
    d = 1.0 / d;
    let mut h = d;
    for m in 1..1000 {
        let m = f64::from(m);
        // Even step, then odd step of the continued fraction.
        for aa in [
            m * (b - m) * x / ((a + 2.0 * m - 1.0) * (a + 2.0 * m)),
            -(a + m) * (a + b + m) * x / ((a + 2.0 * m) * (a + 2.0 * m + 1.0)),
        ] {
            d = 1.0 + aa * d;
            if d.abs() < TINY {
                d = TINY;
            }
            c = 1.0 + aa / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            h *= d * c;
            if (d * c - 1.0).abs() < 1e-15 {
                return h;
            }
        }
    }
    h
}

/// `ln Γ(x)` for x > 0 (Lanczos, g = 7, n = 9).
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        // Reflection; only reached for tiny degrees of freedom.
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + 7.5;
    let s = G[0]
        + G.iter()
            .enumerate()
            .skip(1)
            .map(|(i, g)| g / (x + i as f64))
            .sum::<f64>();
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + s.ln()
}

/// Box plot numbers of the old `BoxPlot::setData`: whiskers are the minimum and maximum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FiveNumbers {
    pub min: f64,
    pub q1: f64,
    pub median: f64,
    pub q3: f64,
    pub max: f64,
}

/// Quartiles by linear interpolation (`gsl_stats_quantile_from_sorted_data`, R `type = 7`).
pub fn five_numbers(values: &[f64]) -> Option<FiveNumbers> {
    let mut x: Vec<f64> = values.iter().copied().filter(|v| !v.is_nan()).collect();
    if x.is_empty() {
        return None;
    }
    x.sort_by(f64::total_cmp);
    let q = |p: f64| {
        let h = p * (x.len() - 1) as f64;
        let i = h.floor() as usize;
        x[i] + (h - i as f64) * (x[(i + 1).min(x.len() - 1)] - x[i])
    };
    Some(FiveNumbers {
        min: x[0],
        q1: q(0.25),
        median: q(0.5),
        q3: q(0.75),
        max: x[x.len() - 1],
    })
}

/// Principal components of the runs, computed as the old program does (`showPCA`).
#[derive(Debug, Clone, PartialEq)]
pub struct Pca {
    /// Columns of the table that are samples (all but the excluded group).
    pub samples: Vec<usize>,
    /// Rows of the table that are variables (finite values, some spread).
    pub rows: Vec<usize>,
    /// Share of the variance of each component, descending, summing to 1.
    pub proportion: Vec<f64>,
    /// `scores[k][i]`: component k of sample i (`samples[i]`).
    pub scores: Vec<Vec<f64>>,
    /// `loadings[k][j]`: weight of variable j (`rows[j]`) in component k; unit length.
    pub loadings: Vec<Vec<f64>>,
}

/// PCA of the z-scored substances over the runs, without the runs of `exclude` (the pool group).
/// With `missing_as_zero` off, a substance with a missing value is left out.
pub fn pca(t: &Table, exclude: &str, missing_as_zero: bool) -> Result<Pca, String> {
    let samples: Vec<usize> = (0..t.runs.len())
        .filter(|c| t.groups[*c] != exclude)
        .collect();
    let n = samples.len();
    if n < 2 {
        return Err("PCA needs at least two runs".into());
    }
    let mut rows = Vec::new();
    let mut z = Vec::new(); // row-major: variable by sample
    for (j, r) in t.rows.iter().enumerate() {
        let x: Option<Vec<f64>> = samples
            .iter()
            .map(|c| r.values[*c].or(missing_as_zero.then_some(0.0)))
            .collect();
        let Some(x) = x.filter(|x| x.iter().all(|v| v.is_finite())) else {
            continue;
        };
        let (m, sd) = (mean(&x), variance(&x).sqrt());
        if sd > 0.0 && sd.is_finite() {
            rows.push(j);
            z.extend(x.iter().map(|v| (v - m) / sd));
        }
    }
    if rows.is_empty() {
        return Err(
            "PCA needs at least one substance with a value in every run and some spread".into(),
        );
    }
    // Samples × variables. Every variable is z-scored, hence centred; the SVD of this matrix has
    // the eigenvectors of the old covariance matrix as right singular vectors.
    let x = DMatrix::from_fn(n, rows.len(), |i, j| z[j * n + i]);
    let svd = x.svd(true, true);
    let (u, vt) = (svd.u.ok_or("SVD failed")?, svd.v_t.ok_or("SVD failed")?);
    let mut order: Vec<usize> = (0..svd.singular_values.len()).collect();
    order.sort_by(|a, b| svd.singular_values[*b].total_cmp(&svd.singular_values[*a]));
    let total: f64 = svd.singular_values.iter().map(|s| s * s).sum();
    let mut out = Pca {
        samples,
        rows,
        proportion: Vec::new(),
        scores: Vec::new(),
        loadings: Vec::new(),
    };
    for k in order {
        let s = svd.singular_values[k];
        let mut load: Vec<f64> = vt.row(k).iter().copied().collect();
        // The sign of a component is arbitrary; make its largest loading positive.
        let big = load
            .iter()
            .copied()
            .max_by(|a, b| a.abs().total_cmp(&b.abs()))
            .unwrap_or(0.0);
        let sign = if big < 0.0 { -1.0 } else { 1.0 };
        load.iter_mut().for_each(|v| *v *= sign);
        out.proportion.push(s * s / total);
        out.scores
            .push(u.column(k).iter().map(|v| sign * s * v).collect());
        out.loadings.push(load);
    }
    Ok(out)
}

/// Writes the comparison as CSV (`.` decimals, `,` separators, readable by R). A missing p-value
/// is an empty cell.
pub fn write_csv(t: &Table, r: &[Biomarker], mut w: impl Write) -> std::io::Result<()> {
    writeln!(
        w,
        "Substance,RT_s,Hits_1,Hits_2,Mean_1,Mean_2,SD_1,SD_2,Ratio,P_value"
    )?;
    for r in r {
        let name = &t.rows[r.row].name;
        let name = match name.contains([',', '"', '\n', '\r']) {
            true => format!("\"{}\"", name.replace('"', "\"\"")),
            false => name.clone(),
        };
        let p = r.p.map(|p| p.to_string()).unwrap_or_default();
        writeln!(
            w,
            "{name},{},{},{},{},{},{},{},{},{p}",
            t.rows[r.row].rt_s,
            r.hits[0],
            r.hits[1],
            r.mean[0],
            r.mean[1],
            r.sd[0],
            r.sd[1],
            r.ratio
        )?;
    }
    w.flush()
}

#[cfg(test)]
// The R reference values are pasted as R prints them (17 digits).
#[allow(clippy::excessive_precision)]
mod tests {
    use super::*;
    use crate::batch::Row;

    /// Reference values: `scripts/r/stats_reference.R` (R 4.5.2, `t.test`, `quantile(type = 7)`,
    /// `prcomp`). Columns: replicate groups c1, c1, c2, c2, d1, d1, d2, d2, $POOL, $POOL.
    fn table() -> Table {
        let s = Some;
        let rows: [[Option<f64>; 10]; 6] = [
            [
                s(107.),
                s(193.8),
                s(93.9),
                s(214.),
                s(100.6),
                s(145.4),
                s(10.2),
                s(120.8),
                s(181.5),
                s(131.9),
            ],
            [
                s(204.8),
                s(112.2),
                None,
                s(202.3),
                s(140.3),
                s(104.5),
                s(185.2),
                s(99.4),
                s(129.1),
                s(167.8),
            ],
            [
                s(190.9),
                s(89.2),
                s(181.2),
                s(141.9),
                s(141.),
                s(214.4),
                s(164.3),
                s(127.3),
                s(203.6),
                s(103.8),
            ],
            [
                s(104.6),
                s(148.6),
                s(112.),
                s(103.2),
                s(159.4),
                s(368.8),
                s(49.7),
                s(112.5),
                s(77.9),
                s(301.),
            ],
            [s(246.8); 10],
            [
                s(174.1),
                s(212.5),
                s(373.1),
                s(126.3),
                s(84.2),
                s(77.9),
                s(149.9),
                s(140.4),
                s(429.5),
                s(207.8),
            ],
        ];
        Table {
            runs: (1..=10).map(|i| format!("run{i}")).collect(),
            groups: [
                "c1", "c1", "c2", "c2", "d1", "d1", "d2", "d2", "$POOL", "$POOL",
            ]
            .map(String::from)
            .to_vec(),
            rows: rows
                .iter()
                .enumerate()
                .map(|(i, v)| Row {
                    compound: i,
                    name: format!("s{i}"),
                    rt_s: 0.0,
                    ions: vec![],
                    values: v.to_vec(),
                    rts_s: vec![None; v.len()],
                })
                .collect(),
            ..Default::default()
        }
    }

    fn groups(samples: Samples, missing_as_zero: bool) -> TwoGroups {
        TwoGroups {
            group1: vec!["c1".into(), "c2".into()],
            group2: vec!["d1".into(), "d2".into()],
            samples,
            missing_as_zero,
        }
    }

    /// Relative difference of 1e-9 or absolute 1e-12.
    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 + 1e-9 * b.abs()
    }

    fn assert_close(got: &[f64], want: &[f64]) {
        assert_eq!(got.len(), want.len());
        for (g, w) in got.iter().zip(want) {
            assert!(close(*g, *w), "{g} != {w} in {got:?}");
        }
    }

    fn ps(samples: Samples, zero: bool) -> Vec<f64> {
        compare(&table(), &groups(samples, zero))
            .iter()
            .map(|r| r.p.unwrap_or(f64::NAN))
            .collect()
    }

    #[test]
    fn summary_follows_the_old_group_columns() {
        // The zero and the NaN are left out of mean and error but a zero would still not be a hit.
        let s = summary(&[2.0, 0.0, 4.0, 6.0, f64::NAN]);
        assert_eq!(s.hits, 3);
        assert_eq!(s.mean, 4.0);
        assert!((s.rel_std_err - 2.0 / 3f64.sqrt() / 4.0).abs() < 1e-12);
        assert!(summary(&[5.0, 0.0]).rel_std_err.is_nan());
    }

    #[test]
    fn student_pools_the_variance() {
        // Hand calculation: means 2 and 5, both variance 1, n = 3: t = 3 / sqrt(1) * sqrt(1.5) = 3.674234
        // with 4 degrees of freedom; R: t.test(c(1,2,3), c(4,5,6), var.equal = TRUE)$p.value.
        let p = student(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]).unwrap();
        assert!((p - 0.02131).abs() < 5e-5, "{p}");
        assert_eq!(student(&[1.0], &[2.0]), None);
    }

    #[test]
    fn welch_matches_r() {
        assert_close(
            &ps(Samples::Independent, true),
            &[
                0.2194457641460999,
                0.96378147903546618,
                0.72807352761907573,
                0.48320004744279366,
                1.0,
                0.13374620164900489,
            ],
        );
        // The missing value is left out: 3 samples against 4.
        assert_close(
            &ps(Samples::Independent, false),
            &[
                0.2194457641460999,
                0.33099226783494695,
                0.72807352761907573,
                0.48320004744279366,
                1.0,
                0.13374620164900489,
            ],
        );
        assert_close(
            &ps(Samples::GroupMeans, true),
            &[
                0.29230840320306134,
                0.94508866368773414,
                0.63352696366480421,
                0.65268479927156542,
                1.0,
                0.12835383220085064,
            ],
        );
    }

    #[test]
    fn paired_matches_r() {
        let p = ps(Samples::Paired, true);
        assert_close(
            &[p[0], p[1], p[2], p[3], p[5]],
            &[
                0.060448931607732896,
                0.97098702627560962,
                0.79677240894854806,
                0.42289479817094922,
                0.11549841945139186,
            ],
        );
        // All differences are 0: the test is undefined.
        assert!(p[4].is_nan());
    }

    #[test]
    fn paired_follows_the_group_order_and_needs_matching_pairs() {
        let t = table();
        let two = |g1: &[&str], g2: &[&str]| TwoGroups {
            group1: g1.iter().map(|g| g.to_string()).collect(),
            group2: g2.iter().map(|g| g.to_string()).collect(),
            samples: Samples::Paired,
            missing_as_zero: true,
        };
        let same = two(&["c1", "c2"], &["d1", "d2"]);
        let swapped = two(&["c1", "c2"], &["d2", "d1"]);
        assert!(pairs_match(&t, &same) && pairs_match(&t, &swapped));
        // Pairing c1 with d2 instead of d1 changes the differences, so the p-values.
        let (a, b) = (compare(&t, &same), compare(&t, &swapped));
        assert!(a.iter().zip(&b).any(|(a, b)| a.p != b.p));
        // One group against two, or sides with different group counts: no pairs.
        let uneven = two(&["c1"], &["d1", "d2"]);
        assert!(!pairs_match(&t, &uneven));
        assert!(compare(&t, &uneven).iter().all(|r| r.p.is_none()));
    }

    #[test]
    fn means_deviations_hits_and_ratio() {
        let r = compare(&table(), &groups(Samples::Independent, true));
        let mean1: Vec<f64> = r.iter().map(|r| r.mean[0]).collect();
        let mean2: Vec<f64> = r.iter().map(|r| r.mean[1]).collect();
        assert_close(
            &mean1,
            &[
                152.17500000000001,
                129.82500000000002,
                150.80000000000001,
                117.09999999999999,
                246.80000000000001,
                221.5,
            ],
        );
        assert_close(
            &mean2,
            &[
                94.25,
                132.34999999999999,
                161.75,
                172.59999999999999,
                246.80000000000001,
                113.10000000000001,
            ],
        );
        assert_close(
            &r.iter().map(|r| r.sd[0]).collect::<Vec<_>>(),
            &[
                60.530233492583413,
                96.676518176166581,
                46.209450692832661,
                21.351971025333157,
                0.0,
                107.04105131521582,
            ],
        );
        assert_eq!(r[1].hits, [3, 4]);
        assert!(close(r[0].ratio, 152.17500000000001 / 94.25));
    }

    #[test]
    fn box_numbers_match_r() {
        let t = table();
        let x = group_values(&t, 0, &["c1".into(), "c2".into()], true);
        let f = five_numbers(&x).unwrap();
        assert_close(
            &[f.min, f.q1, f.median, f.q3, f.max],
            &[93.9, 103.72499999999999, 150.4, 198.85000000000002, 214.0],
        );
        assert_eq!(five_numbers(&[]), None);
        assert_eq!(five_numbers(&[7.0]).map(|f| f.median), Some(7.0));
    }

    #[test]
    fn pca_matches_r() {
        let p = pca(&table(), "$POOL", true).unwrap();
        // The constant substance 4 is dropped; the pool runs are not samples.
        assert_eq!(p.rows, [0, 1, 2, 3, 5]);
        assert_eq!(p.samples, (0..8).collect::<Vec<_>>());
        assert_close(
            &p.proportion[..5],
            &[
                0.34010757483510873,
                0.30119974929897469,
                0.26290499147371144,
                0.08576972355885841,
                0.010017960833346828,
            ],
        );
        assert!((p.proportion.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        let abs = |v: &[f64]| v.iter().map(|x| x.abs()).collect::<Vec<_>>();
        assert_close(
            &abs(&p.scores[0]),
            &[
                0.35326950399601542,
                0.23518043691699736,
                2.8308304317185184,
                1.1134597691280421,
                0.64701290598949956,
                1.3968675557293047,
                0.25026015174465432,
                0.19433871446269166,
            ],
        );
        assert_close(
            &abs(&p.scores[1]),
            &[
                0.14067945246110566,
                1.000756619799239,
                0.98903356389110797,
                1.0371164999188829,
                0.17864729242478034,
                2.5830738563805031,
                0.71859128978178366,
                0.4963162658858199,
            ],
        );
        assert_close(
            &abs(&p.loadings[0]),
            &[
                0.2741283162169148,
                0.54305704799293886,
                0.034635453030580243,
                0.33486373977687417,
                0.71875543131811304,
            ],
        );
        assert_close(
            &abs(&p.loadings[1]),
            &[
                0.071453698328864612,
                0.39991652604431982,
                0.62886713902947577,
                0.66286709436598068,
                0.0097199715769222678,
            ],
        );
        // The sign is fixed: the largest loading of each component is positive.
        assert!(p.loadings[0][4] > 0.0 && p.loadings[1][3] > 0.0);
    }

    #[test]
    fn pca_needs_data() {
        let mut t = table();
        assert!(
            pca(&t, "$POOL", false).is_ok(),
            "row 1 has a missing value and is left out"
        );
        t.groups = vec!["x".into(); 10];
        t.rows.retain(|r| r.name == "s4");
        assert!(pca(&t, "", true).is_err());
        t.rows.clear();
        assert!(pca(&t, "", true).is_err());
    }

    #[test]
    fn small_samples_and_csv() {
        assert_eq!(welch(&[1.0], &[1.0, 2.0]), 1.0);
        assert_eq!(welch(&[1.0, 1.0], &[1.0, 2.0]), 1.0);
        // t = 0 gives p = 1; a huge difference gives p near 0.
        assert!(close(welch(&[1.0, 2.0, 3.0], &[1.0, 2.0, 3.0]), 1.0));
        assert!((welch(&[1.0, 1.1, 0.9], &[100.0, 101.0, 99.0]) - 2.87141e-5).abs() < 1e-10);
        let t = table();
        let r = compare(&t, &groups(Samples::Independent, true));
        let mut out = Vec::new();
        write_csv(&t, &r, &mut out).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.starts_with(
            "Substance,RT_s,Hits_1,Hits_2,Mean_1,Mean_2,SD_1,SD_2,Ratio,P_value\ns0,0,4,4,"
        ));
        assert_eq!(out.lines().count(), 7);
    }
}
