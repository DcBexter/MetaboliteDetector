//! Normalisation of the batch table (F13, `docs/reverse-engineering/algorithms.md` §7c): internal
//! standard, total signal, pool measurement and z-score. Pure functions on [`Table`]. A value
//! becomes `None` where the old program divides by zero (it yields inf/NaN there).

use crate::batch::Table;

/// How the summed quant ion integrals of each run are put on a common scale.
#[derive(Debug, Clone, PartialEq)]
pub enum Normalization {
    /// Raw integrals.
    None,
    /// Divide each run by the value of one substance (index into `Table::rows`) in that run. That
    /// substance's own row becomes 1.
    InternalStandard(usize),
    /// Divide each run by the sum of all substance values in that run.
    TotalSignal,
    /// Divide each substance in a run by its mean or median in the nearest pool runs.
    Pool(Pool),
}

/// Pool measurement settings. Pool runs are the columns whose replicate group is `group`; they keep
/// their raw values.
#[derive(Debug, Clone, PartialEq)]
pub struct Pool {
    /// Replicate group of the pool runs (the old program uses [`POOL_GROUP`]).
    pub group: String,
    /// Acquisition time of each column of the table; only differences matter.
    pub times: Vec<i64>,
    /// Number of nearest pool runs per run.
    pub nearest: usize,
    /// Median instead of mean of those pool runs.
    pub median: bool,
}

pub const POOL_GROUP: &str = "$POOL";

/// The table with `mode` applied; the shape stays the same.
pub fn normalize(t: &Table, mode: &Normalization) -> Table {
    let mut out = t.clone();
    match mode {
        Normalization::None => {}
        Normalization::InternalStandard(row) => {
            let standard = t.rows.get(*row).map(|r| r.values.clone());
            for r in &mut out.rows {
                for (c, v) in r.values.iter_mut().enumerate() {
                    let d = standard.as_ref().and_then(|s| s[c]);
                    *v = divide(*v, d);
                }
            }
        }
        Normalization::TotalSignal => {
            // A missing value counts as 0, like the old table.
            let total: Vec<f64> = (0..t.runs.len())
                .map(|c| t.rows.iter().filter_map(|r| r.values[c]).sum())
                .collect();
            for r in &mut out.rows {
                for (c, v) in r.values.iter_mut().enumerate() {
                    *v = divide(*v, Some(total[c]));
                }
            }
        }
        Normalization::Pool(p) => {
            let near = nearest_pools(t, p);
            for (r, o) in t.rows.iter().zip(&mut out.rows) {
                for (c, v) in o.values.iter_mut().enumerate() {
                    if let Some(near) = &near[c] {
                        let pool: Vec<f64> = near.iter().filter_map(|i| r.values[*i]).collect();
                        *v = divide(*v, centre(pool, p.median));
                    }
                }
            }
        }
    }
    out
}

/// For each column, the columns of its `p.nearest` closest pool runs by acquisition time (a pool
/// run counts itself, at distance 0); a tie goes to the run scanned earlier like the old multimap,
/// so moving columns does not change the result.
pub fn pool_mapping(t: &Table, p: &Pool) -> Vec<Vec<usize>> {
    let at = |i: usize| p.times.get(i).copied().unwrap_or(0);
    let scanned = |i: usize| t.sources.get(i).copied().unwrap_or(i);
    let pools: Vec<usize> = (0..t.runs.len())
        .filter(|c| t.groups[*c] == p.group)
        .collect();
    (0..t.runs.len())
        .map(|c| {
            let mut near = pools.clone();
            near.sort_by_key(|i| (at(*i).abs_diff(at(c)), scanned(*i)));
            near.truncate(p.nearest);
            near
        })
        .collect()
}

/// The old program's "Export pool mapping" file (`exportPoolMapping` 0x4552a0): a header, then one
/// line per run with the names of the pool runs its values were divided by, tab-separated.
pub fn pool_mapping_text(t: &Table, p: &Pool) -> String {
    let mut out = String::from("Mapping\nSample\tPools\n\n");
    for (name, near) in t.runs.iter().zip(pool_mapping(t, p)) {
        out += name;
        for i in near {
            out += &format!("\t{}", t.runs[i]);
        }
        out.push('\n');
    }
    out
}

/// For each column: `None` for a pool run (kept raw), else its nearest pool runs.
fn nearest_pools(t: &Table, p: &Pool) -> Vec<Option<Vec<usize>>> {
    (pool_mapping(t, p).into_iter().enumerate())
        .map(|(c, near)| (t.groups[c] != p.group).then_some(near))
        .collect()
}

/// Mean or median (GSL: the middle of an odd count, else the mean of the two middle values).
/// shortcut: a missing value in a pool run is skipped here, the old program counted it as 0.
fn centre(mut v: Vec<f64>, median: bool) -> Option<f64> {
    if v.is_empty() {
        return None;
    }
    if !median {
        return Some(v.iter().sum::<f64>() / v.len() as f64);
    }
    v.sort_by(f64::total_cmp);
    let n = v.len();
    Some(if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    })
}

fn divide(v: Option<f64>, d: Option<f64>) -> Option<f64> {
    v.zip(d).and_then(|(v, d)| (d != 0.0).then_some(v / d))
}

/// Z-score of each substance over the runs where it has a value: `(x - mean) / sd` with the sample
/// standard deviation (`tools::Filter::zScoreTransformation`, N-1). A substance with fewer than two
/// values or no spread becomes all `None` (the old program yields NaN/inf).
pub fn z_score(t: &Table) -> Table {
    let mut out = t.clone();
    for r in &mut out.rows {
        let x: Vec<f64> = r.values.iter().flatten().copied().collect();
        let n = x.len() as f64;
        let mean = x.iter().sum::<f64>() / n;
        let sd = (x.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        for v in &mut r.values {
            *v = v.and_then(|v| (sd > 0.0).then(|| (v - mean) / sd));
        }
    }
    out
}

/// The normalised table, then optionally its z-scores: the order of the old `getTableData`.
pub fn apply(t: &Table, mode: &Normalization, z: bool) -> Table {
    let t = normalize(t, mode);
    if z { z_score(&t) } else { t }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::Row;

    fn row(name: &str, values: &[Option<f64>]) -> Row {
        Row {
            compound: 0,
            name: name.into(),
            rt_s: 0.0,
            ions: vec![],
            values: values.to_vec(),
            rts_s: vec![None; values.len()],
        }
    }

    /// Columns: a1, a2, pool1, b1, pool2 (group order is irrelevant here).
    fn table() -> Table {
        let s = |v: f64| Some(v);
        Table {
            runs: ["a1", "a2", "p1", "b1", "p2"].map(String::from).to_vec(),
            groups: ["A", "A", "$POOL", "B", "$POOL"].map(String::from).to_vec(),
            rows: vec![
                row("std", &[s(10.0), s(20.0), s(10.0), s(40.0), s(30.0)]),
                row("x", &[s(5.0), None, s(20.0), s(20.0), s(40.0)]),
            ],
            ..Default::default()
        }
    }

    fn values(t: &Table, row: usize) -> Vec<Option<f64>> {
        t.rows[row].values.clone()
    }

    #[test]
    fn internal_standard_divides_by_the_standard_row() {
        let n = normalize(&table(), &Normalization::InternalStandard(0));
        assert_eq!(values(&n, 0), vec![Some(1.0); 5]);
        // 5/10, missing, 20/10, 20/40, 40/30
        assert_eq!(
            values(&n, 1),
            vec![Some(0.5), None, Some(2.0), Some(0.5), Some(40.0 / 30.0)]
        );
    }

    #[test]
    fn standard_that_is_missing_or_zero_gives_none() {
        let mut t = table();
        t.rows[0].values[0] = None;
        t.rows[0].values[1] = Some(0.0);
        let n = normalize(&t, &Normalization::InternalStandard(0));
        assert_eq!(values(&n, 1)[..2], [None, None]);
        assert_eq!(
            values(&normalize(&t, &Normalization::InternalStandard(9)), 1),
            vec![None; 5]
        );
    }

    #[test]
    fn total_signal_divides_by_the_column_sum() {
        let n = normalize(&table(), &Normalization::TotalSignal);
        // Column sums: 15, 20, 30, 60, 70.
        assert_eq!(values(&n, 0)[0], Some(10.0 / 15.0));
        assert_eq!(values(&n, 1)[0], Some(5.0 / 15.0));
        assert_eq!(values(&n, 0)[1], Some(1.0));
        assert_eq!(values(&n, 1)[1], None);
        assert_eq!(values(&n, 1)[4], Some(40.0 / 70.0));
    }

    #[test]
    fn pool_uses_the_nearest_pool_runs() {
        let pool = |nearest, median| {
            Normalization::Pool(Pool {
                group: POOL_GROUP.into(),
                times: vec![0, 10, 2, 8, 9],
                nearest,
                median,
            })
        };
        // a1 (t=0): nearest pool p1 (t=2); a2 (t=10): p2 (t=9); b1 (t=8): p2 (t=9).
        let n = normalize(&table(), &pool(1, false));
        assert_eq!(
            values(&n, 0),
            vec![
                Some(1.0),
                Some(20.0 / 30.0),
                Some(10.0),
                Some(40.0 / 30.0),
                Some(30.0)
            ]
        );
        // Pool runs keep their raw values; x is missing in a2.
        assert_eq!(
            values(&n, 1),
            vec![
                Some(5.0 / 20.0),
                None,
                Some(20.0),
                Some(20.0 / 40.0),
                Some(40.0)
            ]
        );
        // Both pool runs: x has 20 and 40, mean 30, median 30; std has 10 and 30, mean = median = 20.
        let n = normalize(&table(), &pool(2, true));
        assert_eq!(values(&n, 0)[0], Some(10.0 / 20.0));
        assert_eq!(values(&n, 1)[0], Some(5.0 / 30.0));
        // Three values in the reference: median differs from the mean.
        let mut t = table();
        t.runs.push("p3".into());
        t.groups.push("$POOL".into());
        t.rows[0].values.push(Some(100.0));
        t.rows[1].values.push(Some(0.0));
        let mut p = Pool {
            group: POOL_GROUP.into(),
            times: vec![0, 10, 2, 8, 9, 5],
            nearest: 3,
            median: true,
        };
        assert_eq!(
            values(&normalize(&t, &Normalization::Pool(p.clone())), 0)[0],
            Some(10.0 / 30.0)
        );
        p.median = false;
        assert_eq!(
            values(&normalize(&t, &Normalization::Pool(p)), 0)[0],
            Some(10.0 / (140.0 / 3.0))
        );
    }

    #[test]
    fn pool_ties_do_not_depend_on_the_column_order() {
        let t = Table {
            runs: ["p1", "p2", "x"].map(String::from).to_vec(),
            sources: vec![0, 1, 2],
            groups: ["$POOL", "$POOL", "A"].map(String::from).to_vec(),
            rows: vec![row("s", &[Some(1.0), Some(2.0), Some(10.0)])],
        };
        let pool = Normalization::Pool(Pool {
            group: POOL_GROUP.into(),
            times: vec![0; 3],
            nearest: 1,
            median: false,
        });
        assert_eq!(values(&normalize(&t, &pool), 0)[2], Some(10.0));
        // Pool runs swapped, as after the user moved a run: same result.
        let moved = t.select_runs(&[1, 0, 2]);
        assert_eq!(values(&normalize(&moved, &pool), 0)[2], Some(10.0));
    }

    #[test]
    fn the_pool_mapping_file_lists_the_pools_each_value_was_divided_by() {
        let p = Pool {
            group: POOL_GROUP.into(),
            times: vec![0, 10, 2, 8, 9],
            nearest: 1,
            median: false,
        };
        // a1 (t=0): p1; a2 (t=10): p2; p1 itself; b1 (t=8): p2; p2 itself.
        assert_eq!(
            pool_mapping_text(&table(), &p),
            "Mapping\nSample\tPools\n\na1\tp1\na2\tp2\np1\tp1\nb1\tp2\np2\tp2\n"
        );
        let p = Pool { nearest: 2, ..p };
        assert!(pool_mapping_text(&table(), &p).contains("\na2\tp2\tp1\n"));
        // The mapping is what `normalize` divides by: a2's x is missing, a1's is 5 / x(p1) = 5 / 20.
        let n = normalize(&table(), &Normalization::Pool(p));
        assert_eq!(values(&n, 1)[0], Some(5.0 / 30.0));
    }

    #[test]
    fn z_score_uses_the_sample_deviation() {
        // 2, 4, 6: mean 4, variance (4 + 0 + 4) / 2 = 4, sd 2.
        let t = Table {
            runs: vec!["a".into(), "b".into(), "c".into(), "d".into()],
            groups: vec![String::new(); 4],
            rows: vec![
                row("x", &[Some(2.0), Some(4.0), None, Some(6.0)]),
                row("flat", &[Some(1.0), Some(1.0), Some(1.0), Some(1.0)]),
                row("single", &[None, Some(1.0), None, None]),
            ],
            ..Default::default()
        };
        let z = z_score(&t);
        assert_eq!(values(&z, 0), vec![Some(-1.0), Some(0.0), None, Some(1.0)]);
        assert_eq!(values(&z, 1), vec![None; 4]);
        assert_eq!(values(&z, 2), vec![None; 4]);
    }

    #[test]
    fn apply_normalises_before_the_z_score() {
        let t = table();
        let a = apply(&t, &Normalization::InternalStandard(0), true);
        assert_eq!(
            a,
            z_score(&normalize(&t, &Normalization::InternalStandard(0)))
        );
        assert_eq!(apply(&t, &Normalization::None, false), t);
    }
}
