//! Edits of the result table the batch window offers: keep, drop or reorder runs, drop
//! substances, and sums per run. Pure functions on [`Table`]; they leave the cells untouched.

use std::borrow::Cow;
use std::collections::HashMap;

use super::{RunResult, Table};

/// Per run: library index of a substance -> position in `RunResult::found`, so the details behind
/// a cell are found without searching every frame.
pub fn hit_index(runs: &[RunResult]) -> Vec<HashMap<usize, usize>> {
    (runs.iter())
        .map(|r| {
            r.found
                .iter()
                .enumerate()
                .map(|(i, f)| (f.compound, i))
                .collect()
        })
        .collect()
}

/// Picks `v[i]` for each `i` of `idx`, in that order. The one helper that moves run columns, shared by
/// [`Table::select_runs`] and the MID report's run moves; it moves values and computes nothing, so
/// no statistic changes with the order.
pub(crate) fn pick<T: Clone>(v: &[T], idx: &[usize]) -> Vec<T> {
    idx.iter().map(|i| v[*i].clone()).collect()
}

/// Columns and rows removed or moved in the shown table. Runs are edited before the
/// normalisation (a removed run takes no part in z-scores), substances after it.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Edits {
    /// Columns of the table that are kept, in display order; `None` = all.
    pub cols: Option<Vec<usize>>,
    /// Library indices of the removed substances.
    pub rows: Vec<usize>,
}

/// What a click in the shown table asks for; columns are numbered as shown.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Act {
    DropRun(usize),
    MoveRun(usize, usize),
    DropRow(usize),
}

impl Edits {
    pub fn cols(&self, n: usize) -> Vec<usize> {
        self.cols.clone().unwrap_or_else(|| (0..n).collect())
    }

    pub fn is_empty(&self) -> bool {
        self.cols.is_none() && self.rows.is_empty()
    }

    /// "Restore removed": brings back all runs and substances and the original order.
    pub fn restore(&mut self) {
        *self = Self::default();
    }

    /// Applies what the user asked for in the shown table with `n` columns in the base table.
    pub fn apply(&mut self, act: Act, n: usize) {
        let mut cols = self.cols(n);
        match act {
            Act::DropRun(j) if j < cols.len() => {
                cols.remove(j);
            }
            Act::MoveRun(j, to) if j < cols.len() && to < cols.len() => cols.swap(j, to),
            Act::DropRow(compound) => return self.rows.push(compound),
            _ => return,
        }
        self.cols = Some(cols);
    }
}

/// The runs as the user left them (removed and moved), with the file dates that follow them.
pub fn edited_runs<'a>(
    t: &'a Table,
    edits: &Edits,
    times: Option<&[i64]>,
) -> (Cow<'a, Table>, Vec<i64>) {
    let cols = edits.cols.as_deref();
    let base = cols.map_or(Cow::Borrowed(t), |c| Cow::Owned(t.select_runs(c)));
    // The file dates follow the runs that are left.
    let times = match (times, cols) {
        (Some(ts), Some(cols)) => cols
            .iter()
            .map(|c| ts.get(*c).copied().unwrap_or(0))
            .collect(),
        (Some(ts), None) => ts.to_vec(),
        (None, _) => Vec::new(),
    };
    (base, times)
}

impl Table {
    /// The table with only the columns `cols`, in that order (drop = leave out, reorder = permute).
    /// Out-of-range indices are ignored.
    pub fn select_runs(&self, cols: &[usize]) -> Table {
        let cols: Vec<usize> = cols
            .iter()
            .copied()
            .filter(|c| *c < self.runs.len())
            .collect();
        Table {
            runs: pick(&self.runs, &cols),
            sources: (cols.iter())
                .filter_map(|c| self.sources.get(*c).copied())
                .collect(),
            groups: pick(&self.groups, &cols),
            rows: (self.rows.iter())
                .map(|r| super::Row {
                    values: pick(&r.values, &cols),
                    rts_s: pick(&r.rts_s, &cols),
                    ..r.clone()
                })
                .collect(),
        }
    }

    /// The table without the substances whose library index is in `compounds`.
    pub fn without_rows(&self, compounds: &[usize]) -> Table {
        Table {
            rows: (self.rows.iter())
                .filter(|r| !compounds.contains(&r.compound))
                .cloned()
                .collect(),
            ..self.clone()
        }
    }

    /// Sum of the values in each column; a missing value counts as 0 (the old "total sample
    /// intensities").
    pub fn total_intensities(&self) -> Vec<f64> {
        (0..self.runs.len())
            .map(|c| self.rows.iter().filter_map(|r| r.values[c]).sum())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::super::{Found, Row};
    use super::*;

    fn table() -> Table {
        let row = |compound, v: [Option<f64>; 3]| Row {
            compound,
            name: format!("s{compound}"),
            rt_s: 60.0,
            ions: vec![73],
            values: v.to_vec(),
            rts_s: v.map(|v| v.map(|_| 60.0 + compound as f64)).to_vec(),
        };
        Table {
            runs: vec!["a".into(), "b".into(), "c".into()],
            sources: vec![0, 1, 2],
            groups: vec!["x".into(), "x".into(), "y".into()],
            rows: vec![
                row(4, [Some(1.0), None, Some(3.0)]),
                row(7, [Some(10.0), Some(20.0), None]),
            ],
        }
    }

    #[test]
    fn runs_can_be_dropped_and_reordered_with_their_cells() {
        let t = table().select_runs(&[2, 0, 9]);
        assert_eq!(t.runs, ["c", "a"]);
        assert_eq!(t.groups, ["y", "x"]);
        assert_eq!(t.sources, [2, 0]);
        assert_eq!(t.rows[0].values, [Some(3.0), Some(1.0)]);
        assert_eq!(t.rows[1].values, [None, Some(10.0)]);
        assert_eq!(t.rows[1].rts_s, [None, Some(67.0)]);
    }

    #[test]
    fn rows_are_dropped_by_library_index() {
        let t = table().without_rows(&[4, 99]);
        assert_eq!(t.rows.len(), 1);
        assert_eq!(t.rows[0].compound, 7);
    }

    #[test]
    fn totals_count_a_missing_value_as_zero() {
        assert_eq!(table().total_intensities(), [11.0, 20.0, 3.0]);
        assert_eq!(table().select_runs(&[]).total_intensities(), [0.0; 0]);
    }

    #[test]
    fn edits_drop_move_and_restore_columns_and_rows() {
        let mut e = Edits::default();
        e.apply(Act::MoveRun(0, 1), 3);
        assert_eq!(e.cols, Some(vec![1, 0, 2]));
        e.apply(Act::DropRun(1), 3);
        assert_eq!(e.cols, Some(vec![1, 2]));
        // Out of range clicks (a stale frame) change nothing.
        e.apply(Act::DropRun(5), 3);
        e.apply(Act::MoveRun(1, 2), 3);
        assert_eq!(e.cols, Some(vec![1, 2]));
        e.apply(Act::DropRow(7), 3);
        assert_eq!(e.rows, [7]);
        assert!(!e.is_empty());
        e.restore();
        assert!(e.is_empty());
    }

    #[test]
    fn the_file_dates_follow_the_runs_that_are_left() {
        let t = table();
        let times = [10, 20, 30];
        let (base, dates) = edited_runs(&t, &Edits::default(), Some(&times[..]));
        assert!(matches!(base, Cow::Borrowed(_)));
        assert_eq!((base.runs.len(), dates), (3, vec![10, 20, 30]));
        let edits = Edits {
            cols: Some(vec![2, 0]),
            rows: Vec::new(),
        };
        let (base, dates) = edited_runs(&t, &edits, Some(&times[..]));
        assert_eq!(base.runs, ["c", "a"]);
        assert_eq!(dates, [30, 10]);
        assert!(edited_runs(&t, &edits, None).1.is_empty());
    }

    #[test]
    fn a_run_move_picks_the_values_in_the_new_order() {
        assert_eq!(pick(&[1.5, 2.5, 3.5], &[2, 0]), [3.5, 1.5]);
        assert_eq!(pick(&["a", "b"], &[1, 1]), ["b", "b"]);
    }

    #[test]
    fn hit_index_maps_each_substance_to_its_hit_in_every_run() {
        let found = |compound| Found {
            compound,
            score: 0.9,
            rt_s: 60.0,
            ions: Vec::new(),
            rescanned: false,
            rescued: Vec::new(),
            alternatives: Vec::new(),
        };
        let run = |found| RunResult {
            input: 0,
            path: std::path::PathBuf::new(),
            found,
            unknown: Vec::new(),
            calibration: None,
        };
        let index = hit_index(&[run(vec![found(5), found(9)]), run(Vec::new())]);
        assert_eq!(index.len(), 2);
        assert_eq!(index[0], HashMap::from([(5, 0), (9, 1)]));
        assert!(index[1].is_empty());
    }
}
