//! Edits of the result table the batch window offers: keep, drop or reorder runs, drop
//! substances, and sums per run. Pure functions on [`Table`]; they leave the cells untouched.

use super::Table;

impl Table {
    /// The table with only the columns `cols`, in that order (drop = leave out, reorder = permute).
    /// Out-of-range indices are ignored.
    pub fn select_runs(&self, cols: &[usize]) -> Table {
        let cols: Vec<usize> = cols
            .iter()
            .copied()
            .filter(|c| *c < self.runs.len())
            .collect();
        let pick = |v: &Vec<String>| cols.iter().map(|c| v[*c].clone()).collect();
        let pick_cells = |v: &Vec<Option<f64>>| cols.iter().map(|c| v[*c]).collect();
        Table {
            runs: pick(&self.runs),
            groups: pick(&self.groups),
            rows: (self.rows.iter())
                .map(|r| super::Row {
                    values: pick_cells(&r.values),
                    rts_s: pick_cells(&r.rts_s),
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
    use super::super::Row;
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
}
