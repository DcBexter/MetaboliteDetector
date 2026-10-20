//! The statistics of the replicate groups per substance, as the "Group statistics" view shows them
//! and as `Table::write_group_stats_csv` writes them.

use std::io::{self, Write};
use std::path::{Path, PathBuf};

use super::{BOM, Table, text_cell};
use crate::stats::{Summary, anova, group_values, student, summary};

/// Where the group statistics of a results CSV go: `results.csv` → `results.groupstats.csv`.
pub fn group_stats_path(csv: &Path) -> PathBuf {
    csv.with_extension("groupstats.csv")
}

/// The sorted names of the replicate groups of `t` (a run without a group has the empty name).
pub fn group_names(t: &Table) -> Vec<&str> {
    let mut g: Vec<&str> = t.groups.iter().map(String::as_str).collect();
    g.sort_unstable();
    g.dedup();
    g
}

/// The name of a group in a heading: a run without a group belongs to the whole batch.
pub fn group_label(g: &str) -> &str {
    if g.is_empty() { "all runs" } else { g }
}

/// One substance in the groups, a missing value counting as 0 (old `BatchQuantificationStatisticsTableModel`).
#[derive(Debug, Clone, PartialEq)]
pub struct GroupStats {
    /// Per group, in the order given.
    pub groups: Vec<Summary>,
    /// Per group the p-value of its values against those of the blank group (`student`); `None`
    /// for the blank group itself and where no test is defined. Empty without a blank group.
    pub p_blank: Vec<Option<f64>>,
    /// The ANOVA p-value over all groups, blank included.
    pub anova: Option<f64>,
}

impl Table {
    /// The statistics of substance `row` in `groups`; `blank` names the group to test against.
    pub fn group_stats(&self, row: usize, groups: &[&str], blank: Option<&str>) -> GroupStats {
        let values = |g: &str| group_values(self, row, &[g.to_string()], true);
        let all: Vec<Vec<f64>> = groups.iter().map(|g| values(g)).collect();
        let blanks = blank.map(values);
        let p_blank = match &blanks {
            Some(b) => (groups.iter().zip(&all))
                .map(|(g, x)| (blank != Some(g)).then(|| student(x, b)).flatten())
                .collect(),
            None => Vec::new(),
        };
        GroupStats {
            groups: all.iter().map(|x| summary(x)).collect(),
            p_blank,
            anova: anova(&all),
        }
    }

    /// The group statistics of every substance as UTF-8 CSV like `write_csv`:
    /// `Substance,RT_s,Quant_ions`, then `Mean_<g>`, `Rel_Std_Err_<g>`, `Hits_<g>` and with a blank
    /// group `P_blank_<g>` for every group, last `P_ANOVA`. Values are not rounded; a value that is
    /// not defined is an empty cell (R reads `NA`). A run without a group is the group `all runs`.
    pub fn write_group_stats_csv(&self, blank: Option<&str>, mut w: impl Write) -> io::Result<()> {
        let groups = group_names(self);
        let blank = blank.filter(|b| groups.contains(b));
        let mut head = ["Substance", "RT_s", "Quant_ions"]
            .map(String::from)
            .to_vec();
        for g in &groups {
            let stats = ["Mean", "Rel_Std_Err", "Hits"];
            let blank_stat = blank.map(|_| "P_blank");
            head.extend(
                stats
                    .into_iter()
                    .chain(blank_stat)
                    .map(|s| text_cell(&format!("{s}_{}", group_label(g)))),
            );
        }
        head.push("P_ANOVA".into());
        w.write_all(BOM.as_bytes())?;
        writeln!(w, "{}", head.join(","))?;
        let number = |v: f64| {
            if v.is_finite() {
                v.to_string()
            } else {
                String::new()
            }
        };
        let optional = |v: Option<f64>| v.map_or(String::new(), number);
        for (i, r) in self.rows.iter().enumerate() {
            let s = self.group_stats(i, &groups, blank);
            let ions: Vec<String> = r.ions.iter().map(u32::to_string).collect();
            let mut cells = vec![text_cell(&r.name), r.rt_s.to_string(), ions.join(" ")];
            for (k, g) in s.groups.iter().enumerate() {
                cells.extend([number(g.mean), number(g.rel_std_err), g.hits.to_string()]);
                if blank.is_some() {
                    cells.push(optional(s.p_blank[k]));
                }
            }
            cells.push(optional(s.anova));
            writeln!(w, "{}", cells.join(","))?;
        }
        w.flush()
    }

    /// Writes the group statistics next to the results CSV `csv` (`group_stats_path`) and returns
    /// that file; without replicate groups an old file of that name is removed and there is none.
    pub fn write_group_stats_file(
        &self,
        csv: &Path,
        blank: Option<&str>,
    ) -> io::Result<Option<PathBuf>> {
        let path = group_stats_path(csv);
        if !self.has_groups() {
            return match std::fs::remove_file(&path) {
                Err(e) if e.kind() != io::ErrorKind::NotFound => Err(e),
                _ => Ok(None),
            };
        }
        let f = io::BufWriter::new(std::fs::File::create(&path)?);
        self.write_group_stats_csv(blank, f)?;
        Ok(Some(path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::Row;

    fn table() -> Table {
        // Groups A (1 2 3), B (2 3 4) and BLANK (5 6 7): ANOVA p = 27 / 4096 (see stats::anova).
        let v = [1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        let row = |name: &str, values: Vec<Option<f64>>| Row {
            compound: 0,
            name: name.into(),
            rt_s: 60.0,
            ions: vec![73, 147],
            values,
            rts_s: vec![Some(60.0); 9],
        };
        Table {
            runs: (0..9).map(|i| format!("r{i}")).collect(),
            groups: ["A", "A", "A", "B", "B", "B", "BLANK", "BLANK", "BLANK"]
                .map(String::from)
                .to_vec(),
            rows: vec![row("one", v.map(Some).to_vec()), row("none", vec![None; 9])],
            ..Default::default()
        }
    }

    #[test]
    fn the_csv_has_the_numbers_of_the_view() {
        let t = table();
        let mut out = Vec::new();
        t.write_group_stats_csv(Some("BLANK"), &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(
            lines[0],
            "\u{feff}Substance,RT_s,Quant_ions,Mean_A,Rel_Std_Err_A,Hits_A,P_blank_A,\
             Mean_B,Rel_Std_Err_B,Hits_B,P_blank_B,\
             Mean_BLANK,Rel_Std_Err_BLANK,Hits_BLANK,P_blank_BLANK,P_ANOVA"
        );
        let c: Vec<&str> = lines[1].split(',').collect();
        assert_eq!(c[..4], ["one", "60", "73 147", "2"]);
        assert_eq!(c[5], "3"); // hits of A
        assert_eq!(c[14], ""); // the blank group has no p-value against itself
        let p: f64 = c[15].parse().unwrap();
        assert!((p - 27.0 / 4096.0).abs() < 1e-12, "{p}");
        // No value at all: the means are 0, the spread and the tests are not defined.
        let c: Vec<&str> = lines[2].split(',').collect();
        assert_eq!(c[3..7], ["0", "", "0", ""]);
        assert_eq!(c[15], "");
    }

    #[test]
    fn without_a_blank_group_there_is_no_blank_column() {
        let t = table();
        let mut out = Vec::new();
        t.write_group_stats_csv(Some("water"), &mut out).unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(!text.contains("P_blank"), "{text}");
        assert!(text.lines().next().unwrap().ends_with(",P_ANOVA"));
    }

    #[test]
    fn the_file_is_written_with_groups_only() {
        let dir = std::env::temp_dir().join(format!("md-core-gstats-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let csv = dir.join("r.csv");
        let t = table();
        let p = t.write_group_stats_file(&csv, None).unwrap().unwrap();
        assert_eq!(p, dir.join("r.groupstats.csv"));
        let ungrouped = Table {
            groups: vec![String::new(); 9],
            ..t
        };
        // An old file would be read by a script as the statistics of this table.
        assert_eq!(ungrouped.write_group_stats_file(&csv, None).unwrap(), None);
        assert!(!p.exists());
        std::fs::remove_dir_all(&dir).ok();
    }
}
