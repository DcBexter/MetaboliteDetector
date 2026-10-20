//! Sample table: replicate groups that the file names do not tell, e.g. the old `samples.txt` with the
//! columns `CDF_FILE` and `TIME_POINT`. CSV, TSV or semicolon-separated, first line = column names.

use std::collections::BTreeMap;
use std::path::Path;

/// Run key (lower-case file name without extension) → replicate group.
pub type SampleGroups = BTreeMap<String, String>;

/// Key of a run or table entry: path and extension do not matter, nor the case (Windows).
pub fn run_key(name: &str) -> String {
    let file = name.rsplit(['/', '\\']).next().unwrap_or(name);
    let stem = Path::new(file).file_stem().unwrap_or_default();
    stem.to_string_lossy().to_lowercase()
}

#[derive(Debug, Clone, PartialEq)]
pub struct SampleTable {
    pub columns: Vec<String>,
    rows: Vec<Vec<String>>,
}

impl SampleTable {
    pub fn parse(text: &str) -> Result<Self, String> {
        let text = text.trim_start_matches('\u{feff}');
        let mut lines = text.lines().filter(|l| !l.trim().is_empty());
        let head = lines.next().ok_or("the table is empty")?;
        // shortcut: no quoted delimiters inside cells, upgrade when a real table needs them.
        let sep = ['\t', ';', ','].into_iter().find(|c| head.contains(*c));
        let cells = |l: &str| -> Vec<String> {
            l.split(sep.unwrap_or('\t'))
                .map(|c| c.trim().trim_matches('"').trim().to_owned())
                .collect()
        };
        Ok(Self {
            columns: cells(head),
            rows: lines.map(cells).collect(),
        })
    }

    /// Index of the column called `name` (case does not matter).
    pub fn column(&self, name: &str) -> Result<usize, String> {
        (self.columns.iter())
            .position(|c| c.eq_ignore_ascii_case(name))
            .ok_or_else(|| {
                format!(
                    "no column {name}; the table has {}",
                    self.columns.join(", ")
                )
            })
    }

    /// Group of each run named in column `file`, taken from column `group`. Rows without a file or a
    /// group are skipped; a run named twice with different groups is an error.
    pub fn groups(&self, file: usize, group: usize) -> Result<SampleGroups, String> {
        let mut out = SampleGroups::new();
        for (i, row) in self.rows.iter().enumerate() {
            let (Some(f), Some(g)) = (row.get(file), row.get(group)) else {
                return Err(format!("row {}: too few columns", i + 2));
            };
            if f.is_empty() || g.is_empty() {
                continue;
            }
            if let Some(old) = out.insert(run_key(f), g.clone())
                && old != *g
            {
                return Err(format!("{f} is in groups {old} and {g}"));
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reads_groups_from_a_tab_separated_table() {
        let text = "\u{feff}CDF_FILE\tTIME_POINT\tNOTE\n7235eg01.CDF\t0\tx\n\"C:\\d\\7235eg02.cdf\"\t0\n7235eg03.cdf\t3\n\n\
                    blank.cdf\t\n";
        let t = SampleTable::parse(text).unwrap();
        assert_eq!(t.column("time_point"), Ok(1));
        let g = t.groups(0, 1).unwrap();
        assert_eq!(g.len(), 3);
        assert_eq!(g["7235eg02"], "0");
        assert_eq!(g[&run_key("/runs/7235EG03.cdf")], "3");
    }

    #[test]
    fn reads_csv_and_semicolon_tables() {
        for sep in [",", ";"] {
            let t = SampleTable::parse(&format!("File{sep}Group\na.cdf{sep}x\n")).unwrap();
            assert_eq!(t.groups(0, 1).unwrap()["a"], "x");
        }
    }

    #[test]
    fn rejects_bad_tables() {
        assert!(SampleTable::parse("  \n").is_err());
        let t = SampleTable::parse("File,Group\na.cdf,x\nb.cdf\na.cdf,y\n").unwrap();
        assert!(t.column("nope").unwrap_err().contains("File, Group"));
        assert!(t.groups(0, 1).unwrap_err().starts_with("row 3"));
        let t = SampleTable::parse("File,Group\na.cdf,x\nA.CDF,y\n").unwrap();
        assert!(t.groups(0, 1).unwrap_err().contains("groups x and y"));
    }
}
