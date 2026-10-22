//! The files of a batch export next to the results CSV, built in memory: the results table, the
//! replicate groups, the group statistics, the settings and the quant ion integrals. The caller
//! writes them, converts the German spreadsheet variant and removes the stale ones; no file is
//! touched here.

use std::io;
use std::path::{Path, PathBuf};

use super::{QuantSettings, RunResult, Table, TableSettings};
use super::{group_stats_path, groups_path, ions_path, settings_path};

/// Which file of the export it is.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExportPart {
    Results,
    Groups,
    GroupStats,
    Settings,
    Ions,
}

/// A file to write. `contents` is the error that stops it, e.g. ion integrals of runs that do not
/// match the table.
pub struct ExportFile {
    pub part: ExportPart,
    pub path: PathBuf,
    pub contents: io::Result<Vec<u8>>,
    /// Leading text columns the German spreadsheet variant keeps as text.
    pub text_cols: usize,
}

/// One step of an export: a file to write, or a file of an earlier export that must go because
/// this one has no such part (a script would read its old numbers).
pub enum ExportEntry {
    File(ExportFile),
    Stale(ExportPart, PathBuf),
}

/// What an export writes besides the table itself.
pub struct ExportOptions<'a> {
    pub quant: &'a QuantSettings,
    pub settings: &'a TableSettings,
    /// Settings of the shown table that `TableSettings` does not hold.
    pub extra: &'a [(&'static str, String)],
    /// The runs in the order of the table, when the quant ion integrals are written.
    pub ions: Option<&'a [&'a RunResult]>,
    /// Write the group statistics against `blank` (if the table has groups).
    pub group_stats: bool,
    pub blank: &'a str,
}

fn bytes(write: impl FnOnce(&mut Vec<u8>) -> io::Result<()>) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    write(&mut out)?;
    Ok(out)
}

/// The export of `t` next to the results CSV `csv`, in the order the files are written.
pub fn export_set(t: &Table, csv: &Path, o: &ExportOptions) -> Vec<ExportEntry> {
    let file = |part, path: PathBuf, contents, text_cols| {
        ExportEntry::File(ExportFile {
            part,
            path,
            contents,
            text_cols,
        })
    };
    let mut out = vec![file(
        ExportPart::Results,
        csv.to_path_buf(),
        bytes(|w| t.write_csv(w)),
        1,
    )];
    let groups = t.has_groups();
    out.push(match groups {
        true => file(
            ExportPart::Groups,
            groups_path(csv),
            bytes(|w| t.write_groups_csv(w)),
            2,
        ),
        false => ExportEntry::Stale(ExportPart::Groups, groups_path(csv)),
    });
    out.push(match o.group_stats && groups {
        true => file(
            ExportPart::GroupStats,
            group_stats_path(csv),
            bytes(|w| t.write_group_stats_csv(Some(o.blank), w)),
            1,
        ),
        false => ExportEntry::Stale(ExportPart::GroupStats, group_stats_path(csv)),
    });
    out.push(file(
        ExportPart::Settings,
        settings_path(csv),
        bytes(|w| o.settings.write_csv(o.quant, o.extra, w)),
        2,
    ));
    out.push(match o.ions {
        Some(runs) => file(
            ExportPart::Ions,
            ions_path(csv),
            bytes(|w| t.write_ions_csv(runs, o.settings, w)),
            2,
        ),
        None => ExportEntry::Stale(ExportPart::Ions, ions_path(csv)),
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn table(groups: &[&str]) -> Table {
        Table {
            runs: (0..groups.len()).map(|i| format!("r{i}")).collect(),
            sources: (0..groups.len()).collect(),
            groups: groups.iter().map(|g| g.to_string()).collect(),
            rows: Vec::new(),
        }
    }

    fn parts(set: &[ExportEntry]) -> Vec<(ExportPart, Option<PathBuf>)> {
        set.iter()
            .map(|e| match e {
                ExportEntry::File(f) => (f.part, Some(f.path.clone())),
                ExportEntry::Stale(p, _) => (*p, None),
            })
            .collect()
    }

    fn options<'a>(quant: &'a QuantSettings, settings: &'a TableSettings) -> ExportOptions<'a> {
        ExportOptions {
            quant,
            settings,
            extra: &[],
            ions: None,
            group_stats: true,
            blank: "x",
        }
    }

    #[test]
    fn a_table_with_groups_writes_every_companion_next_to_the_csv() {
        let (quant, settings) = (QuantSettings::default(), TableSettings::default());
        let t = table(&["x", "x", "y"]);
        let set = export_set(&t, Path::new("out.csv"), &options(&quant, &settings));
        let names: Vec<_> = parts(&set);
        assert_eq!(
            names,
            [
                (ExportPart::Results, Some(PathBuf::from("out.csv"))),
                (ExportPart::Groups, Some(PathBuf::from("out.groups.csv"))),
                (
                    ExportPart::GroupStats,
                    Some(PathBuf::from("out.groupstats.csv"))
                ),
                (
                    ExportPart::Settings,
                    Some(PathBuf::from("out.settings.csv"))
                ),
                // No ions asked: an older ion file must go.
                (ExportPart::Ions, None),
            ]
        );
        let ExportEntry::File(results) = &set[0] else {
            panic!("the results file comes first");
        };
        let mut direct = Vec::new();
        t.write_csv(&mut direct).unwrap();
        assert_eq!(results.contents.as_ref().unwrap(), &direct);
        assert_eq!(results.text_cols, 1);
    }

    #[test]
    fn without_groups_or_group_statistics_their_old_files_are_removed() {
        let (quant, settings) = (QuantSettings::default(), TableSettings::default());
        let t = table(&["", ""]);
        let o = ExportOptions {
            group_stats: false,
            ..options(&quant, &settings)
        };
        let set = export_set(&t, Path::new("out.csv"), &o);
        let names: Vec<_> = parts(&set);
        assert_eq!(names[1], (ExportPart::Groups, None));
        assert_eq!(names[2], (ExportPart::GroupStats, None));
        assert!(matches!(&set[1], ExportEntry::Stale(_, p) if p.ends_with("out.groups.csv")));
    }

    #[test]
    fn ion_integrals_need_one_scanned_run_per_column() {
        let (quant, settings) = (QuantSettings::default(), TableSettings::default());
        let t = table(&["x", "x", "y"]);
        let none: Vec<&RunResult> = Vec::new();
        let o = ExportOptions {
            ions: Some(&none),
            ..options(&quant, &settings)
        };
        let set = export_set(&t, Path::new("out.csv"), &o);
        let ExportEntry::File(ions) = &set[4] else {
            panic!("the ion file is written, its contents may fail");
        };
        assert_eq!(ions.part, ExportPart::Ions);
        assert_eq!(
            ions.contents.as_ref().unwrap_err().kind(),
            io::ErrorKind::InvalidInput
        );
    }
}
