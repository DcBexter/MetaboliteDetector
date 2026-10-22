//! What lies behind a cell of the result table: the integral of each quant ion, the other library
//! hits of the component, and the files that document an export (all ions, settings).

use std::collections::HashMap;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use super::{
    BOM, Found, IonPeak, QuantSettings, RunResult, Table, TableSettings, path_key, run_key,
    text_cell,
};
use crate::library::LibraryCompound;

/// The largest peak of each of `ions` in `f` that has at least `min_sn`; a duplicate peak (§4.5)
/// counts once. `None` where the ion has no such peak.
pub(super) fn quant_peaks<'a>(f: &'a Found, ions: &[u32], min_sn: f32) -> Vec<Option<&'a IonPeak>> {
    (ions.iter())
        .map(|m| (f.ions.iter()).find(|p| p.mass == *m && p.signal_to_noise >= min_sn))
        .collect()
}

/// The integral of each of `ions` in `f`, in `ions` order; `None` where the ion was not found.
pub fn ion_integrals(f: &Found, ions: &[u32], min_sn: f32) -> Vec<Option<f64>> {
    quant_peaks(f, ions, min_sn)
        .into_iter()
        .map(|p| p.map(|p| p.integral))
        .collect()
}

/// The integral of each of `ions` in `f` as far as the cell of `f` counts under `s`, else all `None`.
/// Independent of normalisation: whether a substance was found does not depend on its shown value.
pub fn cell_ions(f: &Found, ions: &[u32], s: &TableSettings) -> Vec<Option<f64>> {
    let need = s.min_ions.clamp(1, ions.len().max(1));
    match super::cell(f, ions, need, s) {
        Some(_) => ion_integrals(f, ions, s.min_sn),
        None => vec![None; ions.len()],
    }
}

/// The library hits of the component `f` came from, the assigned one first (the old tooltip shows
/// it in red), each with its overall score. Indices past `lib` address `extra`.
pub fn hit_list(
    f: &Found,
    lib: &[LibraryCompound],
    extra: &[LibraryCompound],
) -> Vec<(String, f64)> {
    let name = |i: usize| {
        (lib.get(i).or_else(|| extra.get(i.checked_sub(lib.len())?)))
            .map_or_else(|| format!("#{i}"), |c| c.name.clone())
    };
    std::iter::once((f.compound, f.score))
        .chain(f.alternatives.iter().copied())
        .map(|(i, score)| (name(i), score))
        .collect()
}

/// Where the integrals of all quant ions go next to the results CSV `csv`:
/// `results.csv` → `results.ions.csv`.
pub fn ions_path(csv: &Path) -> PathBuf {
    csv.with_extension("ions.csv")
}

/// Where the settings of an export go next to the results CSV `csv`:
/// `results.csv` → `results.settings.csv`.
pub fn settings_path(csv: &Path) -> PathBuf {
    csv.with_extension("settings.csv")
}

impl Table {
    /// The integral of every quant ion, one line per substance and ion:
    /// `Substance,Quant_ion,<one column per run>`, same encoding as `write_csv`. `runs` are the
    /// scanned runs of the table's columns, in order; a value is empty where the ion or the whole
    /// cell was not found. The integrals are the raw ones, whatever normalisation the table shows.
    pub fn write_ions_csv(
        &self,
        runs: &[&RunResult],
        s: &TableSettings,
        mut w: impl Write,
    ) -> io::Result<()> {
        if runs.len() != self.runs.len() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "one scanned run per column needed",
            ));
        }
        w.write_all(BOM.as_bytes())?;
        let head = ["Substance", "Quant_ion"].map(String::from);
        let head: Vec<String> = (head.into_iter())
            .chain(self.runs.iter().map(|r| text_cell(r)))
            .collect();
        writeln!(w, "{}", head.join(","))?;
        // The hit of each run, found once per run instead of once per cell.
        let by_compound: Vec<HashMap<usize, &Found>> = runs
            .iter()
            .map(|r| r.found.iter().map(|f| (f.compound, f)).collect())
            .collect();
        for r in &self.rows {
            let per_run: Vec<Vec<Option<f64>>> = by_compound
                .iter()
                .map(|hits| match hits.get(&r.compound) {
                    Some(f) => cell_ions(f, &r.ions, s),
                    None => vec![None; r.ions.len()],
                })
                .collect();
            for (k, ion) in r.ions.iter().enumerate() {
                let cells = per_run
                    .iter()
                    .map(|v| v[k].map_or(String::new(), |v| v.to_string()));
                let line = [text_cell(&r.name), ion.to_string()]
                    .into_iter()
                    .chain(cells);
                writeln!(w, "{}", line.collect::<Vec<_>>().join(","))?;
            }
        }
        w.flush()
    }
}

impl TableSettings {
    /// Where the replicate groups come from, for the settings file: runs set by hand, the sample
    /// table (`table` describes it), the file name pattern.
    /// Only sources that group some of `runs` are named.
    pub fn group_source(&self, table: Option<&str>, runs: &[PathBuf]) -> String {
        let by_hand = |p: &PathBuf| self.manual_groups.contains_key(&path_key(p));
        let by_table = |p: &PathBuf| {
            self.sample_groups
                .contains_key(&run_key(&p.to_string_lossy()))
        };
        let n = runs.iter().filter(|p| by_hand(p)).count();
        let rest = || runs.iter().filter(|p| !by_hand(p));
        let hand = (n > 0).then(|| format!("by hand ({n} runs)"));
        let table = table
            .filter(|_| rest().any(by_table))
            .map(|t| format!("sample table {t}"));
        let pattern = (!self.replicates.is_empty() && rest().any(|p| !by_table(p)))
            .then(|| "file name pattern".to_owned());
        let parts: Vec<String> = [hand, table, pattern].into_iter().flatten().collect();
        if parts.is_empty() {
            "none".into()
        } else {
            parts.join(", ")
        }
    }

    /// The settings the table was made with as `Setting,Value`, so a result file can be read
    /// later; `extra` are settings of the caller (e.g. the normalisation) written after them.
    pub fn write_csv(
        &self,
        quant: &QuantSettings,
        extra: &[(&str, String)],
        mut w: impl Write,
    ) -> io::Result<()> {
        w.write_all(BOM.as_bytes())?;
        writeln!(w, "Setting,Value")?;
        let q = quant;
        let excluded: Vec<String> = (q.excluded.iter())
            .map(|(a, b)| format!("{a}-{b}"))
            .collect();
        let rows = [
            ("QUANT_ION_NUMBER", q.ion_number.to_string()),
            ("QUANT_MIN_DIST_IONS", q.min_dist.to_string()),
            ("QUANT_QUALITY_CUTOFF", q.quality_cutoff.to_string()),
            ("QUANT_EXCL_IONS", excluded.join(";")),
            ("REPLICATE_PATTERN", self.replicates.clone()),
            ("MIN_REPRODUCIBILITY", self.min_reproducibility.to_string()),
            ("MIN_QUANT_IONS", self.min_ions.to_string()),
            ("MIN_SN", self.min_sn.to_string()),
            ("MAX_DISCREPANCY", self.max_discrepancy.to_string()),
        ];
        let extra = extra.iter().map(|(k, v)| (*k, v.clone()));
        for (k, v) in rows.into_iter().chain(extra) {
            writeln!(w, "{},{}", text_cell(k), text_cell(&v))?;
        }
        w.flush()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::Row;

    fn found(compound: usize, ions: &[(u32, f64, f32)]) -> Found {
        Found {
            compound,
            score: 0.9,
            rt_s: 60.0,
            ions: (ions.iter())
                .map(|&(mass, integral, signal_to_noise)| IonPeak {
                    mass,
                    integral,
                    quality: 0.5,
                    signal_to_noise,
                })
                .collect(),
            rescanned: false,
            rescued: Vec::new(),
            alternatives: vec![(1, 0.8), (7, 0.6)],
        }
    }

    #[test]
    fn ion_integrals_skip_missing_and_noisy_ions() {
        let f = found(0, &[(73, 500.0, 10.0), (147, 20.0, 1.0), (73, 5.0, 10.0)]);
        // The first (largest) peak of 73 counts, 147 is below the S/N limit, 200 is absent.
        assert_eq!(
            ion_integrals(&f, &[73, 147, 200], 3.0),
            [Some(500.0), None, None]
        );
        assert_eq!(ion_integrals(&f, &[147], 0.0), [Some(20.0)]);
    }

    #[test]
    fn hit_list_has_the_assigned_hit_first_and_names_the_pooled_unknowns() {
        let c = LibraryCompound::new;
        let (lib, extra) = ([c("A"), c("B")], [c("Unknown 1")]);
        let mut f = found(0, &[]);
        f.alternatives = vec![(1, 0.8), (2, 0.6), (9, 0.1)];
        assert_eq!(
            hit_list(&f, &lib, &extra),
            [
                ("A".to_string(), 0.9),
                ("B".to_string(), 0.8),
                ("Unknown 1".to_string(), 0.6),
                ("#9".to_string(), 0.1)
            ]
        );
    }

    #[test]
    fn ions_csv_has_a_line_per_ion_and_leaves_unquantified_cells_empty() {
        let two = TableSettings {
            min_ions: 2,
            ..Default::default()
        };
        let run = |f: Found| RunResult {
            input: 0,
            path: "a.cdf".into(),
            found: vec![f],
            ..Default::default()
        };
        let (a, b) = (
            run(found(0, &[(73, 500.0, 10.0), (147, 100.0, 10.0)])),
            run(found(0, &[(73, 50.0, 10.0)])),
        );
        let t = Table {
            runs: vec!["a.cdf".into(), "=b.cdf".into()],
            groups: vec![String::new(); 2],
            rows: vec![Row {
                compound: 0,
                name: "Glucose".into(),
                rt_s: 60.0,
                ions: vec![73, 147],
                // Run b is below the table's minimum of two ions: the whole cell is not quantified.
                values: vec![Some(600.0), None],
                rts_s: vec![Some(60.0), None],
            }],
            ..Default::default()
        };
        let mut out = Vec::new();
        t.write_ions_csv(&[&a, &b], &two, &mut out).unwrap();
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "\u{feff}Substance,Quant_ion,a.cdf,'=b.cdf\nGlucose,73,500,\nGlucose,147,100,\n"
        );
        assert!(t.write_ions_csv(&[&a], &two, Vec::new()).is_err());
        // A normalisation that blanks the shown value (e.g. z-scores) does not blank the integrals.
        let mut z = t.clone();
        z.rows[0].values = vec![None, None];
        let mut out = Vec::new();
        z.write_ions_csv(&[&a, &b], &two, &mut out).unwrap();
        assert!(
            String::from_utf8(out)
                .unwrap()
                .contains("\nGlucose,73,500,\n")
        );
        assert_eq!(
            ions_path(Path::new("d/results.csv")),
            Path::new("d/results.ions.csv")
        );
    }

    #[test]
    fn the_group_source_names_hand_table_and_pattern() {
        let runs: Vec<PathBuf> = ["/x/a.cdf", "/x/b.cdf", "/y/b.cdf"].map(Into::into).into();
        let mut s = TableSettings::default();
        assert_eq!(s.group_source(None, &runs), "none");
        s.replicates = "^(.+)_".into();
        s.manual_groups = crate::batch::SampleGroups::from([(path_key(&runs[0]), "blank".into())]);
        s.sample_groups = crate::batch::SampleGroups::from([("b".into(), "t1".into())]);
        let table = Some("samples.txt");
        assert_eq!(
            s.group_source(table, &runs),
            "by hand (1 runs), sample table samples.txt"
        );
        // A removed run no longer counts, and a source that groups nothing is not named.
        assert_eq!(s.group_source(table, &runs[..1]), "by hand (1 runs)");
        assert_eq!(
            s.group_source(table, &runs[1..]),
            "sample table samples.txt"
        );
        s.sample_groups.clear();
        assert_eq!(s.group_source(table, &runs[1..]), "file name pattern");
    }

    #[test]
    fn settings_csv_lists_the_filters_and_the_callers_extras() {
        let mut s = TableSettings::default();
        let q = QuantSettings {
            excluded: vec![(73.0, 75.0), (147.0, 147.0)],
            ..Default::default()
        };
        s.replicates = "(\\w+)_\\d".into();
        let mut out = Vec::new();
        s.write_csv(&q, &[("NORMALISATION", "total signal".into())], &mut out)
            .unwrap();
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with("\u{feff}Setting,Value\nQUANT_ION_NUMBER,3\n"));
        assert!(text.contains("\nQUANT_EXCL_IONS,73-75;147-147\n"));
        assert!(text.contains("\nREPLICATE_PATTERN,(\\w+)_\\d\n"));
        assert!(text.ends_with("\nNORMALISATION,total signal\n"));
    }
}
