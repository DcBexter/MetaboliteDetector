//! "Export MIDs" of the batch table (old `BatchQuantificationResultWidget::exportMIDs`,
//! algorithms.md §7g.5): the `FRAG` feature of a substance, e.g. `(174/8)(261/5)`, names fragments
//! by first ion and number of ions, and the export lists their uncorrected isotopomer shares per run.

use std::fmt::Write;

use regex::Regex;

use super::MAX_FRAGMENT_IONS;
use crate::batch::{Found, RunResult};
use crate::library::LibraryCompound;

/// `(174/8)(261/5)` → `[(174, 8), (261, 5)]`: (first ion, number of ions), sorted by ion; repeated
/// ions come latest first, as in the old `QMultiMap`. Empty when nothing parses. Fragments of more
/// than [`MAX_FRAGMENT_IONS`] ions are dropped; the old program would loop over all of them.
pub fn parse_fragments(text: &str) -> Vec<(u32, u32)> {
    let re = Regex::new(r"\((\d+)/(\d+)\)").expect("valid pattern");
    let mut out: Vec<(u32, u32)> = (re.captures_iter(text))
        // Old `QString::toInt`: up to i32::MAX, so `ion + size` cannot overflow a u32.
        .filter_map(|c| {
            let n = |i: usize| c[i].parse::<i32>().ok().map(|v| v as u32);
            Some((n(1)?, n(2).filter(|&s| s <= MAX_FRAGMENT_IONS)?))
        })
        .collect();
    out.reverse();
    out.sort_by_key(|f| f.0);
    out
}

/// The export file and the messages the old program shows after writing it.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct FragExport {
    pub csv: String,
    /// Substances without `FRAG`.
    pub warnings: Vec<String>,
    /// `FRAG` values that do not parse.
    pub errors: Vec<String>,
}

/// Writes the three tables of the old export for `compounds` (library substance and its index in
/// `Found::compound`, in table order): retention times in minutes (-1 where a run lacks the
/// substance), total ion counts of each fragment, and each ion's share of its fragment's total.
pub fn frag_export(runs: &[RunResult], compounds: &[(usize, &LibraryCompound)]) -> FragExport {
    let names: Vec<String> = (runs.iter())
        .map(|r| {
            let name = r.path.file_name().unwrap_or_default().to_string_lossy();
            name.strip_suffix(".bin").unwrap_or(&name).to_owned()
        })
        .collect();
    let n = compounds.iter().map(|c| c.0 + 1).max().unwrap_or(0);
    let hits = super::series::hits(runs, n);
    let found = |i: usize| &hits[i];
    let integral = |f: Option<&Found>, m: u32| {
        f.and_then(|f| f.ions.iter().find(|p| p.mass == m))
            .map_or(0.0, |p| p.integral)
    };
    let frag = |c: &LibraryCompound| c.features.get("FRAG").filter(|v| !v.is_empty()).cloned();
    let mut out = FragExport::default();
    let csv = &mut out.csv;

    csv.push_str("Retention Times (Min)\nMetabolite\t");
    names.iter().for_each(|n| _ = write!(csv, "{n}\t"));
    csv.push('\n');
    for &(i, c) in compounds.iter().filter(|(_, c)| frag(c).is_some()) {
        csv.push_str(&c.name);
        for f in found(i) {
            _ = write!(csv, "\t{}", f.map_or(-1.0, |f| f.rt_s / 60.0));
        }
        csv.push('\n');
    }

    csv.push_str("\nTotal ion counts\nFragment\t");
    names.iter().for_each(|n| _ = write!(csv, "{n}\t"));
    csv.push('\n');
    // Per substance: per fragment, the total of each run.
    let mut totals: Vec<Vec<Vec<f64>>> = Vec::new();
    for &(i, c) in compounds {
        let Some(text) = frag(c) else {
            out.warnings.push(format!(
                "Warning: No fragment (FRAG) defined for compound {}",
                c.name
            ));
            totals.push(Vec::new());
            continue;
        };
        let parsed = parse_fragments(&text);
        if parsed.is_empty() {
            out.errors.push(format!(
                "Error: Can not parse \"{text}\". SYNTAX(e.g): (174/8)(261/5)"
            ));
        }
        let hits = found(i);
        let mut mine = Vec::new();
        for (ion, size) in parsed {
            _ = write!(csv, "{}_{ion}", c.name);
            let sums: Vec<f64> = (hits.iter())
                .map(|f| (ion..ion + size).map(|m| integral(*f, m)).sum())
                .collect();
            sums.iter().for_each(|s| _ = write!(csv, "\t{s}"));
            csv.push('\n');
            mine.push(sums);
        }
        totals.push(mine);
    }

    csv.push_str("\nMass Isotopomer Abundances\nIon");
    names.iter().for_each(|n| _ = write!(csv, "\t{n}"));
    csv.push('\n');
    for (&(i, c), sums) in compounds.iter().zip(&totals) {
        let Some(text) = frag(c) else { continue };
        let hits = found(i);
        for ((ion, size), sums) in parse_fragments(&text).into_iter().zip(sums) {
            for m in ion..ion + size {
                _ = write!(csv, "{}_{m} (M{})", c.name, m - ion);
                for (f, total) in hits.iter().zip(sums) {
                    _ = write!(csv, "\t{}", integral(*f, m) / total);
                }
                csv.push('\n');
            }
            csv.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::IonPeak;
    use std::path::PathBuf;

    #[test]
    fn fragments_parse_like_the_old_program() {
        assert_eq!(parse_fragments("(174/8)(261/5)"), [(174, 8), (261, 5)]);
        assert_eq!(parse_fragments("x (261/5) (174/8)"), [(174, 8), (261, 5)]);
        // Repeated ions: the later entry first.
        assert_eq!(parse_fragments("(5/1)(5/2)"), [(5, 2), (5, 1)]);
        assert!(parse_fragments("174/8").is_empty());
        assert!(parse_fragments("(99999999999/1)").is_empty());
        assert!(parse_fragments("(2147483648/1)").is_empty());
        assert!(parse_fragments("(1/1001)").is_empty());
        assert_eq!(
            parse_fragments("(2147483647/1000)"),
            [(i32::MAX as u32, 1000)]
        );
    }

    #[test]
    fn export_lists_totals_and_shares() {
        let ion = |mass, integral| IonPeak {
            mass,
            integral,
            quality: 0.0,
            signal_to_noise: 0.0,
        };
        let run = |name: &str, found: Vec<Found>| RunResult {
            path: PathBuf::from(name),
            found,
            unknown: Vec::new(),
            calibration: None,
        };
        let hit = |ions| Found {
            compound: 0,
            score: 1.0,
            rt_s: 90.0,
            ions,
            rescanned: false,
        };
        let mut ala = LibraryCompound::new("Ala");
        ala.features.insert("FRAG".into(), "(10/2)".into());
        let gly = LibraryCompound::new("Gly");
        let runs = [
            run("a.bin", vec![hit(vec![ion(10, 3.0), ion(11, 1.0)])]),
            run("b.cdf", vec![]),
        ];
        let e = frag_export(&runs, &[(0, &ala), (1, &gly)]);
        assert_eq!(
            e.csv,
            "Retention Times (Min)\nMetabolite\ta\tb.cdf\t\nAla\t1.5\t-1\n\
             \nTotal ion counts\nFragment\ta\tb.cdf\t\nAla_10\t4\t0\n\
             \nMass Isotopomer Abundances\nIon\ta\tb.cdf\nAla_10 (M0)\t0.75\tNaN\nAla_11 (M1)\t0.25\tNaN\n\n"
        );
        assert_eq!(
            e.warnings,
            ["Warning: No fragment (FRAG) defined for compound Gly"]
        );
        assert!(e.errors.is_empty());
    }
}
