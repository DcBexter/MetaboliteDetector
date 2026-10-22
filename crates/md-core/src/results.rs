//! Result table of one run: one row per component with its library match, the old `DataTableModel`
//! (columns Compound Name, RT, RI, Peak Area, Overall, Spec. and RI Score, Quant. Ions) and
//! `CompoundMapper::findCompound`, see docs/reverse-engineering/algorithms.md §7f.

use std::io::{self, Write};

use crate::csv::quote;
use crate::deconv::Component;
use crate::library::LibraryCompound;
use crate::search::Identification;

/// One component of a run with the library compound assigned to it.
#[derive(Debug, Clone, PartialEq)]
pub struct ResultRow {
    /// Index into the components of the run.
    pub component: usize,
    /// Library index of the assigned compound; `None` for a component without match.
    pub compound: Option<usize>,
    /// The library name, else `No match: ` and the component's own name.
    pub name: String,
    pub rt_s: f64,
    /// `None` while the run is not RI-calibrated.
    pub ri: Option<f64>,
    /// Summed integral of the quant ions; `None` without match (old value -1).
    pub area: Option<f64>,
    pub overall: Option<f64>,
    pub spec_score: Option<f64>,
    pub ri_score: Option<f64>,
    /// Quant ions of the library compound, empty when it has none.
    pub quant_ions: Vec<f32>,
}

/// Rows of all components in component order. `peak_masses[p]` is the ion mass of
/// peak `p` of the peak list the components were deconvoluted from; `ris` as in `map_compounds_to_library`.
/// An assigned hit below `cutoff` (`LS_CUTOFF`) is shown as no match, like the old `DataTableModel`.
pub fn result_rows(
    comps: &[Component],
    scan_times_s: &[f64],
    ris: Option<&[f64]>,
    peak_masses: &[u32],
    id: &Identification,
    lib: &[LibraryCompound],
    cutoff: f64,
) -> Vec<ResultRow> {
    comps
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let hit = (id.assigned.get(i).copied().flatten()).filter(|h| h.overall >= cutoff);
            let l = hit.map(|h| &lib[h.compound]);
            let ions = l.map_or(&[][..], |l| match l.quant_ions.is_empty() {
                true => &l.auto_quant_ions,
                false => &l.quant_ions,
            });
            ResultRow {
                component: i,
                compound: hit.map(|h| h.compound),
                // shortcut: the old component name is not decoded (a `Compound` has none in `.cmp`), so
                // the unmatched row shows its number; replace when .cmp reading (#13) knows the name.
                name: l.map_or_else(|| format!("No match: {}", i + 1), |l| l.name.clone()),
                rt_s: scan_times_s[c.apex_scan as usize],
                ri: ris.and_then(|r| r.get(i).copied()),
                area: hit.map(|_| area(c, ions, peak_masses)),
                overall: hit.map(|h| h.overall),
                spec_score: hit.map(|h| h.spec_sim),
                ri_score: hit.and_then(|h| h.ri_sim),
                quant_ions: l.map_or_else(Vec::new, |_| ions.to_vec()),
            }
        })
        .collect()
}

/// The one apex rule: component indices by apex scan, stable, so equal apexes keep the deconvolution
/// order. The list sorts by it and ↓/↑ with the list hidden steps along it; the run keeps its order.
pub fn apex_order(comps: &[Component]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..comps.len()).collect();
    order.sort_by_key(|&k| comps[k].apex_scan);
    order
}

/// The rows of the component list in [`apex_order`] of their components (rows come in component order).
pub fn sort_by_apex(rows: &mut [ResultRow], comps: &[Component]) {
    rows.sort_by_key(|r| comps[r.component].apex_scan);
}

/// Old `getIntegral`: the sum over the quant ions of the largest peak of that ion (a duplicate peak
/// counts once, §4.5). Without quant ions the old code asks the component for three ions; here these are
/// the three largest peaks (Inferred).
pub fn area(c: &Component, ions: &[f32], peak_masses: &[u32]) -> f64 {
    let mut by_mass: Vec<(u32, f64)> = c
        .peaks
        .iter()
        .map(|f| (peak_masses[f.peak], f.integral))
        .collect();
    by_mass.sort_by(|a, b| a.0.cmp(&b.0).then(b.1.total_cmp(&a.1)));
    by_mass.dedup_by_key(|p| p.0);
    if ions.is_empty() {
        by_mass.sort_by(|a, b| b.1.total_cmp(&a.1));
        return by_mass.iter().take(3).map(|p| p.1).sum();
    }
    let mut masses: Vec<u32> = ions.iter().map(|m| m.round() as u32).collect();
    masses.sort_unstable();
    masses.dedup();
    masses
        .iter()
        .filter_map(|m| by_mass.binary_search_by_key(m, |p| p.0).ok())
        .map(|k| by_mass[k].1)
        .sum()
}

/// RFC 4180 CSV of the table, one line per row; `.` decimals, `,` separators, retention time in minutes.
/// Missing values are empty cells (the old export wrote -1), which R reads as `NA`. The old program wrote
/// the same columns separated by tabs. Quant ions are space-separated.
pub fn write_csv(mut w: impl Write, rows: &[ResultRow]) -> io::Result<()> {
    let head = "Compound Name,RT (min),RI,Peak Area,Overall Score,Spec. Score,RI Score,Quant. Ions";
    writeln!(w, "{head}")?;
    let num = |v: Option<f64>| v.map(|v| v.to_string()).unwrap_or_default();
    for r in rows {
        let ions: Vec<String> = r.quant_ions.iter().map(f32::to_string).collect();
        writeln!(
            w,
            "{},{},{},{},{},{},{},{}",
            quote(&r.name),
            r.rt_s / 60.0,
            num(r.ri),
            num(r.area),
            num(r.overall),
            num(r.spec_score),
            num(r.ri_score),
            ions.join(" ")
        )?;
    }
    w.flush()
}

/// Old `CompoundMapper::findCompound`: the components that have a hit of at least `cutoff`
/// (`LS_CUTOFF`) with the same name as library compound `compound`, best score first. Returns
/// (component index, overall score).
pub fn find_substance(
    id: &Identification,
    lib: &[LibraryCompound],
    compound: usize,
    cutoff: f64,
) -> Vec<(usize, f64)> {
    let name = &lib[compound].name;
    let mut found: Vec<(usize, f64)> = id
        .hits
        .iter()
        .enumerate()
        .flat_map(|(i, hits)| {
            hits.iter()
                .filter(|h| h.overall >= cutoff && lib[h.compound].name == *name)
                .map(move |h| (i, h.overall))
        })
        .collect();
    found.sort_by(|a, b| b.1.total_cmp(&a.1));
    found
}

/// A column of the component list, in the order of the screen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Col {
    Name,
    Rt,
    Ri,
    Area,
    Score,
    Spectrum,
    RiScore,
    Sn,
    Ions,
    QuantIons,
}

/// Without a library there are no scores and no quant ions.
pub const PLAIN: [Col; 6] = [Col::Name, Col::Rt, Col::Ri, Col::Area, Col::Sn, Col::Ions];
pub const WITH_LIBRARY: [Col; 10] = [
    Col::Name,
    Col::Rt,
    Col::Ri,
    Col::Area,
    Col::Score,
    Col::Spectrum,
    Col::RiScore,
    Col::Sn,
    Col::Ions,
    Col::QuantIons,
];

/// A row of the component list: the result row and what the list shows beyond the CSV.
#[derive(Debug, Clone, PartialEq)]
pub struct ListRow {
    pub res: ResultRow,
    /// Also without a match: then the three largest peaks (see [`area`]).
    pub area: f64,
    /// Mean S/N of the ion peaks; `None` for a component without peaks (saved by the old program).
    pub sn: Option<f64>,
    pub ions: usize,
}

/// Click order of a column header: ascending, descending, off.
pub fn next_sort(sort: Option<(Col, bool)>, col: Col) -> Option<(Col, bool)> {
    match sort {
        Some((s, false)) if s == col => Some((col, true)),
        Some((s, true)) if s == col => None,
        _ => Some((col, false)),
    }
}

/// The rows the list shows, in their order: those whose shown name (`shown`, lower-cased) contains
/// `filter` (case-insensitive; empty keeps all), sorted by `sort` when given. The shown name is
/// the one on screen, so the name sort and the filter follow the language.
pub fn order(
    rows: &[ListRow],
    filter: &str,
    shown: impl Fn(&ResultRow) -> String,
    sort: Option<(Col, bool)>,
) -> Vec<usize> {
    let want = filter.to_lowercase();
    let names: Vec<String> = rows.iter().map(|r| shown(&r.res).to_lowercase()).collect();
    let mut order: Vec<usize> = (0..rows.len())
        .filter(|i| want.is_empty() || names[*i].contains(&want))
        .collect();
    if let Some((col, desc)) = sort {
        order.sort_by(|a, b| {
            let o = compare(&rows[*a], &rows[*b], col, &names[*a], &names[*b]);
            if desc { o.reverse() } else { o }
        });
    }
    order
}

/// Missing values sort first. `na`/`nb` are the lower-cased shown names of `a` and `b`.
fn compare(a: &ListRow, b: &ListRow, col: Col, na: &str, nb: &str) -> std::cmp::Ordering {
    use crate::view::missing_first as num;
    let (x, y) = (&a.res, &b.res);
    match col {
        Col::Name => na.cmp(nb),
        Col::Rt => x.rt_s.total_cmp(&y.rt_s),
        Col::Ri => num(x.ri).total_cmp(&num(y.ri)),
        Col::Area => a.area.total_cmp(&b.area),
        Col::Score => num(x.overall).total_cmp(&num(y.overall)),
        Col::Spectrum => num(x.spec_score).total_cmp(&num(y.spec_score)),
        Col::RiScore => num(x.ri_score).total_cmp(&num(y.ri_score)),
        Col::Sn => num(a.sn).total_cmp(&num(b.sn)),
        Col::Ions => a.ions.cmp(&b.ions),
        Col::QuantIons => x.quant_ions.len().cmp(&y.quant_ions.len()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::deconv::FittedPeak;

    fn row(name: &str, rt_s: f64, ri: Option<f64>, area: f64) -> ListRow {
        ListRow {
            res: ResultRow {
                component: 0,
                compound: Some(0),
                name: name.into(),
                rt_s,
                ri,
                area: Some(area),
                overall: Some(0.5),
                spec_score: None,
                ri_score: None,
                quant_ions: Vec::new(),
            },
            area,
            sn: None,
            ions: 1,
        }
    }

    #[test]
    fn a_header_click_cycles_ascending_descending_off() {
        let up = next_sort(None, Col::Rt);
        assert_eq!(up, Some((Col::Rt, false)));
        let down = next_sort(up, Col::Rt);
        assert_eq!(down, Some((Col::Rt, true)));
        assert_eq!(next_sort(down, Col::Rt), None);
        assert_eq!(next_sort(up, Col::Ri), Some((Col::Ri, false)));
    }

    #[test]
    fn the_list_filters_by_shown_name_and_sorts_missing_values_first() {
        let rows = vec![
            row("Glycine", 120.0, None, 3.0),
            row("Alanine", 60.0, Some(1500.0), 15.0),
            row("No match: 3", 90.0, Some(1400.0), 1.0),
        ];
        let shown = |r: &ResultRow| r.name.clone();
        assert_eq!(order(&rows, "", shown, None), vec![0, 1, 2]);
        assert_eq!(order(&rows, "ALA", shown, None), vec![1]);
        assert_eq!(
            order(&rows, "", shown, Some((Col::Name, false))),
            vec![1, 0, 2]
        );
        assert_eq!(
            order(&rows, "", shown, Some((Col::Ri, false))),
            vec![0, 2, 1]
        );
        assert_eq!(
            order(&rows, "", shown, Some((Col::Area, true))),
            vec![1, 0, 2]
        );
        // The translated name is what the filter and the name sort see.
        let german = |r: &ResultRow| r.name.replace("No match", "Kein Treffer");
        assert_eq!(order(&rows, "treffer", german, None), vec![2]);
        assert_eq!(
            order(&rows, "", german, Some((Col::Name, false))),
            vec![1, 0, 2]
        );
    }
    use crate::search::LibraryHit;

    fn comp(apex: u32, peaks: &[(usize, f64)]) -> Component {
        Component {
            first_scan: 0,
            apex_scan: apex,
            last_scan: 3,
            model: Vec::new(),
            spectrum: Vec::new(),
            peaks: peaks
                .iter()
                .map(|&(peak, integral)| FittedPeak {
                    peak,
                    model_share: 1.0,
                    integral,
                })
                .collect(),
            exact_spectrum: Vec::new(),
        }
    }

    fn hit(compound: usize, overall: f64) -> LibraryHit {
        LibraryHit {
            compound,
            spec_sim: overall - 0.1,
            ri_sim: None,
            overall,
        }
    }

    fn setup() -> (Vec<Component>, Identification, Vec<LibraryCompound>) {
        let mut a = LibraryCompound::new("Alanine, 2TMS");
        a.quant_ions = vec![73.0, 116.0];
        let b = LibraryCompound::new("Glycine");
        // Peak 0 and 1 are ion 73 (the duplicate counts once), peak 2 is ion 116, peak 3 is ion 147.
        let comps = vec![
            comp(1, &[(0, 10.0), (1, 4.0), (2, 5.0), (3, 100.0)]),
            comp(2, &[(0, 1.0)]),
        ];
        let id = Identification {
            hits: vec![vec![hit(0, 0.9), hit(1, 0.85)], vec![hit(1, 0.82)]],
            assigned: vec![Some(hit(0, 0.9)), None],
        };
        (comps, id, vec![a, b])
    }

    #[test]
    fn rows_follow_the_old_table() {
        let (comps, id, lib) = setup();
        let rows = result_rows(
            &comps,
            &[0.0, 60.0, 120.0],
            None,
            &[73, 73, 116, 147],
            &id,
            &lib,
            0.8,
        );
        assert_eq!(rows[0].name, "Alanine, 2TMS");
        assert_eq!(rows[0].rt_s, 60.0);
        assert_eq!(rows[0].area, Some(15.0));
        assert_eq!(rows[0].quant_ions, vec![73.0, 116.0]);
        assert_eq!(rows[1].name, "No match: 2");
        assert_eq!((rows[1].area, rows[1].overall), (None, None));
        // Below the cutoff the assigned hit counts as no match.
        let rows = result_rows(&comps, &[0.0, 60.0, 120.0], None, &[73; 4], &id, &lib, 0.95);
        assert_eq!(
            (rows[0].name.as_str(), rows[0].overall),
            ("No match: 1", None)
        );
    }

    #[test]
    fn rows_sort_by_apex_and_keep_deconvolution_order_on_ties() {
        let comps = vec![comp(2, &[]), comp(1, &[]), comp(2, &[]), comp(0, &[])];
        let id = Identification {
            hits: vec![Vec::new(); 4],
            assigned: vec![None; 4],
        };
        let mut rows = result_rows(&comps, &[0.0, 60.0, 120.0], None, &[], &id, &[], 0.8);
        sort_by_apex(&mut rows, &comps);
        let order: Vec<usize> = rows.iter().map(|r| r.component).collect();
        assert_eq!(order, vec![3, 1, 0, 2]);
        assert_eq!(apex_order(&comps), vec![3, 1, 0, 2]);
    }

    #[test]
    fn csv_has_empty_cells_for_missing_values() {
        let (comps, id, lib) = setup();
        let rows = result_rows(
            &comps,
            &[0.0, 60.0, 120.0],
            Some(&[1500.0, 1600.0]),
            &[73, 73, 116, 147],
            &id,
            &lib,
            0.0,
        );
        let mut out = Vec::new();
        write_csv(&mut out, &rows).unwrap();
        let text = String::from_utf8(out).unwrap();
        let mut lines = text.lines().skip(1);
        assert_eq!(
            lines.next(),
            Some("\"Alanine, 2TMS\",1,1500,15,0.9,0.8,,73 116")
        );
        assert_eq!(lines.next(), Some("No match: 2,2,1600,,,,,"));
    }

    #[test]
    fn find_lists_every_component_with_that_name_best_first() {
        let (_, id, lib) = setup();
        assert_eq!(
            find_substance(&id, &lib, 1, 0.8),
            vec![(0, 0.85), (1, 0.82)]
        );
        assert_eq!(find_substance(&id, &lib, 1, 0.83), vec![(0, 0.85)]);
        assert_eq!(find_substance(&id, &lib, 0, 0.8), vec![(0, 0.9)]);
    }

    #[test]
    fn without_quant_ions_the_three_largest_peaks_count() {
        let c = comp(0, &[(0, 1.0), (1, 2.0), (2, 3.0), (3, 4.0)]);
        assert_eq!(area(&c, &[], &[10, 20, 30, 40]), 9.0);
    }
}
