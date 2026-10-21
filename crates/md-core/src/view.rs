//! View model of the tables: the rules that decide which number the user sees (missing values,
//! NaN ordering, digits, filters). The app draws what these functions return and rebuilds a
//! [`BatchView`] only when the table generation or the view changes.

use crate::batch::{Table, group_label, group_names};
use crate::deconv::Component;
use crate::detail::component_detail;
use crate::peaks::Peak;
use crate::results::area;
use crate::stats::Biomarker;

/// What the batch table lists.
#[derive(Clone, Copy, PartialEq, Eq, Default, Debug)]
pub enum Mode {
    /// Summed quant ion integrals per run.
    #[default]
    Values,
    /// Retention time of the hit per run, in minutes.
    Times,
    /// Mean, relative standard error, hits and p-value against the blanks, per replicate group.
    Groups,
}

/// A statistic that is undefined (NaN) shows as a dash.
pub fn number(v: f64, decimals: usize) -> String {
    if v.is_finite() {
        format!("{v:.decimals$}")
    } else {
        "–".into()
    }
}

/// A mean with 6 significant digits, like the old table (Qt shows a double as `'g'` with 6), so
/// normalised fractions stay readable; without exponent.
pub fn significant(v: f64) -> String {
    if v == 0.0 || !v.is_finite() {
        return number(v, 0);
    }
    if v.abs() < 1e-10 {
        return format!("{v:.5e}");
    }
    let decimals = (5 - v.abs().log10().floor() as i32).clamp(0, 15) as usize;
    let s = format!("{v:.decimals$}");
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.').into()
    } else {
        s
    }
}

/// Retention time in minutes, two decimals.
pub fn minutes(rt_s: f64) -> String {
    format!("{:.2}", rt_s / 60.0)
}

/// A p-value with four decimals, a dash where the test is undefined.
pub fn p_fixed(p: Option<f64>) -> String {
    p.map_or("–".into(), |p| format!("{p:.4}"))
}

/// A p-value of the two-group comparison: tiny ones in scientific notation.
pub fn p_text(p: Option<f64>) -> String {
    match p {
        None => "–".into(),
        Some(p) if p < 1e-4 => format!("{p:.1e}"),
        Some(p) => format!("{p:.4}"),
    }
}

/// A missing value counts as 0, unless the table is z-scored: there 0 means "average", not
/// "not found".
pub fn missing_as_zero(wanted: bool, z_scored: bool) -> bool {
    wanted && !z_scored
}

/// The index an axis mark stands for. Plots hand the tick position over as a float; only a mark
/// within 1e-6 of a non-negative whole number is an index (else the tick is between two).
pub fn axis_index(v: f64) -> Option<usize> {
    let i = v.round();
    ((v - i).abs() < 1e-6 && i >= 0.0).then_some(i as usize)
}

/// Rows of the two-group comparison to show: p at most `p_max` (1 or more keeps all, also those
/// without p), ordered by p if `sort_by_p`, undefined p last (as 2).
pub fn biomarker_rows(rows: &[Biomarker], p_max: f64, sort_by_p: bool) -> Vec<&Biomarker> {
    let mut shown: Vec<&Biomarker> = rows
        .iter()
        .filter(|r| p_max >= 1.0 || r.p.is_some_and(|p| p <= p_max))
        .collect();
    if sort_by_p {
        shown.sort_by(|a, b| a.p.unwrap_or(2.0).total_cmp(&b.p.unwrap_or(2.0)));
    }
    shown
}

/// The ten cells of a row of the two-group comparison.
pub fn biomarker_cells(t: &Table, r: &Biomarker) -> [String; 10] {
    let s = &t.rows[r.row];
    [
        s.name.clone(),
        minutes(s.rt_s),
        r.hits[0].to_string(),
        r.hits[1].to_string(),
        significant(r.mean[0]),
        significant(r.mean[1]),
        significant(r.sd[0]),
        significant(r.sd[1]),
        format!("{:.4}", r.ratio),
        p_text(r.p),
    ]
}

/// The area and the mean S/N the component list shows for a component. The area is the match's,
/// else (no match) that of the three largest peaks; S/N is `None` for a component without peaks
/// (saved by the old program).
pub fn component_area_sn(
    matched_area: Option<f64>,
    c: &Component,
    peaks: &[Peak],
    peak_masses: &[u32],
) -> (f64, Option<f64>) {
    let d = component_detail(c, peaks, peak_masses);
    (
        matched_area.unwrap_or_else(|| area(c, &[], peak_masses)),
        (!d.peaks.is_empty()).then_some(d.mean_signal_to_noise),
    )
}

/// The sort key of a column where a missing value sorts first.
pub fn missing_first(v: Option<f64>) -> f64 {
    v.unwrap_or(f64::NEG_INFINITY)
}

/// The statistics a group has in "Group statistics", in column order.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Stat {
    Mean,
    RelStdErr,
    Hits,
    PBlank,
}

/// What a heading of the batch table is called, before translation.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Title {
    /// File name of a run.
    Run(String),
    Stat(Stat),
    /// The ANOVA column over all groups.
    Anova,
}

/// A column heading: title and the group it belongs to ("all runs" for a run without group).
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Head {
    pub title: Title,
    pub group: String,
}

/// One line of the table body.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Line {
    Total,
    Substance(usize),
    /// Quant ion `k` of substance `row`.
    Ion(usize, usize),
}

/// Headings and lines of the batch table, and the data the cells need.
#[derive(Debug, Default)]
pub struct BatchView {
    pub mode: Mode,
    pub groups: Vec<String>,
    /// Group statistics also have the blank group's p-value column.
    pub blank: Option<String>,
    pub heads: Vec<Head>,
    pub lines: Vec<Line>,
}

impl BatchView {
    /// `blanks` is the replicate group of the blank runs; `with_total` adds the total line (values
    /// only); `ion_rows` a line per quant ion (values only).
    pub fn new(t: &Table, mode: Mode, ion_rows: bool, with_total: bool, blanks: &str) -> Self {
        let groups: Vec<String> = group_names(t).into_iter().map(String::from).collect();
        let blank = (mode == Mode::Groups && groups.iter().any(|g| g == blanks))
            .then(|| blanks.to_string());
        let heads = heads(t, mode, &groups, blank.is_some());
        let lines = body_lines(t, mode, ion_rows, with_total && mode == Mode::Values);
        Self {
            mode,
            groups,
            blank,
            heads,
            lines,
        }
    }

    /// Cells of substance `row` after "Quant ions".
    pub fn cells(&self, t: &Table, row: usize) -> Vec<String> {
        let groups: Vec<&str> = self.groups.iter().map(String::as_str).collect();
        cells(t, self.mode, &groups, self.blank.as_deref(), row)
    }
}

/// What a [`BatchView`] was built from.
type Key = (u64, Mode, bool, bool, String);

/// Keeps the [`BatchView`] of the shown table: rebuilt only when the table generation or the view
/// options change, else every frame reuses it.
#[derive(Debug, Default)]
pub struct BatchViewCache {
    key: Option<Key>,
    view: BatchView,
    /// How often the view was built, for the test that frames do not rebuild it.
    pub rebuilds: usize,
}

impl BatchViewCache {
    pub fn get(
        &mut self,
        generation: u64,
        t: &Table,
        (mode, ion_rows, with_total): (Mode, bool, bool),
        blanks: &str,
    ) -> &BatchView {
        let key = (generation, mode, ion_rows, with_total, blanks.to_string());
        if self.key.as_ref() != Some(&key) {
            self.view = BatchView::new(t, mode, ion_rows, with_total, blanks);
            self.key = Some(key);
            self.rebuilds += 1;
        }
        &self.view
    }
}

/// Headings of the columns after "Quant ions".
pub fn heads(t: &Table, mode: Mode, groups: &[String], blank: bool) -> Vec<Head> {
    if mode != Mode::Groups {
        return (t.runs.iter().zip(&t.groups))
            .map(|(r, g)| Head {
                title: Title::Run(r.clone()),
                group: group_label(g).into(),
            })
            .collect();
    }
    let stats = [Stat::Mean, Stat::RelStdErr, Stat::Hits, Stat::PBlank];
    let per_group = if blank { 4 } else { 3 };
    (groups.iter())
        .flat_map(|g| {
            stats[..per_group].iter().map(move |s| Head {
                title: Title::Stat(*s),
                group: group_label(g).into(),
            })
        })
        .chain([Head {
            title: Title::Anova,
            group: String::new(),
        }])
        .collect()
}

/// The lines of the table body.
pub fn body_lines(t: &Table, mode: Mode, ion_rows: bool, with_total: bool) -> Vec<Line> {
    let ions = ion_rows && mode == Mode::Values;
    let mut lines = Vec::with_capacity(t.rows.len() + 1);
    lines.extend(with_total.then_some(Line::Total));
    for (i, r) in t.rows.iter().enumerate() {
        lines.push(Line::Substance(i));
        if ions {
            lines.extend((0..r.ions.len()).map(|k| Line::Ion(i, k)));
        }
    }
    lines
}

/// The cells of one substance in "Group statistics": per group mean, relative standard error,
/// hits and (with a blank group) the p-value against the blanks; last the ANOVA p-value over all
/// the groups, blanks included, like the old "p-value" column.
pub fn group_cells(t: &Table, row: usize, groups: &[&str], blank: Option<&str>) -> Vec<String> {
    let stats = t.group_stats(row, groups, blank);
    let mut out = Vec::new();
    for (k, s) in stats.groups.iter().enumerate() {
        out.push(significant(s.mean));
        out.push(number(s.rel_std_err, 3));
        out.push(s.hits.to_string());
        if blank.is_some() {
            out.push(p_fixed(stats.p_blank[k]));
        }
    }
    out.push(p_fixed(stats.anova));
    out
}

/// Cells of substance `row` after "Quant ions".
pub fn cells(
    t: &Table,
    mode: Mode,
    groups: &[&str],
    blank: Option<&str>,
    row: usize,
) -> Vec<String> {
    let r = &t.rows[row];
    match mode {
        Mode::Values => (r.values.iter())
            .map(|v| v.map_or("–".into(), |v| format!("{v:.0}")))
            .collect(),
        Mode::Times => (r.rts_s.iter())
            .map(|rt| rt.map_or("–".into(), minutes))
            .collect(),
        Mode::Groups => group_cells(t, row, groups, blank),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::batch::Row;

    fn table(values: &[f64], groups: &[&str]) -> Table {
        Table {
            runs: (0..values.len()).map(|i| i.to_string()).collect(),
            groups: groups.iter().map(|g| g.to_string()).collect(),
            rows: vec![Row {
                compound: 0,
                name: "s0".into(),
                rt_s: 60.0,
                ions: vec![73, 147],
                values: values.iter().map(|v| Some(*v)).collect(),
                rts_s: vec![Some(90.0); values.len()],
            }],
            ..Default::default()
        }
    }

    #[test]
    fn digits() {
        assert_eq!(significant(1_234_567.0), "1234567");
        assert_eq!(significant(-0.0123456789), "-0.0123457");
        assert_eq!(significant(f64::NAN), "–");
        assert_eq!(significant(-1e-20), "-1.00000e-20");
        assert_eq!(significant(0.0), "0");
        assert_eq!(number(f64::INFINITY, 3), "–");
        assert_eq!(minutes(90.0), "1.50");
        assert_eq!(p_text(Some(1e-5)), "1.0e-5");
        assert_eq!(p_text(Some(0.5)), "0.5000");
        assert_eq!(p_text(None), "–");
        assert_eq!(p_fixed(Some(1e-5)), "0.0000");
    }

    #[test]
    fn missing_value_counts_as_zero_unless_z_scored() {
        assert!(missing_as_zero(true, false));
        assert!(!missing_as_zero(true, true));
        assert!(!missing_as_zero(false, false));
    }

    #[test]
    fn axis_index_tolerates_float_noise_only() {
        assert_eq!(axis_index(2.0), Some(2));
        assert_eq!(axis_index(2.0 + 5e-7), Some(2));
        assert_eq!(axis_index(2.0 - 5e-7), Some(2));
        assert_eq!(axis_index(2.0 + 2e-6), None);
        assert_eq!(axis_index(1.5), None);
        assert_eq!(axis_index(-1.0), None);
        assert_eq!(axis_index(0.0), Some(0));
    }

    fn bm(row: usize, p: Option<f64>) -> Biomarker {
        Biomarker {
            row,
            hits: [1, 1],
            mean: [1.0, 1.0],
            sd: [0.0, 0.0],
            ratio: 1.0,
            p,
        }
    }

    #[test]
    fn biomarkers_filter_by_p_and_sort_undefined_p_last() {
        let rows = [
            bm(0, Some(0.5)),
            bm(1, None),
            bm(2, Some(0.01)),
            bm(3, Some(f64::NAN)),
        ];
        let order = |p_max, sort| -> Vec<usize> {
            biomarker_rows(&rows, p_max, sort)
                .iter()
                .map(|r| r.row)
                .collect()
        };
        // 1 or more keeps rows without p; below it they are dropped (NaN <= p_max is false).
        assert_eq!(order(1.0, false), [0, 1, 2, 3]);
        assert_eq!(order(0.05, false), [2]);
        assert_eq!(order(1.0, true)[..3], [2, 0, 1]);
        // NaN sorts after the undefined p (total_cmp puts positive NaN last).
        assert_eq!(order(1.0, true)[3], 3);
    }

    #[test]
    fn group_cells_count_missing_as_zero_and_an_anova_stands_last() {
        let t = table(
            &[1.0, 2.0, 3.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0],
            &["A", "A", "A", "B", "B", "B", "C", "C", "C"],
        );
        let groups = ["A", "B", "C"];
        let cells = group_cells(&t, 0, &groups, None);
        assert_eq!(cells.len(), 3 * 3 + 1);
        assert_eq!(cells.last().unwrap(), "0.0066");
        let cells = group_cells(&t, 0, &groups, Some("C"));
        assert_eq!(cells.len(), 3 * 4 + 1);
        assert_eq!(cells.last().unwrap(), "0.0066");
        let g: Vec<String> = groups.map(String::from).to_vec();
        assert_eq!(heads(&t, Mode::Groups, &g, true).len(), cells.len());
        assert_eq!(heads(&t, Mode::Groups, &g, false).len(), 3 * 3 + 1);
    }

    #[test]
    fn group_cells_are_readable() {
        let t = table(&[1.0, 2.0, 10.0], &["A"; 3]);
        assert_eq!(
            group_cells(&t, 0, &["A"], None)[..3],
            ["4.33333", "0.657", "3"]
        );
        // One group: nothing to compare, p = 1 like the old code.
        assert_eq!(group_cells(&t, 0, &["A"], None)[3], "1.0000");
    }

    #[test]
    fn lines_follow_the_mode() {
        let t = table(&[1.0], &[""]);
        assert_eq!(
            body_lines(&t, Mode::Values, false, false),
            [Line::Substance(0)]
        );
        assert_eq!(
            body_lines(&t, Mode::Values, true, true),
            [
                Line::Total,
                Line::Substance(0),
                Line::Ion(0, 0),
                Line::Ion(0, 1)
            ]
        );
        assert_eq!(body_lines(&t, Mode::Times, true, false).len(), 1);
    }

    #[test]
    fn view_has_blank_only_in_group_mode_and_total_only_for_values() {
        let t = table(&[1.0, 2.0], &["A", "BLANK"]);
        let v = BatchView::new(&t, Mode::Groups, false, true, "BLANK");
        assert_eq!(v.blank.as_deref(), Some("BLANK"));
        assert_eq!(v.heads.len(), 2 * 4 + 1);
        assert_eq!(v.lines, [Line::Substance(0)]);
        let v = BatchView::new(&t, Mode::Values, false, true, "BLANK");
        assert_eq!(v.blank, None);
        assert_eq!(v.lines, [Line::Total, Line::Substance(0)]);
        assert_eq!(v.cells(&t, 0), ["1", "2"]);
        assert_eq!(
            BatchView::new(&t, Mode::Times, false, true, "BLANK").cells(&t, 0),
            ["1.50", "1.50"]
        );
    }

    #[test]
    fn the_view_is_rebuilt_only_when_generation_or_view_change() {
        let t = table(&[1.0, 2.0], &["A", "BLANK"]);
        let mut c = BatchViewCache::default();
        for _ in 0..5 {
            c.get(1, &t, (Mode::Values, false, false), "BLANK");
        }
        assert_eq!(c.rebuilds, 1);
        c.get(2, &t, (Mode::Values, false, false), "BLANK");
        c.get(2, &t, (Mode::Groups, false, false), "BLANK");
        c.get(2, &t, (Mode::Groups, true, false), "BLANK");
        c.get(2, &t, (Mode::Groups, true, false), "OTHER");
        assert_eq!(c.rebuilds, 5);
    }

    #[test]
    fn missing_sorts_first() {
        assert!(missing_first(None) < missing_first(Some(-1e300)));
    }
}
