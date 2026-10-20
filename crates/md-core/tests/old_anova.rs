//! `stats::anova` against `tools::Statistics::oneWayANOVA` of MetaboliteDetector 3.3, run through
//! `docs/reverse-engineering/tools/anova/old_anova.cpp`. A few cases are inline;
//! `MD_ANOVA_DIR` (cases.txt and old.tsv from make-anova.sh) checks all of them.

use std::path::PathBuf;

use md_core::stats::anova;

const CASES: &str = "\
case hand
group 1 2 3
group 2 3 4
group 5 6 7
end
case unequal
group 1 2 3
group 2 4
group 5 6 7 8
group 3 3 4
end
case zeros
group 0 0 5
group 0 0 0
group 0 1 0
end
case singletons_and_pair
group 1
group 2 3
end
case separated
group 1 1
group 2 2
end
case equal
group 3 3
group 3 3
end
case n_equals_k
group 1
group 2
group 3
end
";

const OLD: &str = "\
anova	hand	0.0065917968750000069
anova	unequal	0.0031605828418441971
anova	zeros	0.489078108912819
anova	singletons_and_pair	0.33333333333333304
anova	separated	0
anova	equal	nan
error	n_equals_k
";

fn parse_cases(text: &str) -> Vec<(String, Vec<Vec<f64>>)> {
    let mut cases = Vec::new();
    for line in text.lines() {
        let mut w = line.split_whitespace();
        match w.next() {
            Some("case") => cases.push((w.next().unwrap().to_string(), Vec::new())),
            Some("group") => cases
                .last_mut()
                .unwrap()
                .1
                .push(w.map(|v| v.parse().unwrap()).collect()),
            _ => {}
        }
    }
    cases
}

/// What the old program printed for each case: `Some(p)`, `None` for nan (0 / 0) and for a case in
/// which it aborts (its GSL error handler: one or no degrees of freedom).
fn parse_old(text: &str) -> Vec<(String, Option<f64>)> {
    (text.lines())
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let p = if f[0] == "anova" {
                f[2].parse().ok()
            } else {
                None
            };
            (f[1].to_string(), p.filter(|p: &f64| p.is_finite()))
        })
        .collect()
}

fn check(cases: &str, old: &str) {
    let (cases, old) = (parse_cases(cases), parse_old(old));
    assert_eq!(cases.len(), old.len());
    for ((name, groups), (old_name, want)) in cases.iter().zip(&old) {
        assert_eq!(name, old_name);
        let got = anova(groups);
        // Where the old program stops or divides 0 by 0 the port says "no p-value".
        match (got, want) {
            (Some(g), Some(w)) => assert!(
                (g - w).abs() <= 1e-9 * w.abs().max(1e-300) + 1e-14,
                "{name}: {g} against {w}"
            ),
            (None, None) => {}
            _ => panic!("{name}: {got:?} against {want:?}"),
        }
    }
}

#[test]
fn matches_the_old_one_way_anova() {
    check(CASES, OLD);
}

/// `getANOVAPValue` (0x63ac2c) does not call `oneWayANOVA` for one group or none, which would stop it.
#[test]
fn one_group_or_none_is_p_1() {
    assert_eq!(anova(&[vec![1.0, 2.0, 3.0]]), Some(1.0));
    assert_eq!(anova(&[]), Some(1.0));
}

#[test]
fn all_cases_of_the_old_program() {
    let Some(dir) = std::env::var_os("MD_ANOVA_DIR").map(PathBuf::from) else {
        return;
    };
    let read = |f: &str| std::fs::read_to_string(dir.join(f)).unwrap();
    // The cases for which getANOVAPValue answers before calling oneWayANOVA are not in the old output.
    let cases = read("cases.txt");
    let old = read("old.tsv");
    let skip = ["one_group", "none"];
    let keep = |name: &str| !skip.contains(&name);
    let cases: String = {
        let mut out = String::new();
        let mut on = true;
        for l in cases.lines() {
            if let Some(n) = l.strip_prefix("case ") {
                on = keep(n);
            }
            if on {
                out += l;
                out += "\n";
            }
        }
        out
    };
    let old: String = old
        .lines()
        .filter(|l| keep(l.split('\t').nth(1).unwrap()))
        .map(|l| format!("{l}\n"))
        .collect();
    check(&cases, &old);
}
