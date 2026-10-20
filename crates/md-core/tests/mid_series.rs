//! `mid::determine` against `gcms::MIDCompoundSeries` of MetaboliteDetector 3.3, run through
//! `docs/reverse-engineering/tools/isotopes/old_midseries.cpp` on synthetic integrals. Two cases are
//! inline; `MD_MID_SERIES_DIR` (cases.txt and old.tsv from make-midseries.sh) checks all of them, and
//! with `MD_LBR_SAMPLE` (lib_aa.lbr) the library spectrum comes from the old library itself. Group
//! statistics and the text exports of the result window go against old_midstats.cpp the same way
//! (one inline case, all with old-stats.txt and the library).

use std::collections::BTreeMap;
use std::path::PathBuf;

use md_core::batch::{Found, IonPeak, RunResult};
use md_core::isotope::IsotopeTable;
use md_core::library::LibraryCompound;
use md_core::mid::{FragmentMid, Method, Mids, determine, group_stat};

const CASES: &str = "\
case lib_n 0 N 1 2
frag 71 74 - 0
int 0 71 145613.305
int 0 72 122308.913
int 0 73 28085.554
int 0 74 11210.64
int 1 71 88489.312
int 1 72 84550.038
int 1 73 11421.418
int 1 74 82027.455
end
case form_none 0 - 0 2
frag 71 74 C4H7O 2
int 0 71 391170.364
int 0 72 192938.593
int 0 73 84893.295
int 0 74 58258.149
int 1 71 19507.217
int 1 72 38439.034
int 1 73 12330.078
int 1 74 14935.102
end
";

const OLD: &str = "\
norm	lib_n	71	729
norm	lib_n	72	52
norm	lib_n	73	0
norm	lib_n	74	0
mid	lib_n	0	71	0	0.50778202686474416	-1	1
mid	lib_n	0	71	1	0.35574295086316721	-1	1
mid	lib_n	0	72	0	0.38895530251679089	-1	1
mid	lib_n	0	72	1	0.31345163561847994	-1	1
mid	lib_n	0	73	0	0.071043797884228516	-1	1
mid	lib_n	0	73	1	0.024467743749871401	-1	1
mid	lib_n	0	74	0	0.034161709073323811	-1	1
mid	lib_n	0	74	1	0.3248101955382976	-1	1
norm	form_none	71	729
norm	form_none	72	52
mid	form_none	0	71	0	0.56384211512365345	0.18305373215719994	0.86774409631138538
mid	form_none	0	71	1	0.24004312037735365	0.29858077845914321	0.20694827790615045
mid	form_none	0	72	0	0.25774159534516505	0.18305373303442762	0.86774409631138538
mid	form_none	0	72	1	0.46834809386819642	0.29858077988999815	0.20694827790615045
";

struct Case {
    name: String,
    compound: usize,
    tracer: String,
    library: bool,
    runs: usize,
    frags: Vec<(u32, u32, String, u32)>,
    ints: Vec<(usize, u32, f64)>,
}

fn cases(text: &str) -> Vec<Case> {
    let mut out: Vec<Case> = Vec::new();
    for line in text.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let dash = |s: &str| if s == "-" { String::new() } else { s.into() };
        match w.first().copied() {
            Some("case") => out.push(Case {
                name: w[1].into(),
                compound: w[2].parse().unwrap(),
                tracer: dash(w[3]),
                library: w[4] == "1",
                runs: w[5].parse().unwrap(),
                frags: Vec::new(),
                ints: Vec::new(),
            }),
            Some("frag") => out.last_mut().unwrap().frags.push((
                w[1].parse().unwrap(),
                w[2].parse().unwrap(),
                dash(w[3]),
                w[4].parse().unwrap(),
            )),
            Some("int") => out.last_mut().unwrap().ints.push((
                w[1].parse().unwrap(),
                w[2].parse().unwrap(),
                w[3].parse().unwrap(),
            )),
            _ => {}
        }
    }
    out
}

fn old_table() -> IsotopeTable {
    let text = include_str!("../../../docs/reverse-engineering/tools/isotopes/old-isotopes.tsv");
    let rows: Vec<(&str, u32, f64, f64)> = text
        .lines()
        .filter_map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Some((
                f[0],
                f[1].parse().ok()?,
                f[2].parse().ok()?,
                f[3].trim().parse().ok()?,
            ))
        })
        .collect();
    IsotopeTable::new(rows)
}

/// The MIDs of one case and its library compound; `lib` gives the library compounds, else the old
/// output's `norm` lines make up the spectrum.
fn mids_of(
    k: &Case,
    old: &str,
    lib: Option<&[LibraryCompound]>,
    table: &IsotopeTable,
) -> (Mids, LibraryCompound) {
    let mut c = match lib {
        Some(lib) => lib[k.compound].clone(),
        None => LibraryCompound::new(&k.name),
    };
    let norms: BTreeMap<u32, i32> = old
        .lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|f| f[0] == "norm" && f[1] == k.name)
        .map(|f| (f[2].parse().unwrap(), f[3].parse().unwrap()))
        .collect();
    if lib.is_none() {
        c.spectrum.masses = norms.keys().map(|m| *m as f32).collect();
        c.spectrum.values = norms.values().copied().collect();
    }
    // The old program normalises the stored spectrum to 1000; check that it is stored that way.
    for (m, v) in &norms {
        let s = &c.spectrum;
        let got = (s.masses.iter().position(|x| *x == *m as f32)).map_or(0, |i| s.values[i]);
        assert_eq!(got, *v, "{} ion {m}", k.name);
    }
    c.formulas = k.frags.iter().map(|f| (f.0, f.1, f.2.clone())).collect();
    c.mid_sizes = k.frags.iter().map(|f| (f.0, f.1, f.3)).collect();
    let runs: Vec<RunResult> = (0..k.runs)
        .map(|r| RunResult {
            path: PathBuf::from(format!("run{r}")),
            found: vec![Found {
                compound: 0,
                score: 1.0,
                rt_s: 0.0,
                ions: (k.ints.iter())
                    .filter(|i| i.0 == r)
                    .map(|i| IonPeak {
                        mass: i.1,
                        integral: i.2,
                        quality: 0.0,
                        signal_to_noise: 0.0,
                    })
                    .collect(),
                rescanned: false,
                rescued: Vec::new(),
            }],
            unknown: Vec::new(),
            calibration: None,
        })
        .collect();
    let method = if k.library {
        Method::Library
    } else {
        Method::Formula
    };
    let mids = determine(&runs, &[c.clone()], &[], &[], method, &k.tracer, table);
    assert!(mids.errors.is_empty(), "{}: {:?}", k.name, mids.errors);
    (mids, c)
}

/// Compares every case with the old `mid` lines; a NaN MID (run without integrals) means no fit.
fn check(cases_text: &str, old: &str, lib: Option<&[LibraryCompound]>) -> usize {
    let table = old_table();
    let mut compared = 0;
    for k in cases(cases_text) {
        let frags = &mids_of(&k, old, lib, &table).0.compounds[0].fragments;
        let mut rows = 0;
        for line in old.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f[0] != "mid" || f[1] != k.name {
                continue;
            }
            let (frag, ion, run): (usize, u32, usize) = (
                f[2].parse().unwrap(),
                f[3].parse().unwrap(),
                f[4].parse().unwrap(),
            );
            let [mid, ci, r2]: [f64; 3] = [5, 6, 7].map(|i| f[i].parse().unwrap());
            compared += 1;
            if mid.is_nan() {
                assert!(frags[frag].fits[run].is_none(), "{line}");
                continue;
            }
            rows += 1;
            let fit = frags[frag].fits[run].as_ref().expect("fit");
            let n = (ion - frags[frag].lo) as usize;
            let close = |a: f64, e: f64| (a - e).abs() <= 1e-9 * e.abs().max(1.0);
            assert!(close(fit.mids[n], mid), "{line}: mid {}", fit.mids[n]);
            assert!(
                close(fit.cis[n].unwrap_or(-1.0), ci),
                "{line}: ci {:?}",
                fit.cis[n]
            );
            assert!(close(fit.r2.unwrap_or(-1.0), r2), "{line}: r2 {:?}", fit.r2);
        }
        let ours: usize = (frags.iter())
            .map(|f| f.fits.iter().flatten().map(|x| x.mids.len()).sum::<usize>())
            .sum();
        assert_eq!(ours, rows, "{}: number of values", k.name);
    }
    compared
}

/// Compares group statistics and both text exports with the output of old_midstats.cpp (`stat`
/// lines and `begin data|stats NAME` … `end` blocks); `old` gives the spectra as in [`check`].
fn check_exports(cases_text: &str, old: &str, stats: &str, lib: &[LibraryCompound]) -> usize {
    let table = old_table();
    let mut compared = 0;
    for k in cases(cases_text) {
        let (mids, c) = mids_of(&k, old, Some(lib), &table);
        let method = if k.library {
            Method::Library
        } else {
            Method::Formula
        };
        let isotopomers: Vec<(&FragmentMid, usize)> = (mids.compounds[0].fragments.iter())
            .flat_map(|f| (0..f.isotopomers).map(move |n| (f, n)))
            .collect();
        let mut lines = 0;
        for line in stats.lines() {
            let f: Vec<&str> = line.split('\t').collect();
            if f[0] != "stat" || f[1] != k.name {
                continue;
            }
            let (frag, n) = isotopomers[f[2].parse::<usize>().unwrap()];
            let (mean, se) = group_stat(frag, n);
            for (got, want) in [(mean, f[3]), (se, f[4])] {
                let want: f64 = want.parse().unwrap();
                let same = (got.is_nan() && want.is_nan())
                    || (got - want).abs() <= 1e-9 * want.abs().max(1.0);
                assert!(same, "{line}: {mean} {se}");
            }
            lines += 1;
        }
        assert_eq!(lines, isotopomers.len(), "{}: number of statistics", k.name);
        let block = |what: &str| {
            let start = format!("begin\t{what}\t{}\n", k.name);
            let rest = &stats[stats.find(&start).expect(&start) + start.len()..];
            &rest[..rest.find("\nend\n").unwrap()]
        };
        let lib = [c];
        let columns: Vec<String> = (0..k.runs).map(|r| format!("run{r}")).collect();
        assert_eq!(
            mids.run_text(&lib, method, &columns),
            block("data"),
            "{} data",
            k.name
        );
        assert_eq!(mids.stats_text(&lib), block("stats"), "{} stats", k.name);
        compared += lines;
    }
    compared
}

#[test]
fn matches_old_mid_compound_series() {
    assert_eq!(check(CASES, OLD, None), 12);
}

#[test]
fn matches_old_mid_compound_series_on_all_cases() {
    let Some(dir) = std::env::var_os("MD_MID_SERIES_DIR").map(PathBuf::from) else {
        eprintln!("MD_MID_SERIES_DIR not set; skipped");
        return;
    };
    let read = |f: &str| std::fs::read_to_string(dir.join(f)).unwrap();
    let lib = std::env::var_os("MD_LBR_SAMPLE").map(|p| md_core::library::read_lbr(p).unwrap());
    let n = check(&read("cases.txt"), &read("old.tsv"), lib.as_deref());
    eprintln!("{n} values compared");
    assert!(n > 0);
    let Some(lib) = lib else {
        eprintln!("MD_LBR_SAMPLE not set; exports skipped");
        return;
    };
    let n = check_exports(
        &read("cases.txt"),
        &read("old.tsv"),
        &read("old-stats.txt"),
        &lib,
    );
    eprintln!("{n} statistics and their exports compared");
    assert!(n > 0);
}

/// A formula case whose second run lacks the substance (NaN in the old program).
const GAP_CASE: &str = "\
case form_gap 0 C 0 4
frag 71 74 C4H7O 3
int 0 71 462384.312
int 0 72 250263.652
int 0 73 124338.632
int 0 74 55029.149
int 2 71 315563.847
int 2 72 159253.979
int 2 73 164264.62
int 2 74 129922.909
int 3 71 518241.638
int 3 72 167105.913
int 3 73 110235.734
int 3 74 95507.873
end
";

const GAP_STATS: &str = concat!(
    "stat\tform_gap\t0\t0.5276592366180165\t0.05239675641355724\n",
    "stat\tform_gap\t1\t0.20993057988932737\t0.029539106749154075\n",
    "stat\tform_gap\t2\t0.15698711796564299\t0.029285145438480854\n",
    "begin\tdata\tform_gap\n",
    "\trun0\trun1\trun2\trun3\t\tCI (run0)\tCI (run1)\tCI (run2)\tCI (run3)",
    "\t\tR2 (run0)\tR2 (run1)\tR2 (run2)\tR2 (run3)\n",
    "Hexatriacontane;3599.99;1.04257e+06\n",
    "71 (M0)\t0.543308\tnan\t0.430099\t0.609571\t\t0.187033\tnan\t0.527726\t0.335421",
    "\t\t0.972136\tnan\t0.24501\t0.928015\n",
    "72 (M1)\t0.266686\tnan\t0.195758\t0.167348\t\t0.185199\tnan\t0.522551\t0.332132",
    "\t\t0.972136\tnan\t0.24501\t0.928015\n",
    "73 (M2)\t0.133799\tnan\t0.21516\t0.122002\t\t0.183067\tnan\t0.516536\t0.328309",
    "\t\t0.972136\tnan\t0.24501\t0.928015\n",
    "\n\n\nend\n",
    "begin\tstats\tform_gap\n",
    "\tGL\t\tStdErr (GL)\tMetabolite name\tFragment ion\tMass isotopomer\n",
    "Hexatriacontane;3599.99;1.04257e+06\n",
    "71 (M0)\t0.527659\t\t0.0523968\tHexatriacontane\t71\tM0\n",
    "72 (M1)\t0.209931\t\t0.0295391\tHexatriacontane\t71\tM1\n",
    "73 (M2)\t0.156987\t\t0.0292851\tHexatriacontane\t71\tM2\n",
    "\n\n\nend\n",
);

#[test]
fn matches_old_statistics_and_exports() {
    // Substance 0 of lib_aa.lbr as far as the exports show it.
    let mut c = LibraryCompound::new("Hexatriacontane");
    c.spectrum.ri = 3599.99;
    c.spectrum.rt_ms = 1042570.0;
    assert_eq!(check_exports(GAP_CASE, "", GAP_STATS, &[c]), 3);
}
