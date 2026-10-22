//! `mid::determine` against `gcms::MIDCompoundSeries` of MetaboliteDetector 3.3, run through
//! `docs/reverse-engineering/tools/isotopes/old_midseries.cpp` on synthetic integrals. Two cases are
//! inline; `MD_MID_SERIES_DIR` (cases.txt and old.tsv from make-midseries.sh) checks all of them, and
//! with `MD_LBR_SAMPLE` (lib_aa.lbr) the library spectrum comes from the old library itself. Group
//! statistics and the text exports of the result window go against old_midstats.cpp the same way
//! (one inline case, all with old-stats.txt and the library), the group methods 0 and 1 and their
//! export against old_midgroups.cpp (three inline cases, all with groups.txt and old-groups.txt).

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
    method: Method,
    /// Runs of the groups GL and GUL (group methods only).
    labelled: Vec<usize>,
    unlabelled: Vec<usize>,
    runs: usize,
    frags: Vec<(u32, u32, String, u32)>,
    ints: Vec<(usize, u32, f64)>,
}

/// Cases of old_midseries.cpp (`case NAME COMPOUND TRACER FLAG RUNS`) or, with `group`, of
/// old_midgroups.cpp (`case NAME COMPOUND METHOD M1 GL GUL RUNS`).
fn cases(text: &str, group: bool) -> Vec<Case> {
    let mut out: Vec<Case> = Vec::new();
    for line in text.lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        let dash = |s: &str| if s == "-" { String::new() } else { s.into() };
        let runs =
            |s: &str| -> Vec<usize> { s.split(',').filter_map(|r| r.parse().ok()).collect() };
        match w.first().copied() {
            Some("case") if group => out.push(Case {
                name: w[1].into(),
                compound: w[2].parse().unwrap(),
                tracer: String::new(),
                method: match w[3] {
                    "0" => Method::Uncorrected,
                    _ => Method::Reference {
                        m1: w[4].parse().unwrap(),
                    },
                },
                labelled: runs(w[5]),
                unlabelled: runs(w[6]),
                runs: w[7].parse().unwrap(),
                frags: Vec::new(),
                ints: Vec::new(),
            }),
            Some("case") => out.push(Case {
                name: w[1].into(),
                compound: w[2].parse().unwrap(),
                tracer: dash(w[3]),
                method: match w[4] {
                    "1" => Method::Library,
                    _ => Method::Formula,
                },
                labelled: Vec::new(),
                unlabelled: Vec::new(),
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
            input: r,
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
                alternatives: Vec::new(),
            }],
            unknown: Vec::new(),
            calibration: None,
        })
        .collect();
    let mids = determine(
        &runs,
        &[c.clone()],
        &k.labelled,
        &k.unlabelled,
        k.method,
        &k.tracer,
        table,
    );
    assert!(mids.errors.is_empty(), "{}: {:?}", k.name, mids.errors);
    (mids, c)
}

/// Compares every case with the old `mid` lines; a NaN MID (run without integrals) means no fit.
fn check(cases_text: &str, old: &str, lib: Option<&[LibraryCompound]>) -> usize {
    let table = old_table();
    let mut compared = 0;
    for k in cases(cases_text, false) {
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
    for k in cases(cases_text, false) {
        let (mids, c) = mids_of(&k, old, Some(lib), &table);
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
            mids.run_text(&lib, k.method, &columns),
            block("data"),
            "{} data",
            k.name
        );
        assert_eq!(mids.stats_text(&lib), block("stats"), "{} stats", k.name);
        compared += lines;
    }
    compared
}

/// Compares the group methods 0 and 1 with the output of old_midgroups.cpp: `gmid`/`gr2` lines
/// per fragment key and the export between `begin group NAME` and `end`. The old export writes no
/// rows for a reversed key `(b, a)`; the port does (algorithms.md §7g.7), so those rows are left
/// out before comparing. A case the old program fails on (`error`, method 0 without any labelled
/// row) must give no MID here.
fn check_groups(cases_text: &str, old: &str, lib: Option<&[LibraryCompound]>) -> usize {
    let table = old_table();
    let mut compared = 0;
    for k in cases(cases_text, true) {
        let (mids, c) = mids_of(&k, "", lib, &table);
        let fields = |tag: &str| -> Vec<Vec<&str>> {
            (old.lines())
                .map(|l| l.split('\t').collect::<Vec<_>>())
                .filter(|f| f[0] == tag && f[1] == k.name)
                .collect()
        };
        let fits: Vec<Option<_>> = (mids.compounds.iter())
            .flat_map(|c| &c.fragments)
            .map(|f| f.fits[0].as_ref())
            .collect();
        // The port reports no errors for the group methods, so the old one must not either.
        assert!(fields("gerror").is_empty(), "{}", k.name);
        if !fields("error").is_empty() {
            assert!(fits.iter().all(Option::is_none), "{}: {fits:?}", k.name);
            continue;
        }
        let frags = &mids.compounds[0].fragments;
        let close = |a: f64, e: f64| (a - e).abs() <= 1e-9 * e.abs().max(1.0);
        for f in fields("gr2") {
            let (frag, r2): (usize, f64) = (f[2].parse().unwrap(), f[3].parse().unwrap());
            let gmid: Vec<[f64; 2]> = (fields("gmid").iter())
                .filter(|g| g[2] == f[2])
                .map(|g| [4, 5].map(|i| g[i].parse().unwrap()))
                .collect();
            compared += gmid.len();
            let Some(fit) = &frags[frag].fits[0] else {
                // The old solver's empty result.
                assert_eq!((gmid, r2), (vec![[-1.0, -1.0]], -1.0), "{} {frag}", k.name);
                continue;
            };
            assert_eq!(fit.mids.len(), gmid.len(), "{} {frag}", k.name);
            assert!(
                close(fit.r2.unwrap_or(-1.0), r2),
                "{} {frag}: r2 {:?}",
                k.name,
                fit.r2
            );
            for (n, [mid, ci]) in gmid.into_iter().enumerate() {
                assert!(
                    close(fit.mids[n], mid),
                    "{} {frag} M{n}: {}",
                    k.name,
                    fit.mids[n]
                );
                let got = fit.cis[n].unwrap_or(-1.0);
                assert!(close(got, ci), "{} {frag} M{n}: ci {got}", k.name);
            }
        }
        let mut shown = mids.clone();
        for (f, key) in shown.compounds[0]
            .fragments
            .iter_mut()
            .zip(fragment_keys(&k))
        {
            if key.0 > key.1 {
                f.fits.clear();
            }
        }
        let start = format!("begin\tgroup\t{}\n", k.name);
        let rest = &old[old.find(&start).expect(&start) + start.len()..];
        assert_eq!(
            shown.group_text(&[c]),
            &rest[..rest.find("\nend\n").unwrap()],
            "{} export",
            k.name
        );
    }
    compared
}

/// Fragment keys of a case in the old map order.
fn fragment_keys(k: &Case) -> Vec<(u32, u32)> {
    let mut keys: Vec<(u32, u32)> = k.frags.iter().map(|f| (f.0, f.1)).collect();
    keys.sort_unstable();
    keys
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
    let n = check_groups(&read("groups.txt"), &read("old-groups.txt"), Some(&lib));
    eprintln!("{n} values of the group methods and their exports compared");
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

/// Group cases: method 0 with a reversed key, method 1 with M1 correction, and method 1 whose
/// reference run lacks the substance (the old solver's empty result).
const GROUP_CASES: &str = "\
case m0_rev 0 0 0 0,1,2 - 3
frag 99 97 - 0
frag 57 58 - 0
int 0 97 6840.373
int 0 98 5156.745
int 0 99 6993.943
int 1 97 118349.324
int 1 98 76121.941
int 1 99 34373.447
int 2 97 436361.024
int 2 98 219773.602
int 2 99 43200.077
int 0 57 614958.851
int 0 58 71311.902
int 1 57 421034.15
int 1 58 66007.769
int 2 57 592328.252
int 2 58 316825.716
end
case m1_m1 0 1 0.011 0,1,2 3,4 5
frag 57 60 - 0
int 0 57 719705.585
int 0 58 95354.932
int 0 59 181884.679
int 0 60 140326.155
int 1 57 394120.126
int 1 58 220766.366
int 1 59 127405.981
int 1 60 112701.757
int 2 57 112625.419
int 2 58 93254.73
int 2 59 8949.519
int 2 60 28715.831
int 3 57 102520.481
int 3 58 27542.047
int 3 59 40726.816
int 3 60 12921.254
int 4 57 159923.825
int 4 58 137828.143
int 4 59 108964.203
int 4 60 50546.176
end
case m1_noref 0 1 0.011 0,1 2 3
frag 57 59 - 0
int 0 57 87522.67
int 0 58 101263.838
int 0 59 28742.928
int 1 57 308824.623
int 1 58 180046.262
int 1 59 61481.868
end
";

const OLD_GROUPS: &str = concat!(
    "gmid\tm0_rev\t0\t0\t0.80402526214245984\t0.15052587657174457\n",
    "gmid\tm0_rev\t0\t1\t0.19597473785754016\t0.15052587657174457\n",
    "gr2\tm0_rev\t0\t-1\n",
    "gmid\tm0_rev\t1\t0\t0.50043832771605301\t0.15014275629775259\n",
    "gmid\tm0_rev\t1\t1\t0.30614404639917336\t0.035474163199280859\n",
    "gmid\tm0_rev\t1\t2\t0.19341762588477365\t0.17851601649619903\n",
    "gr2\tm0_rev\t1\t-1\n",
    "begin\tgroup\tm0_rev\n",
    "Hexatriacontane_57 (M0)\t0.804025\t0.150526\t-1\n",
    "Hexatriacontane_58 (M1)\t0.195975\t0.150526\t-1\n",
    "\n\n\nend\n",
    "gmid\tm1_m1\t0\t0\t1.0397111386668942\t0.13211151180386915\n",
    "gmid\tm1_m1\t0\t1\t0.029990961785648917\t0.14322119978946884\n",
    "gmid\tm1_m1\t0\t2\t-0.25407307419108482\t0.14618003078698663\n",
    "gmid\tm1_m1\t0\t3\t0.16184218124501695\t0.15002564642316882\n",
    "gr2\tm1_m1\t0\t0.64651089044429044\n",
    "begin\tgroup\tm1_m1\n",
    "Hexatriacontane_57 (M0)\t1.03971\t0.132112\t0.646511\n",
    "Hexatriacontane_58 (M1)\t0.029991\t0.143221\t0.646511\n",
    "Hexatriacontane_59 (M2)\t-0.254073\t0.14618\t0.646511\n",
    "Hexatriacontane_60 (M3)\t0.161842\t0.150026\t0.646511\n",
    "\n\nend\n",
    "gmid\tm1_noref\t0\t0\t-1\t-1\n",
    "gr2\tm1_noref\t0\t-1\n",
    "begin\tgroup\tm1_noref\n",
    "Hexatriacontane_57 (M0)\t-1\t-1\t-1\n",
    "\n\nend\n",
);

#[test]
fn matches_old_group_methods_and_export() {
    let c = LibraryCompound::new("Hexatriacontane");
    assert_eq!(check_groups(GROUP_CASES, OLD_GROUPS, Some(&[c])), 10);
}

#[test]
fn matches_old_statistics_and_exports() {
    // Substance 0 of lib_aa.lbr as far as the exports show it.
    let mut c = LibraryCompound::new("Hexatriacontane");
    c.spectrum.ri = 3599.99;
    c.spectrum.rt_ms = 1042570.0;
    assert_eq!(check_exports(GAP_CASE, "", GAP_STATS, &[c]), 3);
}
