//! NIST text import and search on `data/synthetic-nist.sdf` (invented names and peaks, no NIST data).
//! The expected values are what the old 3.3 binary printed for the same file and query
//! (`docs/reverse-engineering/tools/nist/make-nist-dump.sh`).

use md_core::library::{import_nist, is_nist};
use md_core::search::{PreparedLibrary, SearchSettings, library_hits};

const ALPHA: [(u32, u32); 12] = [
    (71, 100),
    (72, 12),
    (77, 40),
    (78, 9),
    (91, 80),
    (92, 15),
    (105, 60),
    (106, 7),
    (119, 50),
    (133, 40),
    (147, 30),
    (161, 20),
];

fn fixture() -> std::path::PathBuf {
    std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/synthetic-nist.sdf")
}

#[test]
fn imports_like_the_old_import_nist() {
    assert!(is_nist(fixture()));
    let lib = import_nist(fixture()).unwrap();
    let names: Vec<_> = lib.iter().map(|c| c.name.as_str()).collect();
    // delta has only peaks at m/z 50 and 60 and is skipped ("No intensity, skipping...").
    assert_eq!(
        names,
        [
            "Synthetic alpha",
            "Synthetic beta",
            "Synthetic gamma",
            "Synthetic \"quoted\" alpha twin",
            "Synthetic eta",
            "Synthetic theta",
            "Synthetic iota",
            "Synthetic kappa",
            "Synthetic lambda",
        ]
    );
    let a = &lib[0];
    assert_eq!(a.spectrum.ri, 1450.0);
    assert_eq!(a.spectrum.rt_ms, -1.0);
    assert_eq!(a.spectrum.masses.len(), 12);
    assert_eq!(a.spectrum.values[..2], [100, 12]);
    for (k, v) in [
        ("ID", "9001"),
        ("CASNO", "12345"),
        ("NISTNO", "1"),
        ("FORMULA", "C7H8"),
        ("MW", "92"),
        ("RI ESTIMATION ERROR", "20"),
    ] {
        assert_eq!(a.features[k], v);
    }
    assert!(a.sources[0].ends_with("synthetic-nist.sdf"));
    // Absent RI: the old -1.
    assert_eq!(lib[1].spectrum.ri, -1.0);
    // Old table `spectrum` holds 73 peak rows for these records.
    assert_eq!(
        lib.iter().map(|c| c.spectrum.masses.len()).sum::<usize>(),
        73 + 8
    );
}

#[test]
fn search_uses_the_old_fragment_prefilter() {
    let lib = import_nist(fixture()).unwrap();
    let s = SearchSettings::default();
    let prepared = PreparedLibrary::new(&lib, &s).with_fragment_filter(true);
    let hits = library_hits(&ALPHA, Some(1455.0), &prepared, &s);
    // Old: alpha and twin (printed twice as alpha, see the notes), then kappa. The old scores are
    // single-precision floats in part, so they agree to about 1e-8.
    let want = [
        ("Synthetic alpha", 1.0, 0.969233234, 0.989637399),
        (
            "Synthetic \"quoted\" alpha twin",
            1.0,
            0.969233234,
            0.989637399,
        ),
        ("Synthetic kappa", 0.773023519, 1.0, 0.842290311),
    ];
    assert_eq!(hits.len(), want.len());
    for (h, (name, spec, ri, overall)) in hits.iter().zip(want) {
        assert_eq!(lib[h.compound].name, name);
        assert!((h.spec_sim - spec).abs() < 1e-7, "{name} {h:?}");
        assert!((h.ri_sim.unwrap() - ri).abs() < 1e-7, "{name} {h:?}");
        assert!((h.overall - overall).abs() < 1e-7, "{name} {h:?}");
    }

    // Old potential hits: alpha, twin, eta, theta, kappa. theta then fails the RI gate.
    let all = SearchSettings { cutoff: 0.0, ..s };
    let names = |p: &PreparedLibrary| -> Vec<&str> {
        let mut n: Vec<_> = library_hits(&ALPHA, Some(1455.0), p, &all)
            .iter()
            .map(|h| lib[h.compound].name.as_str())
            .collect();
        n.sort();
        n
    };
    assert_eq!(
        names(&prepared),
        [
            "Synthetic \"quoted\" alpha twin",
            "Synthetic alpha",
            "Synthetic eta",
            "Synthetic kappa"
        ]
    );
    // Without the pre-filter iota (3 shared fragments), gamma (none) and others are scored too.
    let plain = PreparedLibrary::new(&lib, &all);
    assert!(names(&plain).contains(&"Synthetic iota"));
    assert!(names(&plain).contains(&"Synthetic gamma"));
}
