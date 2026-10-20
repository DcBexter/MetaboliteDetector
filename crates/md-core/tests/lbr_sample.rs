//! Reads the sample library `lib_aa.lbr` from the old 3.3 package (never committed).
//! Run with `MD_LBR_SAMPLE=path/to/lib_aa.lbr cargo test -p md-core --test lbr_sample -- --nocapture`;
//! the printed dump matches `docs/reverse-engineering/tools/dump_lbr.py` field by field.
//! Passes without checking anything when the variable is not set.

use md_core::library::{
    LibraryCompound, LibraryPeak, parse_msl, read_lbr, read_mdlib, save_mdlib, write_mdlib,
    write_msl,
};

#[test]
fn reads_sample_library() {
    let Ok(path) = std::env::var("MD_LBR_SAMPLE") else {
        eprintln!("skipped: MD_LBR_SAMPLE not set");
        return;
    };
    let lib = read_lbr(&path).unwrap();
    assert_eq!(lib.len(), 27);
    assert_eq!(lib[0].name, "Hexatriacontane");
    assert_eq!(format!("{:.3}", lib[0].spectrum.ri), "3599.987");

    let bits = |v: f64| format!("{:016x}", v.to_bits());
    let list = |v: &[String]| v.join(",");
    for c in &lib {
        let s = &c.spectrum;
        let opt = |v: &Option<Vec<f64>>| {
            v.as_ref().map_or("null".into(), |v| {
                list(&v.iter().map(|x| bits(*x)).collect::<Vec<_>>())
            })
        };
        let f32s = |v: &[f32]| list(&v.iter().map(|x| bits(*x as f64)).collect::<Vec<_>>());
        println!(
            "{}|{}|{}|{}|{}|{}|{:?}|{}|{}|{}|{}|{:?}|{}|{}|{}|{:?}|{:?}|{}|{}|{}|{:?}",
            c.name,
            s.scan,
            bits(s.rt_ms),
            bits(s.ri),
            bits(s.bin_width as f64),
            f32s(&s.masses),
            s.values,
            opt(&c.fwhm),
            opt(&c.integrals),
            f32s(&c.quant_ions),
            f32s(&c.auto_quant_ions),
            c.features.iter().collect::<Vec<_>>(),
            bits(c.drt),
            c.time_stamp,
            c.id,
            c.comb_spec
                .iter()
                .map(|(m, v)| (bits(*m as f64), v))
                .collect::<Vec<_>>(),
            c.ion_freq
                .iter()
                .map(|(m, v)| (bits(*m as f64), v))
                .collect::<Vec<_>>(),
            c.spec_count,
            list(&c.ris.iter().map(|x| bits(*x)).collect::<Vec<_>>()),
            list(
                &c.blocks
                    .iter()
                    .map(|b| b.iter().map(|x| format!("{x:02x}")).collect())
                    .collect::<Vec<_>>()
            ),
            c.sources,
        );
    }

    // .lbr → .mdlib → reading back keeps every field (text equality also covers NaN).
    let path = std::env::temp_dir().join(format!("lbr-sample-{}.mdlib", std::process::id()));
    save_mdlib(&path, &lib).unwrap();
    let reopened = read_mdlib(&path);
    std::fs::remove_file(&path).unwrap();
    let reopened = reopened.unwrap();
    assert_eq!(write_mdlib(&reopened), write_mdlib(&lib));
    assert_eq!(reopened, lib);

    let again = parse_msl(&write_msl(&lib)).unwrap();
    assert_eq!(again.len(), lib.len());
    for (a, b) in again.iter().zip(&lib) {
        assert_eq!(
            (&a.name, &a.spectrum.masses, &a.spectrum.values),
            (&b.name, &b.spectrum.masses, &b.spectrum.values)
        );
        assert_eq!(
            (a.spectrum.ri, &a.quant_ions, &a.features),
            (b.spectrum.ri, &b.quant_ions, &b.features)
        );
    }
}

/// Version 5 files written by the old program itself and what it reads back from them (`old.tsv`),
/// made by `docs/reverse-engineering/tools/golden/make-lbr-v5.sh` (never committed). Run with
/// `MD_LBR_V5_DIR=OUT_DIR MD_LBR_SAMPLE=lib_aa.lbr cargo test -p md-core --test lbr_sample`.
#[test]
#[allow(
    clippy::cognitive_complexity,
    reason = "a test of many cases: each assert counts, and Sonar skips tests"
)]
fn reads_old_version_5_libraries() {
    let (Ok(dir), Ok(sample)) = (
        std::env::var("MD_LBR_V5_DIR"),
        std::env::var("MD_LBR_SAMPLE"),
    ) else {
        eprintln!("skipped: MD_LBR_V5_DIR or MD_LBR_SAMPLE not set");
        return;
    };
    let dir = std::path::Path::new(&dir);
    let old = std::fs::read_to_string(dir.join("old.tsv")).unwrap();
    let nums = |f: &[&str]| {
        f.iter()
            .map(|x| x.parse::<f64>().unwrap())
            .collect::<Vec<_>>()
    };
    let fields = |p: &LibraryPeak| {
        vec![
            p.ion as f64,
            p.sharpness,
            p.integral,
            p.quality,
            p.signal_to_noise,
            p.count as f64,
            p.first_scan as f64,
            p.last_scan as f64,
            p.model_parameter,
            p.max_scan as f64,
            p.interpolated_max_scan,
        ]
    };

    // Every file, compound, peak and model value as the old program reads it.
    let mut checked = 0;
    for part in old.split("lib\t").skip(1) {
        let mut lines = part.lines();
        let file = lines.next().unwrap();
        let lib = read_lbr(dir.join(file)).unwrap();
        let mut compounds = lib.iter();
        let mut c = None;
        let mut peaks = Vec::new().into_iter();
        for line in lines {
            let f: Vec<&str> = line.split('\t').collect();
            match f[0] {
                "compound" => {
                    let next = compounds.next().unwrap();
                    assert_eq!(
                        (next.name.as_str(), next.peaks.len()),
                        (f[1], f[2].parse().unwrap())
                    );
                    peaks = next.peaks.clone().into_iter();
                    c = Some(next);
                }
                "peak" => assert_eq!(fields(&peaks.next().unwrap()), nums(&f[1..]), "{file}"),
                "model" => {
                    let m = c.unwrap().model.as_ref().unwrap();
                    let mut want = fields(&m.peak);
                    want.extend([m.noise_factor, m.values.len() as f64]);
                    want.extend(&m.values);
                    assert_eq!(want, nums(&f[1..]), "{file}");
                    checked += 1;
                }
                _ => panic!("unknown line {line}"),
            }
        }
        assert!(compounds.next().is_none(), "{file}");
    }
    assert_eq!(checked, 97);

    // Saved from the version 2 sample: same compounds (saving stamps a new time and id), and each
    // old peak record (f32 ion, f64 sharpness, integral, quality, s/n, u32 count) became a peak.
    let v2 = read_lbr(&sample).unwrap();
    let v5 = read_lbr(dir.join("lib_aa_v5.lbr")).unwrap();
    assert_eq!(v5.len(), v2.len());
    for (a, b) in v5.iter().zip(&v2) {
        assert_eq!(a.peaks.len(), b.blocks.len(), "{}", a.name);
        for (p, block) in a.peaks.iter().zip(&b.blocks) {
            let f64_at = |i: usize| f64::from_le_bytes(block[i..i + 8].try_into().unwrap());
            let old = vec![
                f32::from_le_bytes(block[..4].try_into().unwrap()) as f64,
                f64_at(4),
                f64_at(12),
                f64_at(20),
                f64_at(28),
                u32::from_le_bytes(block[36..40].try_into().unwrap()) as f64,
            ];
            assert_eq!(fields(p)[..6], old, "{}", a.name);
        }
        let mut a = a.clone();
        (a.peaks, a.model, a.blocks) = (Vec::new(), None, b.blocks.clone());
        (a.time_stamp, a.id) = (b.time_stamp, b.id);
        assert_eq!(write_mdlib(&[a]), write_mdlib(std::slice::from_ref(b)));
    }

    let frag = read_lbr(dir.join("fragments_v5.lbr")).unwrap();
    fn formulas(c: &LibraryCompound) -> Vec<(u32, u32, &str)> {
        c.formulas
            .iter()
            .map(|(b, e, f)| (*b, *e, f.as_str()))
            .collect()
    }
    assert_eq!(formulas(&frag[0]), [(5, 6, "N"), (73, 147, "C3H9Si")]);
    assert_eq!(frag[0].mid_sizes, [(5, 6, 0), (73, 147, 2)]);
    assert_eq!(formulas(&frag[1]), [(40, 41, "Ar")]);
    assert_eq!(frag[1].mid_sizes, [(40, 41, 1)]);
    assert!(
        frag[2..]
            .iter()
            .all(|c| c.formulas.is_empty() && c.mid_sizes.is_empty())
    );

    // The old MSL export writes RT in minutes with 4 decimals; its import reads it back the same way.
    let msl = parse_msl(&std::fs::read_to_string(dir.join("lib_aa.msl")).unwrap()).unwrap();
    let imported = read_lbr(dir.join("msl_v5.lbr")).unwrap();
    assert_eq!(msl.len(), frag.len());
    assert_eq!(imported.len(), frag.len());
    for ((m, i), f) in msl.iter().zip(&imported).zip(&frag) {
        assert_eq!((&m.name, &i.name), (&f.name, &f.name));
        let minutes = format!("{:.4}", f.spectrum.rt_ms / 60_000.0);
        assert_eq!(
            format!("{:.4}", m.spectrum.rt_ms / 60_000.0),
            minutes,
            "{}",
            f.name
        );
        assert!(
            (m.spectrum.rt_ms - i.spectrum.rt_ms).abs() < 1e-6,
            "{}",
            f.name
        );
    }
}

/// Compares the three MSL exports with what the old program writes from the same library
/// (docs/reverse-engineering/tools/golden/make-msl-export.sh). Needs `MD_LBR_SAMPLE` and `MD_MSL_GOLDEN`
/// (the OUT_PREFIX of the script); passes without checking anything when one is not set. When
/// `MD_LBR_SAMPLE` is a `.msl` file (the script imports it first) the library is read back from the old
/// program's profile export; that only matches when the input already has a base peak of 1000 (the old
/// import and export both normalise).
#[test]
fn msl_export_matches_the_old_program() {
    use md_core::library::{MslType, write_msl_as};
    let (Ok(lbr), Ok(prefix)) = (
        std::env::var("MD_LBR_SAMPLE"),
        std::env::var("MD_MSL_GOLDEN"),
    ) else {
        eprintln!("skipped: MD_LBR_SAMPLE or MD_MSL_GOLDEN not set");
        return;
    };
    let lib = if lbr.ends_with(".msl") {
        let text = std::fs::read_to_string(format!("{prefix}.profile.msl")).unwrap();
        parse_msl(&text).unwrap()
    } else {
        read_lbr(&lbr).unwrap()
    };
    for (kind, name) in [
        (MslType::Profile, "profile"),
        (MslType::Nominal, "nominal"),
        (MslType::Fwhm, "fwhm"),
    ] {
        let old = std::fs::read_to_string(format!("{prefix}.{name}.msl")).unwrap();
        let new = write_msl_as(&lib, kind);
        for (n, (a, b)) in new.lines().zip(old.lines()).enumerate() {
            assert_eq!(a, b, "{kind:?} line {}", n + 1);
        }
        assert_eq!(new, old, "{kind:?}");
    }
}
