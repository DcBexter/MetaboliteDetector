//! Reads the sample library `lib_aa.lbr` from the old 3.3 package (never committed).
//! Run with `MD_LBR_SAMPLE=path/to/lib_aa.lbr cargo test -p md-core --test lbr_sample -- --nocapture`;
//! the printed dump matches `docs/reverse-engineering/tools/dump_lbr.py` field by field.
//! Passes without checking anything when the variable is not set.

use md_core::library::{parse_msl, read_lbr, read_mdlib, save_mdlib, write_mdlib, write_msl};

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

/// Compares the three MSL exports with what the old program writes from the same library
/// (docs/reverse-engineering/tools/golden/make-msl-export.sh). Needs `MD_LBR_SAMPLE` and `MD_MSL_GOLDEN`
/// (the OUT_PREFIX of the script); passes without checking anything when one is not set.
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
    let lib = read_lbr(&lbr).unwrap();
    // The script names the files after the old ExportType numbers: 0 profile, 1 nominal, 2 FWHM.
    for (kind, name) in [
        (MslType::Profile, "nominal"),
        (MslType::Nominal, "type1"),
        (MslType::Fwhm, "type2"),
    ] {
        let old = std::fs::read_to_string(format!("{prefix}.{name}.msl")).unwrap();
        let new = write_msl_as(&lib, kind);
        for (n, (a, b)) in new.lines().zip(old.lines()).enumerate() {
            assert_eq!(a, b, "{kind:?} line {}", n + 1);
        }
        assert_eq!(new, old, "{kind:?}");
    }
}
