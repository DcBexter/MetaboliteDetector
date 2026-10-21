//! German texts for the error messages of `md-core` (S8, #339). `md-core` has no UI and no language:
//! its errors are English `Display` texts, also for `md-cli` and files. The window translates them where
//! it shows them, with [`text`]. A message is matched against the templates in `errors.txt` (`{}` is a
//! changing part, translated in turn); text before it that ends in `: ` (a file name) stays.
//! shortcut: a message without a template (an operating system text, a new `md-core` wording) stays
//! English; `real_core_messages_are_translated` fails when a covered `md-core` text changes.

use crate::lang;

/// The templates, one `English<TAB>German` line each (`errors.txt`, a data file so the table is not
/// duplicated code); both sides have the same number of `{}`.
const MESSAGES: &str = include_str!("errors.txt");

/// (English `md-core` message, German) in file order: a more specific template comes before a general one.
fn templates() -> impl Iterator<Item = (&'static str, &'static str)> {
    (MESSAGES.lines())
        .filter(|l| !l.starts_with('#'))
        .filter_map(|l| l.split_once('\t'))
}

/// The German text of an `md-core` error message when the window is German, else the message as is.
pub fn text(message: &str) -> String {
    if lang::de() {
        translate(message)
    } else {
        message.to_owned()
    }
}

fn translate(message: &str) -> String {
    // The message, or what follows a `: ` (a file name in front of it stays).
    let starts = std::iter::once(0).chain(message.match_indices(": ").map(|(i, _)| i + 2));
    for at in starts {
        for (en, de) in templates() {
            if let Some(parts) = capture(en, &message[at..]) {
                let parts: Vec<String> = parts.into_iter().map(translate).collect();
                return format!("{}{}", &message[..at], fill(de, &parts));
            }
        }
    }
    message.to_owned()
}

/// The parts of `message` at the `{}` of `template`, if the whole message fits it.
fn capture<'a>(template: &str, message: &'a str) -> Option<Vec<&'a str>> {
    let mut fixed = template.split("{}");
    let mut rest = message.strip_prefix(fixed.next()?)?;
    let fixed: Vec<&str> = fixed.collect();
    let mut parts = Vec::new();
    for (i, lit) in fixed.iter().enumerate() {
        let end = if i + 1 == fixed.len() {
            rest.strip_suffix(lit)?.len()
        } else {
            rest.find(lit)?
        };
        parts.push(&rest[..end]);
        rest = &rest[end + lit.len()..];
    }
    rest.is_empty().then_some(parts)
}

fn fill(template: &str, parts: &[String]) -> String {
    let mut out = String::new();
    for (i, fixed) in template.split("{}").enumerate() {
        if i > 0 {
            out.push_str(&parts[i - 1]);
        }
        out.push_str(fixed);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_core::andi::AndiError;
    use md_core::diskscan::DiskScanError;
    use md_core::fastflight::FastFlightError;
    use md_core::library::{LibraryError, parse_mdlib, parse_msl, parse_targetsearch};
    use md_core::library_add::AddError;
    use md_core::peaks::TooLarge;
    use md_core::ri::{CalibError, RiCalibration};
    use md_core::tof::{TofCalibration, TofError};

    fn german(message: &str) -> String {
        lang::set(lang::Lang::De);
        let de = text(message);
        lang::set(lang::Lang::En);
        de
    }

    #[test]
    fn templates_keep_their_holes() {
        let lines = MESSAGES
            .lines()
            .filter(|l| !l.is_empty() && !l.starts_with('#'));
        assert_eq!(lines.count(), templates().count(), "a line without a tab");
        for (en, de) in templates() {
            assert_eq!(en.matches("{}").count(), de.matches("{}").count(), "{en}");
        }
    }

    #[test]
    fn english_stays_as_it_is() {
        assert_eq!(text("empty formula"), "empty formula");
    }

    #[test]
    fn parts_and_file_names_are_kept() {
        assert_eq!(
            german("C:\\runs\\a.cdf: empty formula"),
            "C:\\runs\\a.cdf: Leere Formel"
        );
        assert_eq!(
            german("MSL line 3: invalid m/z"),
            "MSL-Zeile 3: Ungültiges m/z"
        );
        assert_eq!(
            german("x.cdf: inconsistent ANDI-MS data: the run has no data points"),
            "x.cdf: Widersprüchliche ANDI-MS-Daten: Der Lauf hat keine Datenpunkte"
        );
        assert_eq!(
            german("cannot read x.bin: Disk full"),
            "x.bin kann nicht gelesen werden: Disk full"
        );
        assert_eq!(german("something new"), "something new");
    }

    const ENGLISH: [&str; 11] = [
        "expected `",
        "invalid",
        "missing",
        "twice",
        "unknown",
        "must ",
        "needs",
        "cannot",
        " not ",
        "out of range",
        "but ",
    ];

    fn err<T, E: std::fmt::Display>(result: Result<T, E>) -> String {
        result.err().expect("an error").to_string()
    }

    /// Real texts from `md-core`, so a changed wording there fails here and not silently in the window.
    #[test]
    fn real_core_messages_are_translated() {
        let msl = [
            "NAME: A\nNUM PEAKS: 1\n(1 x)",
            "NAME: A\nNUM PEAKS: 2\n(1 2)",
            "(1 2)",
            "1 2",
            "NAME: A\nNUM PEAKS: 1\n(1 2) x",
            "NAME: A\nNUM PEAKS: 1\n(1",
            "NAME: A\nNUM PEAKS: 1\n(1 2 3)",
            "NAME: A\nRI: x",
            "RI: 1",
            "bad line",
            "NAME: A\nNUM PEAKS: 1\n(1 -2)",
            "NAME: A\nQUANTIFICATION: 1 x",
            "NAME: A\nNUM PEAKS: x",
            "NAME: A\nNUM PEAKS: 1\n71 100 \"C5H11+",
            "NAME: A\nNUM PEAKS: 2\n71 100 72",
            "NAME: A\nNUM PEAKS: 1\n71 x",
        ];
        let mdlib = [
            "end",
            "compound\nend",
            "compound\nname A",
            "compound\nname A\nname B\nend",
            "compound\ncompound\nend",
            "compound\nname A\nscan x\nend",
            "compound\nname A\npeaks 1\nend",
            "compound\nname A\nfeature k\nend",
            "compound\nname A\nfeature k\tv\nfeature k\tv\nend",
            "compound\nname A\nblock 12\nend",
            "compound\nname A\npeak 1 2\nend",
            "compound\nname A\nmodel 1 2\nend",
            "compound\nname A\nformula 1\nend",
            "compound\nname A\nmid_size 1 2\nend",
            "compound\nname A\nwhat 1\nend",
            "compound\nname A\\q\nend",
            "compound\nname A\\\nend",
        ]
        .map(|body| format!("MetaboliteDetector library 2\n{body}"));
        let settings = md_core::settings::GcmsSettings::default();
        let mut bytes = settings.to_bytes();
        let truncated = bytes[..bytes.len() - 1].to_vec();
        let last = bytes.len() - 1;
        bytes[last] ^= 0x40;
        let table = md_core::batch::SampleTable::parse("A\tB\nx\t1\nx\t2").unwrap();
        let invalid_settings = |sf_fwhm| md_core::settings::GcmsSettings {
            sf_fwhm,
            ..Default::default()
        };
        let io = || std::io::Error::other("x");
        let mut errors: Vec<String> = [
            // Typed errors.
            AndiError::MissingVariable("mass_values").to_string(),
            AndiError::Inconsistent("one TIC per scan needed".into()).to_string(),
            DiskScanError::Magic(5).to_string(),
            DiskScanError::Corrupt("scan 3 listed twice".into()).to_string(),
            DiskScanError::Io("a.idx".into(), io()).to_string(),
            FastFlightError::NoCalibration("a".into()).to_string(),
            FastFlightError::Corrupt("scan at 5 ms: end of file inside a record".into())
                .to_string(),
            FastFlightError::Calibration(TofError::NotEnoughPairs).to_string(),
            TofError::Singular.to_string(),
            LibraryError::Io(io()).to_string(),
            LibraryError::Corrupt {
                offset: 3,
                reason: "truncated, 4 more bytes expected".into(),
            }
            .to_string(),
            AddError::NoEntry(2).to_string(),
            AddError::NameTaken("A".into()).to_string(),
            AddError::NoName.to_string(),
            CalibError::Io(io()).to_string(),
            TooLarge {
                masses: 1,
                scans: 2,
            }
            .to_string(),
            // Texts md-core builds itself.
            "could not save a.cmp: disk full".into(),
            "not saved: a.cmp is an input of this batch".into(),
            "not saved: another run is also saved as a.cmp".into(),
            "Warning: No fragment (FRAG) defined for compound A".into(),
            "Error: Can not parse \"x\". SYNTAX(e.g): (174/8)(261/5)".into(),
            // Parsers.
            err(parse_targetsearch("x")),
            err(parse_targetsearch(
                "Name\tRI\tWin_1\tSEL_MASS\tSPECTRUM\nA\t1",
            )),
            err(parse_mdlib("nonsense")),
            err(parse_mdlib("MetaboliteDetector library 99")),
            err(RiCalibration::new([(1.0, 1.0), (1.0 + 1e-12, 2.0)])),
            err(RiCalibration::new([(f64::NAN, 1.0)])),
            err(settings.from_bytes(&[])),
            err(settings.from_bytes(&[0; 200])),
            err(settings.from_bytes(&[0x1f, 0x8b, 0x08])),
            err(settings.from_bytes(&truncated)),
            err(settings.from_bytes(&bytes)),
            err(invalid_settings(2.0).validate()),
            err(invalid_settings(f64::NAN).validate()),
            err(md_core::batch::SampleTable::parse("")),
            err(table.column("C")),
            err(table.groups(0, 0)),
            err(table.groups(0, 1)),
        ]
        .into();
        errors.extend(msl.map(|text| err(parse_msl(text))));
        errors.extend(mdlib.map(|text| err(parse_mdlib(&text))));
        errors.extend(
            ["1 2 3", "1 x", "1 5", "1 5\n2 5"].map(|text| err(RiCalibration::parse_calib(text))),
        );
        errors.extend(
            [
                "",
                "Polynomial order:",
                "Polynomial order: 99",
                "Polynomial order: 1 Parameter: 0 1",
            ]
            .map(|text| err(TofCalibration::parse(text))),
        );
        errors.extend(
            ["(C", "C)", "C-", "", "C99999999999"]
                .map(|text| err(md_core::isotope::parse_formula(text))),
        );
        let untranslated: Vec<_> = errors
            .into_iter()
            .filter(|e| {
                let de = german(e);
                // A tail left in English: a text nothing matched, or a template missing for its reason.
                de == *e || ENGLISH.iter().any(|w| de.contains(w))
            })
            .collect();
        assert!(
            untranslated.is_empty(),
            "no template:\n{}",
            untranslated.join("\n")
        );
    }
}
