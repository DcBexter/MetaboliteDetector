//! German texts for the error messages of `md-core` (S8, #339). `md-core` has no UI and no language:
//! its errors are English `Display` texts, also for `md-cli` and files. The window translates them where
//! it shows them, with [`text`]. A message is matched against the templates in [`MESSAGES`] (`{}` is a
//! changing part, translated in turn); text before it that ends in `: ` (a file name) stays.
//! shortcut: a message without a template (an operating system text, a new `md-core` wording) stays
//! English; `real_core_messages_are_translated` fails when a covered `md-core` text changes.

use crate::lang;

/// (English `md-core` message, German); both have the same number of `{}`. A more specific
/// template comes before a general one that would also match.
const MESSAGES: &[(&str, &str)] = &[
    // andi.rs
    (
        "cannot read NetCDF file: {}",
        "NetCDF-Datei kann nicht gelesen werden: {}",
    ),
    (
        "cannot write NetCDF file: {}",
        "NetCDF-Datei kann nicht geschrieben werden: {}",
    ),
    (
        "missing ANDI-MS variable `{}`",
        "ANDI-MS-Variable `{}` fehlt",
    ),
    (
        "inconsistent ANDI-MS data: {}",
        "Widersprüchliche ANDI-MS-Daten: {}",
    ),
    (
        "{} points do not fit ANDI-MS",
        "{} Punkte passen nicht in ANDI-MS",
    ),
    ("one TIC per scan needed", "Je Scan wird ein TIC benötigt"),
    (
        "the run has no data points",
        "Der Lauf hat keine Datenpunkte",
    ),
    (
        "{} scan times but {} scan indices and {} point counts",
        "{} Scanzeiten, aber {} Scan-Indizes und {} Punktzahlen",
    ),
    (
        "scan times are not ascending",
        "Die Scanzeiten sind nicht aufsteigend",
    ),
    (
        "{} masses but {} intensities",
        "{} Massen, aber {} Intensitäten",
    ),
    (
        "scan {} points to values {}..{} of {}",
        "Scan {} verweist auf die Werte {}..{} von {}",
    ),
    // diskscan.rs
    (
        "not a MetaboliteDetector index file (magic {})",
        "Keine Indexdatei von MetaboliteDetector (Kennung {})",
    ),
    (
        "corrupt index or data file: {}",
        "Beschädigte Index- oder Datendatei: {}",
    ),
    ("scan {} listed twice", "Scan {} ist doppelt aufgeführt"),
    (
        "index lists {} points, data file has {} bytes",
        "Der Index nennt {} Punkte, die Datendatei hat {} Bytes",
    ),
    (
        "record of scan {} does not match the index",
        "Der Datensatz von Scan {} passt nicht zum Index",
    ),
    (
        "unsupported type sizes {}",
        "Nicht unterstützte Typgrößen {}",
    ),
    ("TOF polynomial with {} terms", "TOF-Polynom mit {} Termen"),
    (
        "{} scans do not fit into the index",
        "{} Scans passen nicht in den Index",
    ),
    (
        "scan {} has baseline mean {}",
        "Scan {} hat den Grundlinienmittelwert {}",
    ),
    (
        "offset {} is past the end of the data file",
        "Der Offset {} liegt hinter dem Ende der Datendatei",
    ),
    (
        "file ends at byte {}, needed {} more",
        "Die Datei endet bei Byte {}, es fehlen {} weitere",
    ),
    // fastflight.rs
    (
        "corrupt FastFlight data: {}",
        "Beschädigte FastFlight-Daten: {}",
    ),
    (
        "no .mzcalib TOF calibration file next to {}",
        "Keine .mzcalib-Datei mit der TOF-Kalibrierung neben {}",
    ),
    ("scan at {} ms: {}", "Scan bei {} ms: {}"),
    (
        "TIC mismatch! TICs in the file .FFT and .FFC do not match.",
        "TIC-Abweichung: Die TIC-Werte der Dateien .FFT und .FFC stimmen nicht überein.",
    ),
    (
        "the index points several scans at the same data",
        "Der Index verweist mehrere Scans auf dieselben Daten",
    ),
    ("{}: invalid time factor", "{}: ungültiger Zeitfaktor"),
    (
        "end of file inside a record",
        "Dateiende mitten in einem Datensatz",
    ),
    (
        "expected new frame header",
        "Kopf eines neuen Frames erwartet",
    ),
    (
        "frame shorter than its header",
        "Frame kürzer als sein Kopf",
    ),
    // tof.rs
    (
        "invalid TOF calibration file: {}",
        "Ungültige TOF-Kalibrierdatei: {}",
    ),
    (
        "Not enough values for polynomial regression.",
        "Zu wenige Werte für die Polynomregression.",
    ),
    (
        "the regression has no unique solution (the calibration ions are too similar)",
        "Die Regression hat keine eindeutige Lösung (die Kalibrierionen sind zu ähnlich)",
    ),
    ("missing header", "Kopfzeile fehlt"),
    ("missing polynomial order", "Ordnung des Polynoms fehlt"),
    (
        "polynomial order {} is outside 0..=64",
        "Die Ordnung {} des Polynoms liegt außerhalb von 0..=64",
    ),
    (
        "parameter {} is missing or invalid",
        "Parameter {} fehlt oder ist ungültig",
    ),
    // The general file error of diskscan.rs, fastflight.rs and tof.rs, after the specific "cannot read ...".
    (
        "cannot read library file: {}",
        "Bibliotheksdatei kann nicht gelesen werden: {}",
    ),
    (
        "cannot read calibration file: {}",
        "Kalibrierdatei kann nicht gelesen werden: {}",
    ),
    ("cannot read {}: {}", "{} kann nicht gelesen werden: {}"),
    // library.rs, cmp.rs
    ("not supported: {}", "Nicht unterstützt: {}"),
    (
        "corrupt library at byte {}: {}",
        "Beschädigte Bibliothek bei Byte {}: {}",
    ),
    ("MSL line {}: {}", "MSL-Zeile {}: {}"),
    ("library line {}: {}", "Bibliothekszeile {}: {}"),
    ("library version {}", "Bibliotheksversion {}"),
    (
        "libraries larger than {} bytes uncompressed",
        "Bibliotheken über {} Bytes entpackt",
    ),
    (
        "libraries larger than {} bytes",
        "Bibliotheken über {} Bytes",
    ),
    (
        "compound files larger than {} bytes",
        "Substanzdateien über {} Bytes",
    ),
    ("compound file version {}", "Version {} der Substanzdatei"),
    (
        "library text format version {}",
        "Version {} des Bibliothekstexts",
    ),
    (
        "version 2 and 3 compounds with peak records (layout unknown)",
        "Substanzen der Version 2 und 3 mit Peak-Einträgen (Aufbau unbekannt)",
    ),
    (
        "no NIST record with a peak above m/z 69",
        "Kein NIST-Eintrag mit einem Peak über m/z 69",
    ),
    (
        "{} bytes after the last compound",
        "{} Bytes nach der letzten Substanz",
    ),
    (
        "spectrum check value mismatch",
        "Prüfwert des Spektrums stimmt nicht",
    ),
    (
        "spectrum has different numbers of masses and values",
        "Das Spektrum hat verschieden viele Massen und Werte",
    ),
    (
        "compound check value mismatch",
        "Prüfwert der Substanz stimmt nicht",
    ),
    (
        "model peak check value mismatch",
        "Prüfwert des Modellpeaks stimmt nicht",
    ),
    (
        "truncated, {} more bytes expected",
        "Abgeschnitten, es fehlen {} Bytes",
    ),
    (
        "count {} exceeds the remaining data",
        "Die Anzahl {} übersteigt die restlichen Daten",
    ),
    (
        "{}: component at scan {} ({} ms) is not from this run",
        "{}: Die Komponente bei Scan {} ({} ms) gehört nicht zu diesem Lauf",
    ),
    (
        "{}: scan number {} is not plausible",
        "{}: Die Scannummer {} ist nicht plausibel",
    ),
    (
        "{}: component at scan {} does not fit a run of {} scans",
        "{}: Die Komponente bei Scan {} passt nicht zu einem Lauf mit {} Scans",
    ),
    // library.rs: MSL
    ("expected the header `{}`", "Kopfzeile `{}` erwartet"),
    (
        "expected 5 tab-separated columns",
        "5 durch Tabulator getrennte Spalten erwartet",
    ),
    (
        "peak list before NAME: and NUM PEAKS:",
        "Peakliste vor NAME: und NUM PEAKS:",
    ),
    ("missing `)`", "`)` fehlt"),
    ("text between peak pairs", "Text zwischen den Peakpaaren"),
    ("expected `(m/z intensity)`", "`(m/z Intensität)` erwartet"),
    (
        "unclosed `\"` annotation",
        "Anmerkung mit `\"` nicht geschlossen",
    ),
    (
        "expected `m/z intensity` pairs",
        "Paare `m/z Intensität` erwartet",
    ),
    ("invalid m/z", "Ungültiges m/z"),
    ("invalid intensity", "Ungültige Intensität"),
    (
        "m/z or intensity out of range",
        "m/z oder Intensität außerhalb des Bereichs",
    ),
    (
        "NUM PEAKS Some({}) but {} pairs found",
        "NUM PEAKS {}, aber {} Paar(e) gefunden",
    ),
    (
        "NUM PEAKS None but {} pairs found",
        "NUM PEAKS fehlt, aber {} Paar(e) gefunden",
    ),
    (
        "expected `KEY: value` or `(m/z intensity)`",
        "`SCHLÜSSEL: Wert` oder `(m/z Intensität)` erwartet",
    ),
    (
        "expected NAME: before other keys",
        "NAME: vor den anderen Schlüsseln erwartet",
    ),
    ("invalid number", "Ungültige Zahl"),
    ("invalid peak count", "Ungültige Peakanzahl"),
    ("invalid ion", "Ungültiges Ion"),
    // library.rs: .mdlib
    ("compound without `name`", "Substanz ohne `name`"),
    ("expected `compound`", "`compound` erwartet"),
    ("`compound` before `end`", "`compound` vor `end`"),
    ("`{}` given twice", "`{}` doppelt angegeben"),
    ("compound without `end`", "Substanz ohne `end`"),
    ("invalid number `{}`", "Ungültige Zahl `{}`"),
    (
        "expected `m/z:value`, got `{}`",
        "`m/z:Wert` erwartet, gefunden `{}`",
    ),
    ("expected `KEY<tab>value`", "`SCHLÜSSEL<Tab>Wert` erwartet"),
    ("feature `{}` given twice", "Merkmal `{}` doppelt angegeben"),
    (
        "block must be {} hex digits",
        "Ein Block muss {} Hexziffern haben",
    ),
    ("peak must have 11 fields", "Ein Peak braucht 11 Felder"),
    (
        "model must have at least 12 fields",
        "Ein Modell braucht mindestens 12 Felder",
    ),
    (
        "expected `begin end formula`",
        "`Anfang Ende Formel` erwartet",
    ),
    ("expected `begin end size`", "`Anfang Ende Größe` erwartet"),
    ("unknown key `{}`", "Unbekannter Schlüssel `{}`"),
    ("unknown escape `\\{}`", "Unbekannte Escape-Folge `\\{}`"),
    ("`\\` at the end of a value", "`\\` am Ende eines Werts"),
    // library_add.rs
    (
        "the library has no entry {}",
        "Die Bibliothek hat keinen Eintrag {}",
    ),
    (
        "the library already has an entry named {}",
        "Die Bibliothek hat schon einen Eintrag namens {}",
    ),
    (
        "the compound needs a name",
        "Die Substanz braucht einen Namen",
    ),
    // ri.rs
    (
        "calibration file line {}: {}",
        "Kalibrierdatei, Zeile {}: {}",
    ),
    ("invalid calibration: {}", "Ungültige Kalibrierung: {}"),
    ("not a finite number", "Keine endliche Zahl"),
    (
        "{} point(s), at least 2 needed",
        "{} Punkt(e), mindestens 2 nötig",
    ),
    (
        "RI {} at {} s is not above RI {} at {} s",
        "RI {} bei {} s liegt nicht über RI {} bei {} s",
    ),
    (
        "slope between {} s and {} s is outside 1e-9 to 1e9 RI per second",
        "Die Steigung zwischen {} s und {} s liegt außerhalb von 1e-9 bis 1e9 RI je Sekunde",
    ),
    ("file larger than {} bytes", "Datei größer als {} Bytes"),
    ("not a text file", "Keine Textdatei"),
    (
        "binary calibration is truncated",
        "Die binäre Kalibrierung ist abgeschnitten",
    ),
    (
        "retention time shift is not finite",
        "Die Verschiebung der Retentionszeit ist keine endliche Zahl",
    ),
    ("not a number", "Keine Zahl"),
    (
        "expected two numbers: retention time (ms) and RI",
        "Zwei Zahlen erwartet: Retentionszeit (ms) und RI",
    ),
    // peaks.rs
    (
        "ion index of {} masses x {} scans needs more than {} GiB",
        "Der Ionenindex aus {} Massen x {} Scans braucht mehr als {} GiB",
    ),
    // settings.rs
    ("{}: value out of range", "{}: Wert außerhalb des Bereichs"),
    (
        "{}: must be between {} and {}",
        "{}: muss zwischen {} und {} liegen",
    ),
    // batch.rs, batch/rescan.rs
    ("analysis crashed", "Die Analyse ist abgestürzt"),
    (
        "could not save {}: {}",
        "{} konnte nicht gespeichert werden: {}",
    ),
    (
        "not saved: {} is an input of this batch",
        "Nicht gespeichert: {} ist eine Eingabe dieses Batches",
    ),
    (
        "not saved: another run is also saved as {}",
        "Nicht gespeichert: ein anderer Lauf wird auch als {} gespeichert",
    ),
    // settings.rs: from_bytes
    ("settings file too large", "Einstellungsdatei zu groß"),
    (
        "not a readable gzip file: {}",
        "Keine lesbare gzip-Datei: {}",
    ),
    (
        "settings file too short ({} bytes)",
        "Einstellungsdatei zu kurz ({} Bytes)",
    ),
    (
        "settings file has {} bytes, not 191 + 8 x {}",
        "Die Einstellungsdatei hat {} Bytes statt 191 + 8 x {}",
    ),
    (
        "settings file checksum does not match",
        "Die Prüfsumme der Einstellungsdatei stimmt nicht",
    ),
    // batch/samples.rs
    ("the table is empty", "Die Tabelle ist leer"),
    (
        "no column {}; the table has {}",
        "Keine Spalte {}; die Tabelle hat {}",
    ),
    (
        "choose different columns for the files and the groups",
        "Für die Dateien und die Gruppen müssen verschiedene Spalten gewählt werden",
    ),
    (
        "{} is in groups {} and {}",
        "{} liegt in den Gruppen {} und {}",
    ),
    // isotope.rs, mid/series.rs
    ("count too large: {}", "Anzahl zu groß: {}"),
    ("atom count too large", "Atomanzahl zu groß"),
    ("unmatched )", "`)` ohne `(`"),
    ("unmatched (", "`(` ohne `)`"),
    (
        "unexpected {} at position {}",
        "Unerwartetes Zeichen {} an Position {}",
    ),
    ("empty formula", "Leere Formel"),
    ("unknown element {}", "Unbekanntes Element {}"),
    ("no tracer", "Kein Tracer"),
    ("no sum formula", "Keine Summenformel"),
    ("{} ({}-{}): {}", "{} ({}-{}): {}"),
    (
        "Error: Can not parse \"{}\". SYNTAX(e.g): (174/8)(261/5)",
        "Fehler: \"{}\" ist nicht lesbar. Syntax, zum Beispiel: (174/8)(261/5)",
    ),
    (
        "Warning: No fragment (FRAG) defined for compound {}",
        "Warnung: Für die Substanz {} ist kein Fragment (FRAG) definiert",
    ),
    // stats.rs
    (
        "PCA needs at least two runs",
        "Die PCA braucht mindestens zwei Läufe",
    ),
    (
        "PCA needs at least one substance with a value in every run and some spread",
        "Die PCA braucht mindestens eine Substanz mit einem Wert in jedem Lauf und etwas Streuung",
    ),
    // The .mdlib header, after the `expected ...` texts above that this would also match.
    ("expected `{}`", "`{}` erwartet"),
];

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
        for (en, de) in MESSAGES {
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
        for (en, de) in MESSAGES {
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
