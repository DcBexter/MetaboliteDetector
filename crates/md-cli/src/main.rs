//! Headless batch quantification for big studies: finds the library substances in many runs and writes
//! the substance × run table of the app's batch window as CSV. Depends on `md-core` only.

use std::ffi::OsString;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::time::Instant;

use md_core::batch::{
    RunResult, SampleTable, ShiftReference, TableSettings, complete_scan, ions_path, scan_runs,
    settings_path, table,
};
use md_core::library;
use md_core::ri::RiCalibration;
use md_core::search::{PreparedLibrary, SearchSettings};
use md_core::settings::GcmsSettings;

const USAGE: &str = "\
Usage: md-cli --library FILE [options] RUN...

Finds the substances of a library (.mdlib, .lbr, .msl, TargetSearch .txt, NIST .msp or .sdf) in many runs (.cdf, .nc, old
.bin or saved components .cmp, which are not deconvoluted again and cannot be rescanned) and writes the
substance x run table as CSV. Exit code 2: the table is written, but some runs could not be read
(or saved).

Options:
  -l, --library FILE          library (required)
  -s, --settings FILE         old MetaboliteDetector settings file (default: the old defaults, but the
                              deconvolution width follows each run's scan rate: 0.5 s, at least 1 scan)
  -c, --calibration FILE      RI calibration (.calib); without it the search uses the spectrum only
  -o, --output FILE           CSV file (default: standard output); with groups (--replicates, --sample-table) the group of each
                              run goes to FILE.groups.csv and the statistics of the groups (mean, rel.
                              std. error, hits, p-value against a group BLANK, ANOVA) to FILE.groupstats.csv,
                              the settings used to FILE.settings.csv
  -j, --parallel N            runs in memory at once (default 4)
      --replicates REGEX      replicate group = first capture group of the file name
      --sample-table FILE     replicate groups from a CSV/TSV table (first line = column names); runs it
                              lacks are grouped by --replicates
      --sample-file-column NAME column with the run file names (default: CDF_FILE, else the first)
      --sample-group-column NAME column with the group (default: TIME_POINT, else the last)
      --min-reproducibility F keep a substance only if a group quantifies it in this share of its runs
      --min-ions N            quant ions that must be found for a value
      --ions N                quant ions per substance (default: QUANT_ION_NUMBER of the settings)
      --min-dist F            min. distance of quant ions (default: QUANT_MIN_DIST_IONS)
      --quality-cutoff F      quant ion quality cut-off (default: QUANT_QUALITY_CUTOFF)
      --exclude-ions LIST     m/z ranges not picked as quant ions, e.g. 73,147-149 (default: QUANT_EXCL_IONS)
      --min-sn F              an ion counts only with at least this signal-to-noise ratio
      --max-discrepancy F     a value counts only if the mean peak quality of its ions is at most F
      --ion-integrals         also write the integral of every quant ion to FILE.ions.csv (needs -o)
      --excel-de              write the CSV files for a spreadsheet with German locale: `;` and decimal
                              commas (needs -o)
      --untargeted            also tabulate components without library hit (pooled over all runs)
      --rescan                look for quant ions a run lacks in its ion traces
      --baseline-correction   subtract each scan's baseline of every run (the old wizard's switch; old type 1 .bin/.idx caches only)
      --save-cmp DIR          also save the components of each deconvoluted run as DIR/RUN.cmp
      --shift-reference NAME  shift the calibration per run onto this library substance (needs -c)
      --shift-cutoff F        min. score of the reference in a run (default: the library search cutoff)
  -h, --help                  this text";

struct Args {
    library: PathBuf,
    settings: Option<PathBuf>,
    calibration: Option<PathBuf>,
    output: Option<PathBuf>,
    parallel: usize,
    replicates: String,
    sample_table: Option<PathBuf>,
    sample_file_column: Option<String>,
    sample_group_column: Option<String>,
    min_reproducibility: Option<f64>,
    min_ions: Option<usize>,
    ions: Option<usize>,
    min_dist: Option<f32>,
    quality_cutoff: Option<f32>,
    exclude_ions: Option<Vec<(f32, f32)>>,
    min_sn: Option<f32>,
    max_discrepancy: Option<f32>,
    ion_integrals: bool,
    excel_de: bool,
    untargeted: bool,
    rescan: bool,
    baseline_correction: bool,
    save_cmp: Option<PathBuf>,
    shift_reference: Option<String>,
    shift_cutoff: Option<f64>,
    runs: Vec<PathBuf>,
}

/// `73,147-149` → inclusive ranges; a single mass is a range of its own.
fn ranges(text: &str) -> Result<Vec<(f32, f32)>, String> {
    let mass = |t: &str| {
        t.trim()
            .parse::<f32>()
            .ok()
            .filter(|m| (0.0..=100_000.0).contains(m))
            .ok_or_else(|| format!("not a mass: {t}"))
    };
    text.split(',')
        .filter(|r| !r.trim().is_empty())
        .map(|r| match r.split_once('-') {
            Some((lo, hi)) => Ok((mass(lo)?, mass(hi)?)),
            None => mass(r).map(|m| (m, m)),
        })
        .collect()
}

/// `Ok(None)` for `--help`.
fn parse(args: impl IntoIterator<Item = OsString>) -> Result<Option<Args>, String> {
    let mut args = args.into_iter();
    let mut library = None;
    let mut a = Args {
        library: PathBuf::new(),
        settings: None,
        calibration: None,
        output: None,
        parallel: 4,
        replicates: String::new(),
        sample_table: None,
        sample_file_column: None,
        sample_group_column: None,
        min_reproducibility: None,
        min_ions: None,
        ions: None,
        min_dist: None,
        quality_cutoff: None,
        exclude_ions: None,
        min_sn: None,
        max_discrepancy: None,
        ion_integrals: false,
        excel_de: false,
        untargeted: false,
        rescan: false,
        baseline_correction: false,
        save_cmp: None,
        shift_reference: None,
        shift_cutoff: None,
        runs: Vec::new(),
    };
    let mut only_runs = false;
    while let Some(arg) = args.next() {
        let name = arg.to_string_lossy().into_owned();
        if only_runs || !name.starts_with('-') || name == "-" {
            a.runs.push(arg.into());
            continue;
        }
        let mut value = || args.next().ok_or_else(|| format!("{name}: value missing"));
        let mut text = || {
            value()?
                .into_string()
                .map_err(|_| format!("{name}: not valid text"))
        };
        // The ranges of the app's batch window.
        let mut num = |min: f64, max: f64| {
            let v = text()?;
            v.parse::<f64>()
                .ok()
                .filter(|n| (min..=max).contains(n))
                .ok_or_else(|| format!("{name}: not a number from {min} to {max}: {v}"))
        };
        match name.as_str() {
            "--" => only_runs = true,
            "-h" | "--help" => return Ok(None),
            "-l" | "--library" => library = Some(PathBuf::from(value()?)),
            "-s" | "--settings" => a.settings = Some(value()?.into()),
            "-c" | "--calibration" => a.calibration = Some(value()?.into()),
            "-o" | "--output" => a.output = Some(value()?.into()),
            "-j" | "--parallel" => a.parallel = num(1.0, 64.0)? as usize,
            "--replicates" => a.replicates = text()?,
            "--sample-table" => a.sample_table = Some(value()?.into()),
            "--sample-file-column" => a.sample_file_column = Some(text()?),
            "--sample-group-column" => a.sample_group_column = Some(text()?),
            "--min-reproducibility" => a.min_reproducibility = Some(num(0.0, 1.0)?),
            "--min-ions" => a.min_ions = Some(num(1.0, 20.0)? as usize),
            "--ions" => a.ions = Some(num(1.0, 20.0)? as usize),
            "--min-dist" => a.min_dist = Some(num(0.0, 100.0)? as f32),
            "--quality-cutoff" => a.quality_cutoff = Some(num(0.0, 100.0)? as f32),
            "--exclude-ions" => a.exclude_ions = Some(ranges(&text()?)?),
            "--min-sn" => a.min_sn = Some(num(0.0, 1e6)? as f32),
            "--max-discrepancy" => a.max_discrepancy = Some(num(0.0, 1e6)? as f32),
            "--ion-integrals" => a.ion_integrals = true,
            "--excel-de" => a.excel_de = true,
            "--untargeted" => a.untargeted = true,
            "--rescan" => a.rescan = true,
            "--baseline-correction" => a.baseline_correction = true,
            "--save-cmp" => a.save_cmp = Some(value()?.into()),
            "--shift-reference" => a.shift_reference = Some(text()?),
            "--shift-cutoff" => a.shift_cutoff = Some(num(0.0, 1.0)?),
            _ => return Err(format!("unknown option {name}")),
        }
    }
    a.library = library.ok_or("--library is required")?;
    if a.runs.is_empty() {
        return Err("no runs given".into());
    }
    if a.sample_table.is_none()
        && (a.sample_file_column.is_some() || a.sample_group_column.is_some())
    {
        return Err("--sample-file-column and --sample-group-column need --sample-table".into());
    }
    if a.ion_integrals && a.output.is_none() {
        return Err("--ion-integrals needs -o".into());
    }
    if a.excel_de && a.output.is_none() {
        return Err("--excel-de needs -o".into());
    }
    Ok(Some(a))
}

/// `Ok(false)` if some runs could not be read.
fn run(a: Args) -> Result<bool, String> {
    let settings = match &a.settings {
        Some(p) => {
            let bytes = std::fs::read(p).map_err(|e| format!("{}: {e}", p.display()))?;
            GcmsSettings::default()
                .from_bytes(&bytes)
                .map_err(|e| format!("{}: {e}", p.display()))?
        }
        None => GcmsSettings::default(),
    };
    let calibration = a
        .calibration
        .as_ref()
        .map(|p| RiCalibration::read_calib(p).map_err(|e| format!("{}: {e}", p.display())))
        .transpose()?;
    let compounds =
        library::open(&a.library).map_err(|e| format!("{}: {e}", a.library.display()))?;
    let prepared = PreparedLibrary::new(&compounds, &settings.search())
        .with_fragment_filter(library::is_nist(&a.library));

    let mut ts = TableSettings {
        quant: settings.quant(),
        replicates: a.replicates,
        ..Default::default()
    };
    let mut table_source = None;
    if let Some(p) = &a.sample_table {
        let text = std::fs::read_to_string(p).map_err(|e| format!("{}: {e}", p.display()))?;
        let t = SampleTable::parse(&text).map_err(|e| format!("{}: {e}", p.display()))?;
        let col = |name: &Option<String>, default| match name {
            Some(n) => t.column(n),
            None => Ok(default),
        };
        let (file, group) = t.default_columns();
        let (file, group) = (
            col(&a.sample_file_column, file)?,
            col(&a.sample_group_column, group)?,
        );
        ts.sample_groups = t
            .groups(file, group)
            .map_err(|e| format!("{}: {e}", p.display()))?;
        table_source = Some(format!(
            "{} (files {}, groups {})",
            p.file_name().unwrap_or_default().to_string_lossy(),
            t.columns[file],
            t.columns[group]
        ));
    }
    ts.quant.ion_number = a.ions.unwrap_or(ts.quant.ion_number);
    ts.quant.min_dist = a.min_dist.unwrap_or(ts.quant.min_dist);
    ts.quant.quality_cutoff = a.quality_cutoff.unwrap_or(ts.quant.quality_cutoff);
    ts.min_reproducibility = a.min_reproducibility.unwrap_or(ts.min_reproducibility);
    ts.min_ions = a.min_ions.unwrap_or(ts.min_ions);
    ts.quant.excluded = a.exclude_ions.unwrap_or(ts.quant.excluded);
    ts.min_sn = a.min_sn.unwrap_or(ts.min_sn);
    ts.max_discrepancy = a.max_discrepancy.unwrap_or(ts.max_discrepancy);
    let mut scan = settings.scan(calibration);
    scan.untargeted = a.untargeted;
    scan.baseline_correction = a.baseline_correction;
    if let Some(dir) = a.save_cmp {
        // Before the scan: a wrong folder must not cost a whole scan.
        std::fs::create_dir_all(&dir).map_err(|e| format!("{}: {e}", dir.display()))?;
        scan.save_cmp = Some(md_core::batch::CmpTarget::Dir(dir));
    }
    if let Some(name) = &a.shift_reference {
        if scan.calibration.is_none() {
            return Err("--shift-reference needs a calibration (-c)".into());
        }
        let cutoff = a.shift_cutoff.unwrap_or(SearchSettings::default().cutoff);
        scan.shift = Some(
            ShiftReference::find(&compounds, name, cutoff)
                .ok_or_else(|| format!("--shift-reference: no substance with RI named {name}"))?,
        );
    }
    let rescan = a.rescan;
    // How the table was made, for the settings file next to the output.
    // File names only: the file may be shared and an absolute path would carry the user's directories.
    let file = |p: &Path| {
        p.file_name()
            .map_or(String::new(), |n| n.to_string_lossy().into())
    };
    let named = |p: &Option<PathBuf>| p.as_deref().map_or("none".into(), file);
    let mut how = vec![
        ("LIBRARY", file(&a.library)),
        ("SETTINGS_FILE", named(&a.settings)),
        ("CALIBRATION", named(&a.calibration)),
        (
            "SHIFT_REFERENCE",
            a.shift_reference.clone().unwrap_or_else(|| "none".into()),
        ),
        ("UNTARGETED", a.untargeted.to_string()),
        ("BASELINE_CORRECTION", a.baseline_correction.to_string()),
        ("RESCAN", rescan.to_string()),
        ("RUNS_GIVEN", a.runs.len().to_string()),
    ];

    // Open the output first: a wrong path must not cost a whole scan.
    let mut out: Box<dyn Write> = match &a.output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p).map_err(|e| format!("{}: {e}", p.display()))?,
        )),
        None => Box::new(std::io::stdout().lock()),
    };
    let started = Instant::now();
    let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
    let (mut runs, mut errors) =
        scan_runs(&a.runs, &prepared, &scan, a.parallel, &done, &cancel).ok_or("cancelled")?;
    let more = complete_scan(
        &mut runs, &compounds, &scan, &ts.quant, rescan, a.parallel, &cancel,
    )
    .ok_or("cancelled")?;
    errors.extend(more.errors);
    let table =
        table(&runs, &compounds, &more.extra, &ts).map_err(|e| format!("--replicates: {e}"))?;

    let paths: Vec<_> = runs.iter().map(|r| r.path.clone()).collect();
    how.extend([
        (
            "GROUP_SOURCE",
            ts.group_source(table_source.as_deref(), &paths),
        ),
        ("RT_ONLY", prepared.is_rt_only().to_string()),
        ("SEARCH_CUTOFF", settings.search().cutoff.to_string()),
        ("RUNS_SCANNED", runs.len().to_string()),
        ("UNKNOWN_SUBSTANCES", more.extra.len().to_string()),
        ("IONS_ADDED_BY_RESCAN", more.rescanned.to_string()),
    ]);
    let io = |e: std::io::Error| format!("output: {e}");
    table.write_csv(&mut out).map_err(io)?;
    out.flush().map_err(io)?;
    for e in &errors {
        eprintln!("error: {e}");
    }
    // Groups, settings and ion integrals go next to a CSV file; standard output has no room for more tables.
    if let Some(p) = &a.output {
        // Every file written, with the number of leading name columns `--excel-de` leaves as they are.
        let mut written = vec![(p.clone(), 3)];
        written.extend(table.write_groups_file(p).map_err(io)?.map(|f| (f, 2)));
        let stats = table.write_group_stats_file(p, Some("BLANK")).map_err(io)?;
        written.extend(stats.map(|f| (f, 3)));
        let create = |p: &Path| std::fs::File::create(p).map(std::io::BufWriter::new);
        let f = create(&settings_path(p)).map_err(io)?;
        ts.write_csv(&how, f).map_err(io)?;
        written.push((settings_path(p), 2));
        if !a.ion_integrals {
            // An older ions file would not match this results file.
            std::fs::remove_file(ions_path(p)).ok();
        } else {
            let by_column: Vec<&RunResult> = table.sources.iter().map(|i| &runs[*i]).collect();
            let f = create(&ions_path(p)).map_err(io)?;
            table.write_ions_csv(&by_column, &ts, f).map_err(io)?;
            written.push((ions_path(p), 2));
        }
        if a.excel_de {
            // The files are written in the standard format first; the German variant is a rewrite.
            for (f, text_cols) in &written {
                md_core::csv::excel_de_file(f, *text_cols).map_err(io)?;
            }
        }
    }
    eprintln!(
        "{} runs scanned in {:.2?} ({} failed): {} substances",
        runs.len(),
        started.elapsed(),
        errors.len(),
        table.rows.len()
    );
    Ok(errors.is_empty())
}

fn main() -> ExitCode {
    match parse(std::env::args_os().skip(1)) {
        Ok(None) => println!("{USAGE}"),
        Ok(Some(a)) => match run(a) {
            Ok(true) => {}
            Ok(false) => return ExitCode::from(2),
            Err(e) => {
                eprintln!("error: {e}");
                return ExitCode::FAILURE;
            }
        },
        Err(e) => {
            eprintln!("error: {e}\n\n{USAGE}");
            return ExitCode::FAILURE;
        }
    }
    ExitCode::SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(s: &str) -> Result<Option<Args>, String> {
        parse(s.split_whitespace().map(OsString::from))
    }

    #[test]
    fn parses_options_and_runs() {
        let a = args("-l lib.msl -j 2 --ions 5 a.cdf b.cdf")
            .unwrap()
            .unwrap();
        assert_eq!((a.parallel, a.ions, a.runs.len()), (2, Some(5), 2));
        assert!(args("--help").unwrap().is_none());
        // After `--` a name starting with a dash is a run.
        let a = args("-l lib.msl -- -x.cdf").unwrap().unwrap();
        assert_eq!(a.runs, [PathBuf::from("-x.cdf")]);
    }

    /// Two substances of 12 ions each, apart in time (as in the `md-core` batch test), written to
    /// `dir` as a run, a library, a settings file and a calibration; returns the argument line.
    fn study(dir: &std::path::Path) -> String {
        let (scans, ions) = (200, 24);
        let c = md_core::Chromatogram {
            scan_times_s: (0..scans).map(|s| s as f64 * 0.5).collect(),
            offsets: (0..=scans).map(|s| s * ions).collect(),
            mz: (0..scans * ions)
                .map(|k| 60.0 + (k % ions) as f32)
                .collect(),
            intensity: (0..scans * ions)
                .map(|k| {
                    let (s, i) = ((k / ions) as f32, k % ions);
                    let apex = if i < 12 { 60.0 } else { 140.0 };
                    let noise = ((k as u32).wrapping_mul(2_654_435_761) >> 16) % 41;
                    let peak =
                        5000.0 * (i % 12 + 1) as f32 * (-0.5 * ((s - apex) / 2.5).powi(2)).exp();
                    980.0 + noise as f32 + peak
                })
                .collect(),
            exact_masses: false,
            baseline: Default::default(),
        };
        std::fs::create_dir_all(dir).unwrap();
        md_core::andi::write(dir.join("ctrl_1.cdf"), &c).unwrap();
        let spectrum = |first: usize| -> String {
            (0..12)
                .map(|j| format!("({} {}) ", 60 + first + j, j + 1))
                .collect()
        };
        let msl = format!(
            "NAME: late, \"x\"\nNUM PEAKS: 12\n{}\n\nNAME: early\nNUM PEAKS: 12\n{}\n",
            spectrum(12),
            spectrum(0)
        );
        std::fs::write(dir.join("lib.msl"), msl).unwrap();
        std::fs::write(dir.join("a.settings"), GcmsSettings::default().to_bytes()).unwrap();
        std::fs::write(dir.join("a.calib"), "0 0\n100000 1000\n").unwrap();
        let p = |n: &str| dir.join(n).display().to_string();
        format!(
            "-l {} -s {} -o {} --replicates ^(.+)_1 --ions 2 --min-dist 1 --quality-cutoff 5 --min-ions 1 --min-reproducibility 0.5 {}",
            p("lib.msl"),
            p("a.settings"),
            p("out.csv"),
            p("ctrl_1.cdf")
        )
    }

    #[test]
    fn a_saved_cmp_is_a_run_too() {
        let dir = std::env::temp_dir().join(format!("md-cli-cmp-{}", std::process::id()));
        let line = study(&dir);
        let save = format!("{line} --save-cmp {}", dir.display());
        assert!(run(args(&save).unwrap().unwrap()).unwrap());
        let from_run = std::fs::read_to_string(dir.join("out.csv")).unwrap();

        let line = line.replace("ctrl_1.cdf", "ctrl_1.cmp");
        assert!(run(args(&format!("{line} --rescan")).unwrap().unwrap()).unwrap());
        let from_cmp = std::fs::read_to_string(dir.join("out.csv")).unwrap();
        assert_eq!(from_cmp, from_run.replace("ctrl_1.cdf", "ctrl_1.cmp"));
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn settings_and_ion_integrals_go_next_to_the_csv() {
        let dir = std::env::temp_dir().join(format!("md-cli-files-{}", std::process::id()));
        let line = study(&dir);
        assert!(run(args(&line).unwrap().unwrap()).unwrap());
        let settings = std::fs::read_to_string(dir.join("out.settings.csv")).unwrap();
        assert!(settings.contains("\nQUANT_ION_NUMBER,2\n"), "{settings}");
        assert!(settings.contains("\nUNTARGETED,false\n"), "{settings}");
        assert!(settings.contains("\nCALIBRATION,none\n"), "{settings}");
        assert!(settings.contains("\nLIBRARY,lib.msl\n"), "{settings}");
        assert!(
            settings.contains("\nIONS_ADDED_BY_RESCAN,0\n"),
            "{settings}"
        );
        assert!(!dir.join("out.ions.csv").exists());

        // The ion integrals are written on request: a line per substance and quant ion.
        assert!(run(args(&format!("--ion-integrals {line}")).unwrap().unwrap()).unwrap());
        let ions = std::fs::read_to_string(dir.join("out.ions.csv")).unwrap();
        assert!(
            ions.starts_with("\u{feff}Substance,Quant_ion,ctrl_1.cdf\nearly,"),
            "{ions}"
        );
        assert_eq!(ions.lines().count(), 1 + 2 * 2, "{ions}");
        assert!(args("-l lib.msl --ion-integrals a.cdf").is_err());
        // A later export without the flag removes the ions file, which no longer matches.
        assert!(run(args(&line).unwrap().unwrap()).unwrap());
        assert!(!dir.join("out.ions.csv").exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn excel_de_rewrites_the_whole_export_and_the_default_stays() {
        let dir = std::env::temp_dir().join(format!("md-cli-de-{}", std::process::id()));
        let line = study(&dir);
        assert!(run(args(&format!("--excel-de {line}")).unwrap().unwrap()).unwrap());
        let csv = std::fs::read_to_string(dir.join("out.csv")).unwrap();
        assert!(
            csv.starts_with("\u{feff}Substance;RT_s;Quant_ions;ctrl_1.cdf\n"),
            "{csv}"
        );
        let groups = std::fs::read_to_string(dir.join("out.groups.csv")).unwrap();
        assert_eq!(groups, "\u{feff}Run;Group\nctrl_1.cdf;ctrl\n");
        assert!(args("-l lib.msl --excel-de a.cdf").is_err());
        assert!(run(args(&line).unwrap().unwrap()).unwrap());
        let csv = std::fs::read_to_string(dir.join("out.csv")).unwrap();
        assert!(csv.starts_with("\u{feff}Substance,RT_s,"), "{csv}");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn writes_the_table_of_a_study() {
        let dir = std::env::temp_dir().join(format!("md-cli-study-{}", std::process::id()));
        let line = study(&dir);
        assert!(run(args(&line).unwrap().unwrap()).unwrap());
        let csv = std::fs::read_to_string(dir.join("out.csv")).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], "\u{feff}Substance,RT_s,Quant_ions,ctrl_1.cdf");
        let groups = std::fs::read_to_string(dir.join("out.groups.csv")).unwrap();
        assert_eq!(groups, "\u{feff}Run,Group\nctrl_1.cdf,ctrl\n");
        let stats = std::fs::read_to_string(dir.join("out.groupstats.csv")).unwrap();
        assert!(
            stats.starts_with(
                "\u{feff}Substance,RT_s,Quant_ions,Mean_ctrl,Rel_Std_Err_ctrl,Hits_ctrl,P_ANOVA\nearly,"
            ),
            "{stats}"
        );
        assert!(lines[1].starts_with("early,"), "{csv}");
        assert!(lines[2].starts_with("\"late, \"\"x\"\"\","), "{csv}");
        assert_eq!(lines.len(), 3);

        // With a calibration the library's RI (none here) decides; the run still succeeds.
        let calib = format!("{line} -c {}", dir.join("a.calib").display());
        assert!(run(args(&calib).unwrap().unwrap()).unwrap());

        // An unreadable run is reported (exit code 2), the others are kept.
        let bad = format!("{line} {}", dir.join("missing.cdf").display());
        assert!(!run(args(&bad).unwrap().unwrap()).unwrap());
        assert!(
            run(args(&line.replace("^(.+)_1", "(")).unwrap().unwrap())
                .unwrap_err()
                .contains("--replicates")
        );
        let nowhere = line.replace("out.csv", "no/such/dir/out.csv");
        assert!(
            run(args(&nowhere).unwrap().unwrap())
                .unwrap_err()
                .contains("no/such")
        );
        // Settings, calibration and library are checked.
        let fails = |line: &str| run(args(line).unwrap().unwrap()).is_err();
        std::fs::write(dir.join("a.settings"), b"too short").unwrap();
        assert!(fails(&line));
        std::fs::write(dir.join("a.settings"), GcmsSettings::default().to_bytes()).unwrap();
        std::fs::write(dir.join("a.calib"), "x").unwrap();
        assert!(fails(&calib));
        std::fs::remove_file(dir.join("lib.msl")).unwrap();
        assert!(fails(&line));
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn a_sample_table_sets_the_groups() {
        let dir = std::env::temp_dir().join(format!("md-cli-samples-{}", std::process::id()));
        let line = study(&dir);
        let table = dir.join("samples.txt");
        std::fs::write(&table, "NOTE\tCDF_FILE\tTIME_POINT\nx\tctrl_1.CDF\t5\n").unwrap();
        let with = |opts: &str| format!("{line} --sample-table {} {opts}", table.display());
        // The table wins over --replicates for the runs it names; the columns are chosen by name.
        let cols = "--sample-file-column cdf_file --sample-group-column TIME_POINT";
        assert!(run(args(&with(cols)).unwrap().unwrap()).unwrap());
        let groups = std::fs::read_to_string(dir.join("out.groups.csv")).unwrap();
        assert_eq!(groups, "\u{feff}Run,Group\nctrl_1.cdf,5\n");
        // Without names CDF_FILE and TIME_POINT are used.
        assert!(run(args(&with("")).unwrap().unwrap()).unwrap());
        let groups = std::fs::read_to_string(dir.join("out.groups.csv")).unwrap();
        assert_eq!(groups, "\u{feff}Run,Group\nctrl_1.cdf,5\n");
        let err = run(args(&with("--sample-group-column nope")).unwrap().unwrap()).unwrap_err();
        assert!(err.contains("no column nope"), "{err}");
        assert!(args(&format!("{line} --sample-group-column x")).is_err());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn untargeted_adds_the_unknown_substance_and_options_are_checked() {
        let dir = std::env::temp_dir().join(format!("md-cli-unknown-{}", std::process::id()));
        let line = study(&dir);
        // Only "early" is in the library, so "late" can only come from the pooled unknowns.
        let early: String = (0..12)
            .map(|j| format!("({} {}) ", 60 + j, j + 1))
            .collect();
        std::fs::write(
            dir.join("lib.msl"),
            format!(
                "NAME: early
NUM PEAKS: 12
{early}
"
            ),
        )
        .unwrap();
        let table = |extra: &str| {
            let l = format!("{line} {extra}");
            assert!(run(args(&l).unwrap().unwrap()).unwrap());
            std::fs::read_to_string(dir.join("out.csv")).unwrap()
        };
        assert_eq!(table("").lines().count(), 2, "header and early");
        let csv = table(
            "--untargeted --rescan --baseline-correction --exclude-ions 73,147-149 --min-sn 1 --max-discrepancy 1000",
        );
        let names: Vec<&str> = csv
            .lines()
            .skip(1)
            .map(|l| l.split(',').next().unwrap())
            .collect();
        assert_eq!(names, ["early", "Unknown 1"], "{csv}");
        let settings = std::fs::read_to_string(dir.join("out.settings.csv")).unwrap();
        assert!(settings.contains("\nUNTARGETED,true\n"), "{settings}");
        assert!(
            settings.contains("\nBASELINE_CORRECTION,true\n"),
            "{settings}"
        );
        assert!(settings.contains("\nUNKNOWN_SUBSTANCES,1\n"), "{settings}");
        assert!(settings.contains("\nRUNS_SCANNED,1\n"), "{settings}");

        let shift = |extra: &str| run(args(&format!("{line} {extra}")).unwrap().unwrap());
        assert!(
            shift("--shift-reference early")
                .unwrap_err()
                .contains("calibration")
        );
        let calib = format!("-c {}", dir.join("a.calib").display());
        assert!(
            shift(&format!("{calib} --shift-reference nobody"))
                .unwrap_err()
                .contains("nobody")
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn parses_mass_ranges() {
        assert_eq!(
            ranges("73, 147-149").unwrap(),
            [(73.0, 73.0), (147.0, 149.0)]
        );
        assert!(ranges("").unwrap().is_empty());
        assert!(ranges("x").is_err() && ranges("5-").is_err() && ranges("-1").is_err());
    }

    #[test]
    fn rejects_bad_input() {
        assert!(args("a.cdf").is_err(), "library missing");
        assert!(args("-l lib.msl").is_err(), "no runs");
        assert!(args("-l lib.msl --bogus a.cdf").is_err());
        assert!(args("-l lib.msl -j x a.cdf").is_err());
        assert!(args("-l lib.msl a.cdf --ions").is_err());
        assert!(args("-l lib.msl --ions 0 a.cdf").is_err());
        assert!(args("-l lib.msl --min-reproducibility 5 a.cdf").is_err());
        assert!(args("-l lib.msl -j 1e12 a.cdf").is_err());
    }
}
