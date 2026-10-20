//! Headless batch quantification for big studies: finds the library substances in many runs and writes
//! the substance × run table of the app's batch window as CSV. Depends on `md-core` only.

use std::ffi::OsString;
use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicBool, AtomicUsize};
use std::time::Instant;

use md_core::batch::{TableSettings, scan_runs, table};
use md_core::library;
use md_core::ri::RiCalibration;
use md_core::search::PreparedLibrary;
use md_core::settings::GcmsSettings;

const USAGE: &str = "\
Usage: md-cli --library FILE [options] RUN...

Finds the substances of a library (.mdlib, .lbr, .msl, NIST .msp or .sdf) in many runs (.cdf, .nc or
old .bin) and writes the substance x run table as CSV. Exit code 2: the table is written, but some runs could not be read.

Options:
  -l, --library FILE          library (required)
  -s, --settings FILE         old MetaboliteDetector settings file (default: the old defaults)
  -c, --calibration FILE      RI calibration (.calib); without it the search uses the spectrum only
  -o, --output FILE           CSV file (default: standard output)
  -j, --parallel N            runs in memory at once (default 4)
      --replicates REGEX      replicate group = first capture group of the file name
      --min-reproducibility F keep a substance only if a group quantifies it in this share of its runs
      --min-ions N            quant ions that must be found for a value
      --ions N                quant ions per substance (default: QUANT_ION_NUMBER of the settings)
      --min-dist F            min. distance of quant ions (default: QUANT_MIN_DIST_IONS)
      --quality-cutoff F      quant ion quality cut-off (default: QUANT_QUALITY_CUTOFF)
  -h, --help                  this text";

struct Args {
    library: PathBuf,
    settings: Option<PathBuf>,
    calibration: Option<PathBuf>,
    output: Option<PathBuf>,
    parallel: usize,
    replicates: String,
    min_reproducibility: Option<f64>,
    min_ions: Option<usize>,
    ions: Option<usize>,
    min_dist: Option<f32>,
    quality_cutoff: Option<f32>,
    runs: Vec<PathBuf>,
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
        min_reproducibility: None,
        min_ions: None,
        ions: None,
        min_dist: None,
        quality_cutoff: None,
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
            "--min-reproducibility" => a.min_reproducibility = Some(num(0.0, 1.0)?),
            "--min-ions" => a.min_ions = Some(num(1.0, 20.0)? as usize),
            "--ions" => a.ions = Some(num(1.0, 20.0)? as usize),
            "--min-dist" => a.min_dist = Some(num(0.0, 100.0)? as f32),
            "--quality-cutoff" => a.quality_cutoff = Some(num(0.0, 100.0)? as f32),
            _ => return Err(format!("unknown option {name}")),
        }
    }
    a.library = library.ok_or("--library is required")?;
    if a.runs.is_empty() {
        return Err("no runs given".into());
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
    ts.quant.ion_number = a.ions.unwrap_or(ts.quant.ion_number);
    ts.quant.min_dist = a.min_dist.unwrap_or(ts.quant.min_dist);
    ts.quant.quality_cutoff = a.quality_cutoff.unwrap_or(ts.quant.quality_cutoff);
    ts.min_reproducibility = a.min_reproducibility.unwrap_or(ts.min_reproducibility);
    ts.min_ions = a.min_ions.unwrap_or(ts.min_ions);

    // Open the output first: a wrong path must not cost a whole scan.
    let mut out: Box<dyn Write> = match &a.output {
        Some(p) => Box::new(std::io::BufWriter::new(
            std::fs::File::create(p).map_err(|e| format!("{}: {e}", p.display()))?,
        )),
        None => Box::new(std::io::stdout().lock()),
    };
    let started = Instant::now();
    let (done, cancel) = (AtomicUsize::new(0), AtomicBool::new(false));
    let (runs, errors) = scan_runs(
        &a.runs,
        &prepared,
        &settings.scan(calibration),
        a.parallel,
        &done,
        &cancel,
    )
    .ok_or("cancelled")?;
    let table = table(&runs, &compounds, &ts).map_err(|e| format!("--replicates: {e}"))?;

    let io = |e: std::io::Error| format!("output: {e}");
    table.write_csv(&mut out).map_err(io)?;
    out.flush().map_err(io)?;
    for e in &errors {
        eprintln!("error: {e}");
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
    fn writes_the_table_of_a_study() {
        let dir = std::env::temp_dir().join(format!("md-cli-study-{}", std::process::id()));
        let line = study(&dir);
        assert!(run(args(&line).unwrap().unwrap()).unwrap());
        let csv = std::fs::read_to_string(dir.join("out.csv")).unwrap();
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines[0], "Substance,RT_s,Quant_ions,ctrl_1.cdf");
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
