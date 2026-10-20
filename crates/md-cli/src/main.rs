//! Headless batch quantification for big studies: finds the library substances in many runs and writes
//! the substance × run table of the app's batch window as CSV. Depends on `md-core` only.

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

Finds the substances of a library (.lbr or .msl) in many runs (.cdf, .nc or old .bin) and writes the
substance x run table as CSV. Exit code 2: the table is written, but some runs could not be read.

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
fn parse(args: impl IntoIterator<Item = String>) -> Result<Option<Args>, String> {
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
    while let Some(arg) = args.next() {
        let mut value = || args.next().ok_or_else(|| format!("{arg}: value missing"));
        let mut num = || {
            let v = value()?;
            v.parse::<f64>()
                .ok()
                .filter(|n| n.is_finite() && *n >= 0.0)
                .ok_or_else(|| format!("{arg}: not a non-negative number: {v}"))
        };
        match arg.as_str() {
            "-h" | "--help" => return Ok(None),
            "-l" | "--library" => library = Some(PathBuf::from(value()?)),
            "-s" | "--settings" => a.settings = Some(value()?.into()),
            "-c" | "--calibration" => a.calibration = Some(value()?.into()),
            "-o" | "--output" => a.output = Some(value()?.into()),
            "-j" | "--parallel" => a.parallel = num()? as usize,
            "--replicates" => a.replicates = value()?,
            "--min-reproducibility" => a.min_reproducibility = Some(num()?),
            "--min-ions" => a.min_ions = Some(num()? as usize),
            "--ions" => a.ions = Some(num()? as usize),
            "--min-dist" => a.min_dist = Some(num()? as f32),
            "--quality-cutoff" => a.quality_cutoff = Some(num()? as f32),
            s if s.starts_with('-') && s.len() > 1 => return Err(format!("unknown option {s}")),
            _ => a.runs.push(arg.into()),
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
    let prepared = PreparedLibrary::new(&compounds, &settings.search());

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
    match &a.output {
        Some(p) => {
            let mut w = std::io::BufWriter::new(
                std::fs::File::create(p).map_err(|e| format!("{}: {e}", p.display()))?,
            );
            table.write_csv(&mut w).map_err(io)?;
            w.flush().map_err(io)?;
        }
        None => table.write_csv(&mut std::io::stdout().lock()).map_err(io)?,
    }
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
    match parse(std::env::args().skip(1)) {
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
        parse(s.split_whitespace().map(String::from))
    }

    #[test]
    fn parses_options_and_runs() {
        let a = args("-l lib.msl -j 2 --ions 5 a.cdf b.cdf")
            .unwrap()
            .unwrap();
        assert_eq!((a.parallel, a.ions, a.runs.len()), (2, Some(5), 2));
        assert!(args("--help").unwrap().is_none());
    }

    #[test]
    fn rejects_bad_input() {
        assert!(args("a.cdf").is_err(), "library missing");
        assert!(args("-l lib.msl").is_err(), "no runs");
        assert!(args("-l lib.msl --bogus a.cdf").is_err());
        assert!(args("-l lib.msl -j x a.cdf").is_err());
        assert!(args("-l lib.msl a.cdf --ions").is_err());
    }
}
