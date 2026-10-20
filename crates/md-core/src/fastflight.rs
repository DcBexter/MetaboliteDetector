//! FastFlight TOF raw data (`MsData.FFT` + `MsData.FFC`), read like the old `gcms::FastFlightImporter`.
//! Layout: docs/reverse-engineering/file-formats.md §4.3. The m/z of a point comes from the TOF
//! calibration polynomial (`.mzcalib`, see [`crate::tof`]).

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use std::path::{Path, PathBuf};

use crate::Chromatogram;
use crate::andi;
use crate::tof::{TofCalibration, TofError};

#[derive(Debug, thiserror::Error)]
pub enum FastFlightError {
    #[error("cannot read {0}: {1}")]
    Io(PathBuf, std::io::Error),
    #[error("corrupt FastFlight data: {0}")]
    Corrupt(String),
    #[error("{0}")]
    Calibration(#[from] TofError),
    #[error("no .mzcalib TOF calibration file next to {0}")]
    NoCalibration(PathBuf),
}

const FFT: &str = "MsData.FFT";
const FFC: &str = "MsData.FFC";
/// `FastFlightImporter::initialization`: an FFC time of 0 means this factor.
const DEFAULT_TIME_FACTOR: f64 = 2.5;
const FRAME_MARKER: u32 = 0xFF_FFFF;
const BASE_MARKER: u32 = 0xFF_8000;
const POSITION_MASK: u32 = 0xFF_0000;
/// FFC bytes per scan: 4 skipped, `n` i32, 4 skipped, TIC u64, FFT offset u64, time factor f64.
const FFC_ENTRY: usize = 36;

/// `dir/name` in the letter case the directory has, else as given.
fn named(dir: &Path, name: &str) -> PathBuf {
    std::fs::read_dir(dir)
        .into_iter()
        .flatten()
        .flatten()
        .map(|e| e.path())
        .find(|p| {
            p.file_name()
                .and_then(std::ffi::OsStr::to_str)
                .is_some_and(|n| n.eq_ignore_ascii_case(name))
        })
        .unwrap_or_else(|| dir.join(name))
}

/// True for the two files of a FastFlight run (any letter case).
pub fn is_fastflight(path: &Path) -> bool {
    path.file_name()
        .and_then(std::ffi::OsStr::to_str)
        .is_some_and(|n| n.eq_ignore_ascii_case(FFT) || n.eq_ignore_ascii_case(FFC))
}

/// The `.mzcalib` file lying next to the run (the first by name); the old import dialog asked for it.
pub fn find_calibration(path: &Path) -> Result<PathBuf, FastFlightError> {
    let dir = path.parent().unwrap_or(Path::new("."));
    let mut found: Vec<PathBuf> = std::fs::read_dir(dir)
        .map_err(|e| FastFlightError::Io(dir.to_path_buf(), e))?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            p.extension()
                .is_some_and(|e| e.eq_ignore_ascii_case("mzcalib"))
        })
        .collect();
    found.sort();
    found
        .into_iter()
        .next()
        .ok_or_else(|| FastFlightError::NoCalibration(path.to_path_buf()))
}

/// Reads a run; `path` names `MsData.FFT`, `MsData.FFC` or their directory. Scan times become
/// seconds, m/z comes from `calibration`, and intensities are what `FastFlightImporter::importData`
/// stores: the smallest intensity of a scan subtracted, the first point dropped.
pub fn read(path: &Path, calibration: &TofCalibration) -> Result<Chromatogram, FastFlightError> {
    let (mut times, mut starts, mut counts) = (Vec::new(), Vec::new(), Vec::new());
    let (mut mz, mut intensity) = (Vec::new(), Vec::new());
    for_each_scan(path, |time_ms, bins, values| {
        times.push(time_ms / 1000.0);
        starts.push(mz.len() as f64);
        counts.push(bins.len() as f64);
        mz.extend(bins.iter().map(|b| calibration.mz(*b)));
        intensity.extend(values.iter().map(|v| *v as f32));
    })?;
    let mut run = andi::build(times, &starts, &counts, mz, intensity)
        .map_err(|e| FastFlightError::Corrupt(e.to_string()))?;
    run.exact_masses = true;
    Ok(run)
}

/// The run summed over all scans, per TOF bin: (bin, total intensity) in bin order. This is what
/// the calibration window shows, since the bins are what the calibration maps to m/z.
pub fn bin_spectrum(path: &Path) -> Result<Vec<(f32, f64)>, FastFlightError> {
    let mut sums = std::collections::BTreeMap::<u32, f64>::new();
    for_each_scan(path, |_, bins, values| {
        for (b, v) in bins.iter().zip(&values) {
            *sums.entry(*b as u32).or_default() += f64::from(*v);
        }
    })?;
    Ok(sums.into_iter().map(|(b, v)| (b as f32, v)).collect())
}

/// Walks the scans of a run in time order: `visit(time in ms, TOF bins, intensities)` for each scan
/// with signal, with the lowest intensity subtracted and the first point dropped.
fn for_each_scan(
    path: &Path,
    mut visit: impl FnMut(f64, Vec<f32>, Vec<i32>),
) -> Result<(), FastFlightError> {
    let dir = if path.is_dir() {
        path
    } else {
        path.parent().unwrap_or(Path::new("."))
    };
    let open = |name: &str| {
        let p = named(dir, name);
        File::open(&p).map_err(|e| FastFlightError::Io(p, e))
    };
    let index = read_index(&mut BufReader::new(open(FFC)?))?;
    let fft_file = open(FFT)?;
    // Every point takes a 3-byte word, so more points than that mean index entries share data.
    let max_points = fft_file.metadata().map_or(0, |m| m.len() as usize / 3);
    let mut fft = BufReader::new(fft_file);
    let mut points = 0;
    for scan in index {
        fft.seek(SeekFrom::Start(scan.offset))
            .map_err(|e| FastFlightError::Io(dir.join(FFT), e))?;
        let (mut bins, mut values, tic) = read_spectrum(&mut fft).map_err(|e| match e {
            Fail::Io(e) => FastFlightError::Io(dir.join(FFT), e),
            Fail::Corrupt(s) => {
                FastFlightError::Corrupt(format!("scan at {} ms: {s}", scan.time_ms))
            }
        })?;
        // The old importer compared the trailer with the FFC value and refused the run.
        if f64::from(tic) != scan.tic as f64 {
            return Err(FastFlightError::Corrupt(
                "TIC mismatch! TICs in the file .FFT and .FFC do not match.".into(),
            ));
        }
        // Empty scans carry no signal; the old importer left them out.
        if tic == 0 {
            continue;
        }
        if !bins.is_empty() {
            bins.remove(0);
            values.remove(0);
        }
        points += bins.len();
        if points > max_points {
            return Err(FastFlightError::Corrupt(
                "the index points several scans at the same data".into(),
            ));
        }
        visit(scan.time_ms, bins, values);
    }
    Ok(())
}

struct Scan {
    time_ms: f64,
    offset: u64,
    tic: u64,
}

/// `FastFlightImporter::initialization`, ordered by time. A time seen twice keeps its first scan.
fn read_index(ffc: &mut impl Read) -> Result<Vec<Scan>, FastFlightError> {
    let short = |e: std::io::Error| FastFlightError::Corrupt(format!("{FFC}: {e}"));
    let mut count = [0u8; 4];
    ffc.read_exact(&mut count).map_err(short)?;
    let count = i32::from_le_bytes(count);
    let mut scans = Vec::new();
    let mut seen = BTreeSet::new();
    let mut entry = [0u8; FFC_ENTRY];
    for _ in 0..count.max(0) {
        ffc.read_exact(&mut entry).map_err(short)?;
        let n = i32::from_le_bytes(entry[4..8].try_into().unwrap());
        let tic = u64::from_le_bytes(entry[12..20].try_into().unwrap());
        let offset = u64::from_le_bytes(entry[20..28].try_into().unwrap());
        let factor = f64::from_le_bytes(entry[28..36].try_into().unwrap());
        let factor = if factor == 0.0 {
            DEFAULT_TIME_FACTOR
        } else {
            factor
        };
        let time_ms = factor * f64::from(n);
        if n < 0 || !time_ms.is_finite() {
            return Err(FastFlightError::Corrupt(format!(
                "{FFC}: invalid time factor"
            )));
        }
        if seen.insert(time_ms.to_bits()) {
            scans.push(Scan {
                time_ms,
                offset,
                tic,
            });
        }
    }
    scans.sort_by(|a, b| a.time_ms.total_cmp(&b.time_ms));
    Ok(scans)
}

enum Fail {
    Io(std::io::Error),
    Corrupt(String),
}

fn word(r: &mut impl Read) -> Result<u32, Fail> {
    let mut b = [0u8; 3];
    r.read_exact(&mut b).map_err(|e| match e.kind() {
        std::io::ErrorKind::UnexpectedEof => Fail::Corrupt("end of file inside a record".into()),
        _ => Fail::Io(e),
    })?;
    Ok(u32::from_le_bytes([b[0], b[1], b[2], 0]))
}

/// `FastFlightImporter::readSpectrum<float>` with the importer's settings: (TOF bin, intensity)
/// points of one spectrum and its TIC. A spectrum is one or more frames, the last one flagged final.
/// Words are 24 bit, little endian.
fn read_spectrum(r: &mut impl Read) -> Result<(Vec<f32>, Vec<i32>, i32), Fail> {
    let (mut bins, mut values) = (Vec::new(), Vec::new());
    let (mut base, mut bin) = (0u32, 0u32);
    let mut lowest = i32::MAX;
    loop {
        if word(r)? != FRAME_MARKER {
            return Err(Fail::Corrupt("expected new frame header".into()));
        }
        let header = word(r)?;
        let time_block = header & 0x4000 != 0;
        // The old code tests the final flag only on frames without a time word.
        let last = !time_block && header & 0x8000 != 0;
        let mut skip = 2;
        if time_block {
            // Start of a spectrum: its acquisition time follows (unused here).
            bins.clear();
            values.clear();
            word(r)?;
            skip = 3;
        } else if last {
            skip = 4;
        }
        let words = (header & 0x3FFF)
            .checked_sub(skip)
            .ok_or_else(|| Fail::Corrupt("frame shorter than its header".into()))?;
        let mut i = 0;
        while i < words {
            let w = word(r)?;
            if w == BASE_MARKER {
                base = word(r)?;
                i += 1;
            } else if w & POSITION_MASK == POSITION_MASK {
                bin = w & !POSITION_MASK;
            } else {
                bins.push((base + bin) as f32);
                values.push(w as i32);
                lowest = lowest.min(w as i32);
                bin += 1;
            }
            i += 1;
        }
        if last {
            break;
        }
    }
    let tic = (word(r)? as i32).wrapping_add((word(r)? << 24) as i32);
    values.iter_mut().for_each(|v| *v -= lowest);
    Ok((bins, values, tic))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn words(ws: &[u32]) -> Vec<u8> {
        ws.iter()
            .flat_map(|w| w.to_le_bytes()[..3].to_vec())
            .collect()
    }

    /// A spectrum in two frames: base offset 100, positions 3 and 20, intensities 5 7 | 9 6.
    fn two_frames() -> Vec<u8> {
        words(&[
            // Time frame: header, time word, then 5 data words (marker pair, position, 2 values).
            FRAME_MARKER,
            0x4000 | 8,
            1234,
            BASE_MARKER,
            100,
            POSITION_MASK | 3,
            5,
            7,
            // Final frame: header, 3 data words, then the TIC trailer.
            FRAME_MARKER,
            0x8000 | 7,
            POSITION_MASK | 20,
            9,
            6,
            27,
            0,
        ])
    }

    #[test]
    fn reads_frames_of_one_spectrum() {
        let bytes = two_frames();
        let (bins, values, tic) = read_spectrum(&mut &bytes[..]).ok().unwrap();
        assert_eq!(bins, [103.0, 104.0, 120.0, 121.0]);
        assert_eq!(values, [0, 2, 4, 1]); // lowest intensity 5 subtracted
        assert_eq!(tic, 27);
    }

    #[test]
    fn truncated_data_is_an_error() {
        let bytes = two_frames();
        for len in 0..bytes.len() {
            assert!(read_spectrum(&mut &bytes[..len]).is_err(), "cut at {len}");
        }
        let mut wrong = bytes.clone();
        wrong[0] = 0;
        assert!(read_spectrum(&mut &wrong[..]).is_err());
    }
}
