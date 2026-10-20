//! ANDI-MS (ASTM E2078 / AIA) NetCDF-3 import: the vendor-neutral GC/MS export format
//! that MetaboliteDetector 3.x imported via its "NetCDF import" dialog.

use std::cell::RefCell;
use std::io::{Cursor, Seek, SeekFrom, Write};
use std::path::Path;
use std::rc::Rc;

use netcdf3::{DataSet, DataVector, FileReader, FileWriter, Version};

use crate::Chromatogram;

#[derive(Debug, thiserror::Error)]
pub enum AndiError {
    #[error("cannot read NetCDF file: {0:?}")]
    Read(netcdf3::ReadError),
    #[error("cannot write NetCDF file: {0:?}")]
    Write(netcdf3::WriteError),
    #[error("missing ANDI-MS variable `{0}`")]
    MissingVariable(&'static str),
    #[error("inconsistent ANDI-MS data: {0}")]
    Inconsistent(String),
}

/// Reads an ANDI-MS `.cdf` file into memory.
// shortcut: whole file in RAM (f32, ~8 bytes per data point); add streaming/memory-mapped
// access like the old GCMSDiskScan if real runs exceed workstation memory.
pub fn read(path: impl AsRef<Path>) -> Result<Chromatogram, AndiError> {
    let path = path.as_ref();
    // netcdf3 issues a seek + read per record; on a plain File that is millions of syscalls.
    let bytes = std::fs::read(path).map_err(|e| AndiError::Read(e.into()))?;
    let mut reader = FileReader::open_seek_read(
        &path.to_string_lossy(),
        Box::new(std::io::Cursor::new(bytes)),
    )
    .map_err(AndiError::Read)?;
    let mut var = |name: &'static str| -> Result<(DataVector, f64), AndiError> {
        let scale = reader
            .data_set()
            .get_var_attr(name, "scale_factor")
            .and_then(first_number)
            .unwrap_or(1.0);
        let data = reader
            .read_var(name)
            .map_err(|_| AndiError::MissingVariable(name))?;
        Ok((data, scale))
    };

    let (times, scale) = var("scan_acquisition_time")?;
    let scan_times_s = to_f64(&times, scale);
    let (index, _) = var("scan_index")?;
    let (count, _) = var("point_count")?;
    let (mz, mz_scale) = var("mass_values")?;
    let (intensity, intensity_scale) = var("intensity_values")?;
    build(
        scan_times_s,
        &to_f64(&index, 1.0),
        &to_f64(&count, 1.0),
        to_f32(mz, mz_scale),
        to_f32(intensity, intensity_scale),
    )
}

/// Writes the variables `read` needs, which is the core of ANDI-MS (raw data group).
pub fn write(path: impl AsRef<Path>, c: &Chromatogram) -> Result<(), AndiError> {
    write_with_tic(path, c, &c.tic())
}

/// Like [`write`], with the total intensity of each scan given: a centroided export keeps the
/// TIC of the profile scans it came from.
pub fn write_with_tic(
    path: impl AsRef<Path>,
    c: &Chromatogram,
    tic: &[f64],
) -> Result<(), AndiError> {
    let invalid = |e| AndiError::Inconsistent(format!("{e:?}"));
    let i32s = |v: Vec<usize>| -> Result<Vec<i32>, AndiError> {
        v.into_iter()
            .map(|x| {
                i32::try_from(x)
                    .map_err(|_| AndiError::Inconsistent(format!("{x} points do not fit ANDI-MS")))
            })
            .collect()
    };
    let scans = c.scan_count();
    if tic.len() != scans {
        return Err(AndiError::Inconsistent("one TIC per scan needed".into()));
    }
    // NetCDF-3 forbids fixed dimensions of size 0.
    if scans == 0 || c.mz.is_empty() {
        return Err(AndiError::Inconsistent("the run has no data points".into()));
    }
    let index = i32s(c.offsets[..scans].to_vec())?;
    let count = i32s(c.offsets.windows(2).map(|w| w[1] - w[0]).collect())?;
    let mut ds = DataSet::new();
    ds.add_fixed_dim("scan_number", scans).map_err(invalid)?;
    ds.add_fixed_dim("point_number", c.mz.len())
        .map_err(invalid)?;
    for name in ["scan_acquisition_time", "total_intensity"] {
        ds.add_var_f64(name, &["scan_number"]).map_err(invalid)?;
    }
    for name in ["scan_index", "point_count"] {
        ds.add_var_i32(name, &["scan_number"]).map_err(invalid)?;
    }
    for name in ["mass_values", "intensity_values"] {
        ds.add_var_f32(name, &["point_number"]).map_err(invalid)?;
    }
    // netcdf3 writes value by value; buffer in memory (FileWriter::close would hide a
    // BufWriter's final flush error) and write the file once.
    let buf = SharedBuf::default();
    let mut w = FileWriter::open_seek_write("", Box::new(buf.clone())).map_err(AndiError::Write)?;
    // Classic NetCDF-3 caps file offsets at 2 GiB; the 64-bit-offset variant lifts that.
    w.set_def(&ds, Version::Offset64Bit, 0)
        .map_err(AndiError::Write)?;
    w.write_var_f64("scan_acquisition_time", &c.scan_times_s)
        .and_then(|_| w.write_var_f64("total_intensity", tic))
        .and_then(|_| w.write_var_i32("scan_index", &index))
        .and_then(|_| w.write_var_i32("point_count", &count))
        .and_then(|_| w.write_var_f32("mass_values", &c.mz))
        .and_then(|_| w.write_var_f32("intensity_values", &c.intensity))
        .and_then(|_| w.close())
        .map_err(AndiError::Write)?;
    // Temp file + rename: a failed write never leaves the target (maybe the source .cdf) truncated.
    let path = path.as_ref();
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".tmp");
    let result =
        std::fs::write(&tmp, buf.0.take().into_inner()).and_then(|_| std::fs::rename(&tmp, path));
    if result.is_err() {
        std::fs::remove_file(&tmp).ok();
    }
    result.map_err(|e| AndiError::Write(e.into()))
}

#[derive(Clone, Default)]
struct SharedBuf(Rc<RefCell<Cursor<Vec<u8>>>>);

impl Write for SharedBuf {
    fn write(&mut self, data: &[u8]) -> std::io::Result<usize> {
        self.0.borrow_mut().write(data)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl Seek for SharedBuf {
    fn seek(&mut self, pos: SeekFrom) -> std::io::Result<u64> {
        self.0.borrow_mut().seek(pos)
    }
}

fn first_number(attr: &netcdf3::Attribute) -> Option<f64> {
    attr.get_f64()
        .map(|v| v.first().copied())
        .or_else(|| attr.get_f32().map(|v| v.first().map(|x| *x as f64)))
        .or_else(|| attr.get_i32().map(|v| v.first().map(|x| *x as f64)))
        .or_else(|| attr.get_i16().map(|v| v.first().map(|x| *x as f64)))
        .flatten()
}

fn to_f64(data: &DataVector, scale: f64) -> Vec<f64> {
    let f = |x: f64| x * scale;
    match data {
        DataVector::I8(v) => v.iter().map(|x| f(*x as f64)).collect(),
        DataVector::U8(v) => v.iter().map(|x| f(*x as f64)).collect(),
        DataVector::I16(v) => v.iter().map(|x| f(*x as f64)).collect(),
        DataVector::I32(v) => v.iter().map(|x| f(*x as f64)).collect(),
        DataVector::F32(v) => v.iter().map(|x| f(*x as f64)).collect(),
        DataVector::F64(v) => v.iter().map(|x| f(*x)).collect(),
    }
}

/// Bulk arrays: reuses the file's buffer when it already is unscaled `f32`.
fn to_f32(data: DataVector, scale: f64) -> Vec<f32> {
    match data {
        DataVector::F32(v) if scale == 1.0 => v,
        other => to_f64(&other, scale)
            .into_iter()
            .map(|x| x as f32)
            .collect(),
    }
}

/// Validates the ANDI arrays (untrusted file input) and converts them to offsets.
pub(crate) fn build(
    scan_times_s: Vec<f64>,
    scan_index: &[f64],
    point_count: &[f64],
    mut mz: Vec<f32>,
    mut intensity: Vec<f32>,
) -> Result<Chromatogram, AndiError> {
    let n = scan_times_s.len();
    if scan_index.len() != n || point_count.len() != n {
        return Err(AndiError::Inconsistent(format!(
            "{n} scan times but {} scan indices and {} point counts",
            scan_index.len(),
            point_count.len()
        )));
    }
    if scan_times_s.windows(2).any(|w| w[1] < w[0]) {
        return Err(AndiError::Inconsistent(
            "scan times are not ascending".into(),
        ));
    }
    if mz.len() != intensity.len() {
        return Err(AndiError::Inconsistent(format!(
            "{} masses but {} intensities",
            mz.len(),
            intensity.len()
        )));
    }
    let mut offsets = Vec::with_capacity(n + 1);
    for (s, (start, count)) in scan_index.iter().zip(point_count).enumerate() {
        let (start, end) = (*start as usize, (*start + *count) as usize);
        if *count < 0.0 || offsets.last().is_some_and(|prev| *prev != start) || end > mz.len() {
            return Err(AndiError::Inconsistent(format!(
                "scan {s} points to values {start}..{end} of {}",
                mz.len()
            )));
        }
        if offsets.is_empty() {
            offsets.push(start);
        }
        offsets.push(end);
    }
    if offsets.is_empty() {
        offsets.push(0);
    }
    // Lookups binary-search m/z within a scan; the old app also sorted unsorted scans on load.
    for w in offsets.windows(2) {
        let (m, i) = (&mut mz[w[0]..w[1]], &mut intensity[w[0]..w[1]]);
        if !m.is_sorted() {
            let mut pairs: Vec<(f32, f32)> = m.iter().copied().zip(i.iter().copied()).collect();
            pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
            for (k, (pm, pi)) in pairs.into_iter().enumerate() {
                (m[k], i[k]) = (pm, pi);
            }
        }
    }
    Ok(Chromatogram {
        scan_times_s,
        offsets,
        mz,
        intensity,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use netcdf3::{DataSet, FileWriter, Version};

    #[test]
    fn reads_written_andi_file() {
        let path = std::env::temp_dir().join(format!("md-core-andi-{}.cdf", std::process::id()));
        let mut ds = DataSet::new();
        ds.add_fixed_dim("scan_number", 2).unwrap();
        ds.add_fixed_dim("point_number", 3).unwrap();
        ds.add_var_f64("scan_acquisition_time", &["scan_number"])
            .unwrap();
        ds.add_var_i32("scan_index", &["scan_number"]).unwrap();
        ds.add_var_i32("point_count", &["scan_number"]).unwrap();
        ds.add_var_f32("mass_values", &["point_number"]).unwrap();
        ds.add_var_i32("intensity_values", &["point_number"])
            .unwrap();
        let mut w = FileWriter::create_new(&path).unwrap();
        w.set_def(&ds, Version::Classic, 0).unwrap();
        w.write_var_f64("scan_acquisition_time", &[300.0, 300.5])
            .unwrap();
        w.write_var_i32("scan_index", &[0, 2]).unwrap();
        w.write_var_i32("point_count", &[2, 1]).unwrap();
        w.write_var_f32("mass_values", &[73.0, 147.0, 73.0])
            .unwrap();
        w.write_var_i32("intensity_values", &[100, 50, 80]).unwrap();
        w.close().unwrap();

        let c = read(&path).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(c.offsets, vec![0, 2, 3]);
        assert_eq!(c.tic(), vec![150.0, 80.0]);
        assert_eq!(c.spectrum(1).0, &[73.0]);
    }

    #[test]
    fn written_file_reads_back() {
        let c = build(
            vec![1.5, 2.0],
            &[0.0, 2.0],
            &[2.0, 1.0],
            vec![73.0, 147.5, 73.0],
            vec![10.0, 0.5, 20.0],
        )
        .unwrap();
        let path = std::env::temp_dir().join(format!("md-core-andi-w-{}.cdf", std::process::id()));
        write(&path, &c).unwrap();
        let back = read(&path);
        std::fs::remove_file(&path).ok();
        assert_eq!(back.unwrap(), c);
    }

    #[test]
    fn write_refuses_empty_run_and_keeps_target() {
        let path = std::env::temp_dir().join(format!("md-core-andi-e-{}.cdf", std::process::id()));
        std::fs::write(&path, b"keep").unwrap();
        let empty = build(vec![1.0], &[0.0], &[0.0], vec![], vec![]).unwrap();
        assert!(matches!(
            write(&path, &empty),
            Err(AndiError::Inconsistent(_))
        ));
        assert_eq!(std::fs::read(&path).unwrap(), b"keep");
        std::fs::remove_file(&path).ok();
    }

    #[test]
    fn rejects_scan_pointing_past_data() {
        let err = build(vec![1.0], &[0.0], &[5.0], vec![1.0], vec![1.0]).unwrap_err();
        assert!(matches!(err, AndiError::Inconsistent(_)));
    }

    #[test]
    fn sorts_unsorted_scans() {
        let c = build(
            vec![1.0],
            &[0.0],
            &[3.0],
            vec![147.0, 73.0, 100.0],
            vec![1.0, 2.0, 3.0],
        )
        .unwrap();
        assert_eq!(
            c.spectrum(0),
            (&[73.0, 100.0, 147.0][..], &[2.0, 3.0, 1.0][..])
        );
    }
}
