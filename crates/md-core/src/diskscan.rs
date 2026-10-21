//! Old MetaboliteDetector 3.x raw-data cache: `<run>.bin` (scan records) + `<run>.idx` (index).
//! The old app converted every imported NetCDF into this pair (`GCMSDiskScanWriter`); some users
//! only kept these files. Layout: docs/reverse-engineering/file-formats.md, section 2.

use std::path::{Path, PathBuf};

use crate::Chromatogram;
use crate::andi::{self, AndiError};
use crate::ri::StoredCalibration;
use crate::tof::TofCalibration;

#[derive(Debug, thiserror::Error)]
pub enum DiskScanError {
    #[error("cannot read {0}: {1}")]
    Io(PathBuf, std::io::Error),
    #[error("not a MetaboliteDetector index file (magic {0:#x})")]
    Magic(i64),
    #[error("corrupt index or data file: {0}")]
    Corrupt(String),
}

/// `DataProvider::readIndexFile` refuses any other value.
const MAGIC: i64 = 0x4b49_4104;
/// Record header in the `.bin`: size i32, masses_ptr i64, tic f64, scan u32, rt f64, n u32.
const RECORD_HEADER: usize = 36;

struct Entry {
    scan: u32,
    pos: u64,
    masses_ptr: u64,
    rt_ms: f64,
    n: usize,
    data_type: u32,
    /// `(int)` of the scan's `BaselineStatistics` mean (non-zero only in type 1 caches).
    mean: f32,
}

/// Reads a run from its `.bin` and `.idx` files; `path` may name either (or the bare stem).
// shortcut: both files in RAM like `andi::read`; seek per scan if old caches outgrow memory.
pub fn read(path: impl AsRef<Path>) -> Result<Chromatogram, DiskScanError> {
    let path = path.as_ref();
    let load = |ext| {
        let p = path.with_extension(ext);
        std::fs::read(&p).map_err(|e| DiskScanError::Io(p, e))
    };
    parse(&load("idx")?, &load("bin")?)
}

fn parse(idx: &[u8], bin: &[u8]) -> Result<Chromatogram, DiskScanError> {
    let (_, tof, mut entries) = read_index(idx)?;
    entries.sort_by_key(|e| e.scan);
    if let Some(w) = entries.windows(2).find(|w| w[0].scan == w[1].scan) {
        return corrupt(format!("scan {} listed twice", w[0].scan));
    }
    let total: usize = entries.iter().map(|e| e.n).sum();
    // Each point takes at least 4 bytes in the .bin; caps allocations on a corrupt index.
    if total > bin.len() / 4 {
        return corrupt(format!(
            "index lists {total} points, data file has {} bytes",
            bin.len()
        ));
    }
    let (mut mz, mut intensity) = (Vec::with_capacity(total), Vec::with_capacity(total));
    let (mut index, mut count) = (
        Vec::with_capacity(entries.len()),
        Vec::with_capacity(entries.len()),
    );
    for e in &entries {
        let mut r = Reader::at(bin, e.pos)?;
        let _size = r.u32()?;
        let masses_ptr = r.u64()?;
        let _tic = r.f64()?;
        let scan = r.u32()?;
        let _rt = r.f64()?;
        let n = r.u32()? as usize;
        if scan != e.scan || n != e.n || masses_ptr != e.masses_ptr {
            return corrupt(format!(
                "record of scan {} does not match the index",
                e.scan
            ));
        }
        index.push(mz.len() as f64);
        count.push(n as f64);
        let values = Reader::at(bin, e.pos + RECORD_HEADER as u64)?.take(4 * n)?;
        let masses = Reader::at(bin, e.masses_ptr)?.take(4 * n)?;
        let (values, masses) = (values.as_chunks::<4>().0, masses.as_chunks::<4>().0);
        intensity.extend(values.iter().map(|b| i32::from_le_bytes(*b) as f32));
        // A cache of FastFlight data stores TOF bins; the old reader maps them through the polynomial.
        let masses = masses.iter().map(|b| f32::from_le_bytes(*b));
        match &tof {
            Some(tof) => mz.extend(masses.map(|bin| tof.mz(bin))),
            None => mz.extend(masses),
        }
    }
    let times = entries.iter().map(|e| e.rt_ms / 1000.0).collect();
    let mut run = andi::build(times, &index, &count, mz, intensity).map_err(|e| match e {
        AndiError::Inconsistent(s) => DiskScanError::Corrupt(s),
        other => DiskScanError::Corrupt(other.to_string()),
    })?;
    // Only the CDF import writes DataType 2 (nominal masses); FastFlight writes 0 or 1.
    run.exact_masses = entries.iter().any(|e| e.data_type != 2);
    if entries.iter().any(|e| e.mean != 0.0) {
        run.baseline.means = entries.iter().map(|e| e.mean).collect();
    }
    Ok(run)
}

/// The RI calibration stored in a run's `.idx` (the old RI wizard writes it there); `path` may name
/// the `.bin`, the `.idx` or the bare stem.
pub fn read_calibration(path: impl AsRef<Path>) -> Result<StoredCalibration, DiskScanError> {
    let p = path.as_ref().with_extension("idx");
    let idx = std::fs::read(&p).map_err(|e| DiskScanError::Io(p, e))?;
    Ok(read_index(&idx)?.0)
}

/// `DataProvider::readIndexFile` (0x75bf30), in its read order.
fn read_index(
    idx: &[u8],
) -> Result<(StoredCalibration, Option<TofCalibration>, Vec<Entry>), DiskScanError> {
    let mut r = Reader { buf: idx, pos: 0 };
    let magic = r.u64()? as i64;
    if magic != MAGIC {
        return Err(DiskScanError::Magic(magic));
    }
    // sizeof(int), sizeof(double), sizeof(long long), scan count, sizeof(mass), sizeof(value).
    let sizes = [r.u32()?, r.u32()?, r.u32()?];
    let scans = r.u32()? as usize;
    let (mass_size, value_size) = (r.u32()?, r.u32()?);
    if sizes != [4, 8, 8] || mass_size != 4 || value_size != 4 {
        return corrupt(format!(
            "unsupported type sizes {sizes:?}, {mass_size}, {value_size}"
        ));
    }
    for _ in 0..2 * u64::from(r.u32()?) {
        let len = r.u32()? as usize; // features: key and value strings
        r.take(len)?;
    }
    // TOF calibration polynomial (f64 coefficients); CdfImporter writes none.
    let terms = r.u32()? as usize;
    let tof = match terms {
        0 => None,
        n if n <= 64 => {
            let coefficients = (0..n).map(|_| r.f64()).collect::<Result<_, _>>()?;
            Some(TofCalibration { coefficients })
        }
        n => return corrupt(format!("TOF polynomial with {n} terms")),
    };
    r.take(5 * 8)?; // global BaselineStatistics: 4 x f64 + i64
    let (calibration, len) = StoredCalibration::parse(&r.buf[r.pos..])
        .map_err(|e| DiskScanError::Corrupt(e.to_string()))?;
    r.take(len)?;
    let ions = r.u32()? as usize; // all distinct masses, f32
    r.take(ions.saturating_mul(4))?;
    // 84 bytes per scan; checked before allocating from the untrusted count.
    if scans > r.remaining() / 84 {
        return corrupt(format!("{scans} scans do not fit into the index"));
    }
    (0..scans)
        .map(|_| {
            let (scan, pos, masses_ptr) = (r.u32()?, r.u64()?, r.u64()?);
            let _tic = r.f64()?;
            let (rt_ms, n, data_type) = (r.f64()?, r.u32()? as usize, r.u32()?);
            // Per-scan BaselineStatistics (mean, sd, msd, nf, n); `getSignal` reads `(int) mean`.
            let mean = r.f64()?;
            // A real baseline is far smaller; this keeps the correction exact in f32 (below 2^24).
            if mean.is_nan() || mean.abs() >= 16_777_216.0 {
                return corrupt(format!("scan {scan} has baseline mean {mean}"));
            }
            let mean = mean as i32 as f32;
            r.take(4 * 8)?;
            Ok(Entry {
                scan,
                pos,
                masses_ptr,
                rt_ms,
                n,
                data_type,
                mean,
            })
        })
        .collect::<Result<_, _>>()
        .map(|entries| (calibration, tof, entries))
}

fn corrupt<T>(msg: String) -> Result<T, DiskScanError> {
    Err(DiskScanError::Corrupt(msg))
}

/// Bounds-checked little-endian reads from untrusted bytes.
struct Reader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn at(buf: &'a [u8], pos: u64) -> Result<Self, DiskScanError> {
        match usize::try_from(pos) {
            Ok(pos) if pos <= buf.len() => Ok(Self { buf, pos }),
            _ => corrupt(format!("offset {pos} is past the end of the data file")),
        }
    }

    fn remaining(&self) -> usize {
        self.buf.len() - self.pos
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], DiskScanError> {
        if n > self.remaining() {
            return corrupt(format!(
                "file ends at byte {}, needed {n} more",
                self.buf.len()
            ));
        }
        self.pos += n;
        Ok(&self.buf[self.pos - n..self.pos])
    }

    fn u32(&mut self) -> Result<u32, DiskScanError> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn u64(&mut self) -> Result<u64, DiskScanError> {
        Ok(u64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    fn f64(&mut self) -> Result<f64, DiskScanError> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use super::*;

    /// What the old RI wizard stores in the index of a calibrated run.
    fn wizard_calibration() -> StoredCalibration {
        StoredCalibration {
            points_ms: vec![(1500.0, 1000.0), (2000.0, 1100.0)],
            drt_ms: 100.0,
            source: b"C10-C11".to_vec(),
        }
    }

    /// Index and data file as `GCMSDiskScanWriter` lays them out, scans given as (rt_ms, m/z, values).
    pub(crate) fn files(scans: &[(f64, &[f32], &[i32])]) -> (Vec<u8>, Vec<u8>) {
        files_with(scans, &[])
    }

    /// Same, with the TOF polynomial of a FastFlight cache.
    fn files_with(scans: &[(f64, &[f32], &[i32])], tof: &[f64]) -> (Vec<u8>, Vec<u8>) {
        let (mut idx, mut bin) = (Vec::new(), Vec::new());
        let u32s = |v: &mut Vec<u8>, xs: &[u32]| xs.iter().for_each(|x| v.extend(x.to_le_bytes()));
        idx.extend(MAGIC.to_le_bytes());
        u32s(&mut idx, &[4, 8, 8, scans.len() as u32, 4, 4]);
        u32s(&mut idx, &[1, 1]); // one feature: key "k", value "v"
        idx.push(b'k');
        u32s(&mut idx, &[1]);
        idx.push(b'v');
        u32s(&mut idx, &[tof.len() as u32]);
        tof.iter().for_each(|c| idx.extend(c.to_le_bytes()));
        idx.extend([0; 40]); // global baseline statistics
        idx.extend(wizard_calibration().to_bytes());
        u32s(&mut idx, &[1]); // all ions
        idx.extend(73f32.to_le_bytes());
        for (scan, (rt, mz, values)) in scans.iter().enumerate() {
            let (pos, n) = (bin.len() as u64, mz.len());
            let masses_ptr = pos + (RECORD_HEADER + 4 * n) as u64;
            bin.extend(((RECORD_HEADER + 8 * n) as i32).to_le_bytes());
            bin.extend(masses_ptr.to_le_bytes());
            bin.extend(0f64.to_le_bytes());
            u32s(&mut bin, &[scan as u32]);
            bin.extend(rt.to_le_bytes());
            u32s(&mut bin, &[n as u32]);
            values.iter().for_each(|v| bin.extend(v.to_le_bytes()));
            mz.iter().for_each(|m| bin.extend(m.to_le_bytes()));
            u32s(&mut idx, &[scan as u32]);
            idx.extend(pos.to_le_bytes());
            idx.extend(masses_ptr.to_le_bytes());
            idx.extend(0f64.to_le_bytes());
            idx.extend(rt.to_le_bytes());
            u32s(&mut idx, &[n as u32, 2]);
            idx.extend([0; 40]);
        }
        (idx, bin)
    }

    #[test]
    fn reads_scans() {
        let (idx, bin) = files(&[(1500.0, &[73.0, 147.0], &[10, 5]), (2000.0, &[73.0], &[7])]);
        let c = parse(&idx, &bin).unwrap();
        assert_eq!(c.scan_times_s, vec![1.5, 2.0]);
        assert_eq!(c.offsets, vec![0, 2, 3]);
        assert_eq!(c.spectrum(0), (&[73.0, 147.0][..], &[10.0, 5.0][..]));
        assert!(!c.exact_masses, "DataType 2: the CDF import");
        // A FastFlight cache (DataType 1) keeps exact masses.
        let mut ff = idx.clone();
        let entry = idx.len() - 84;
        ff[entry + 40..entry + 44].copy_from_slice(&1u32.to_le_bytes());
        assert!(parse(&ff, &bin).unwrap().exact_masses);
        let stored = read_index(&idx).unwrap().0;
        assert_eq!(stored, wizard_calibration());
        assert_eq!(stored.calibration().unwrap().ri(1.85), 1050.0);
    }

    #[test]
    fn type_1_caches_fill_and_correct_with_the_scan_baseline() {
        let (mut idx, bin) = files(&[(1500.0, &[73.0, 147.0], &[10, 5]), (2000.0, &[73.0], &[7])]);
        let first = idx.len() - 2 * 84;
        for (s, mean) in [2.7f64, 5.0].into_iter().enumerate() {
            let e = first + 84 * s;
            idx[e + 40..e + 44].copy_from_slice(&1u32.to_le_bytes()); // DataType 1
            idx[e + 44..e + 52].copy_from_slice(&mean.to_le_bytes());
        }
        let mut c = parse(&idx, &bin).unwrap();
        assert_eq!(c.baseline.means, vec![2.0, 5.0], "(int) of the mean");
        let tic = c.tic();
        // Without correction a missing mass reads the scan's mean (old `getSignal`).
        let t = crate::peaks::IonTraces::new(&c, 0.5).unwrap();
        assert_eq!(
            (t.trace(0), t.trace(1)),
            (&[10.0, 7.0][..], &[5.0, 5.0][..])
        );
        c.set_baseline_correction(true);
        assert_eq!(c.spectrum(0).1, &[8.0, 3.0][..]);
        let t = crate::peaks::IonTraces::new(&c, 0.5).unwrap();
        assert_eq!((t.trace(0), t.trace(1)), (&[8.0, 2.0][..], &[3.0, 0.0][..]));
        assert_eq!(c.tic(), tic, "the TIC stays uncorrected");
        c.set_baseline_correction(false);
        assert_eq!(c, parse(&idx, &bin).unwrap());
        // CDF caches store zero means: nothing to keep, nothing to correct.
        let (idx, bin) = files(&[(1500.0, &[73.0], &[10])]);
        let mut cdf = parse(&idx, &bin).unwrap();
        assert!(cdf.baseline.means.is_empty());
        cdf.set_baseline_correction(true);
        assert_eq!(cdf.spectrum(0).1, &[10.0][..]);
        // A NaN or huge mean is a corrupt index, not a baseline.
        for bad in [f64::NAN, 1e9] {
            let mut bad_idx = idx.clone();
            let e = bad_idx.len() - 84;
            bad_idx[e + 44..e + 52].copy_from_slice(&bad.to_le_bytes());
            assert!(parse(&bad_idx, &bin).is_err(), "mean {bad}");
        }
    }

    #[test]
    fn tof_caches_map_bins_to_masses() {
        let (idx, bin) = files_with(&[(1500.0, &[100.0, 200.0], &[10, 5])], &[1.0, 0.5]);
        let c = parse(&idx, &bin).unwrap();
        assert_eq!(c.spectrum(0), (&[51.0, 101.0][..], &[10.0, 5.0][..]));
        // An absurd term count must not allocate or read past the index.
        let mut bad = idx.clone();
        let at = 28 + 2 * 4 + 2; // after the type sizes, scan count, mass/value size and the feature
        bad[at..at + 4].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(parse(&bad, &bin).is_err());
    }

    #[test]
    fn corrupt_files_are_errors_not_panics() {
        let (idx, bin) = files(&[(1500.0, &[73.0, 147.0], &[10, 5]), (2000.0, &[73.0], &[7])]);
        for len in 0..idx.len() {
            assert!(parse(&idx[..len], &bin).is_err(), "idx cut at {len}");
        }
        for len in 0..bin.len() {
            assert!(parse(&idx, &bin[..len]).is_err(), "bin cut at {len}");
        }
        let mut bad = idx.clone();
        bad[0] ^= 1;
        assert!(matches!(parse(&bad, &bin), Err(DiskScanError::Magic(_))));
        // Huge scan count, then a record position past the data file (first scan entry + 4).
        let first_entry = idx.len() - 2 * 84;
        for (at, bytes) in [(20, &[0xff; 4][..]), (first_entry + 4, &[0xff; 8][..])] {
            let mut bad = idx.clone();
            bad[at..at + bytes.len()].copy_from_slice(bytes);
            assert!(parse(&bad, &bin).is_err(), "patched byte {at}");
        }
    }
}
