//! Saved components of a run, the old `.cmp` file (`Compound<int,float>::toDisk`): read and write.
//! Layout: `docs/reverse-engineering/file-formats.md` §3. Only version 4 (MetaboliteDetector 3.3)
//! is decoded; the writer produces the same payload, so the old program opens our files.

use std::collections::BTreeMap;
use std::io::{Read, Write as _};
use std::path::Path;

use crate::deconv::{Component, FittedPeak};
use crate::library::{LibraryError, LibrarySpectrum, Reader, gunzip, read_spectrum, text};

/// Ion peak of a component (`PeakEssential<float>`, ion as m/z) or of its model peak
/// (`PeakEssential<short>`, ion as an integer code), fields in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct PeakEssential<I> {
    pub ion: I,
    pub sharpness: f64,
    pub integral: f64,
    pub quality: f64,
    pub signal_to_noise: f64,
    /// Meaning not recovered; 1 in every sample peak.
    pub count: u32,
    pub first_scan: u32,
    pub last_scan: u32,
    /// Share of the model peak in this ion trace (`c2` of the fit), 1 for the model peak itself.
    pub model_parameter: f64,
    pub apex_scan: u32,
    pub precise_apex: f64,
}

/// Model peak of a component (`Peak<double,short>`).
#[derive(Debug, Clone, PartialEq)]
pub struct ModelPeak {
    pub essential: PeakEssential<i16>,
    pub noise_factor: f64,
    /// Model intensity of each scan in `first_scan..=last_scan`.
    pub values: Vec<f64>,
}

/// One saved component (`Compound<int,float>`), fields in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct CmpCompound {
    /// Deconvoluted spectrum; `rt_ms` is the apex time, `ri` is -1 without a calibration.
    pub spectrum: LibrarySpectrum,
    /// Ion peaks behind the spectrum, ascending ion.
    pub peaks: Vec<PeakEssential<f32>>,
    pub model: ModelPeak,
    pub name: String,
    pub fwhm: Option<Vec<f64>>,
    pub integrals: Option<Vec<f64>>,
    pub quant_ions: Vec<f32>,
    pub auto_quant_ions: Vec<f32>,
    pub features: BTreeMap<String, String>,
    pub drt: f64,
    pub time_stamp: i64,
    pub id: i64,
}

impl CmpCompound {
    /// The component as the deconvolution returns it. Scan range and model come from the model
    /// peak; `FittedPeak::peak` indexes [`Self::peaks`].
    pub fn component(&self) -> Component {
        let m = &self.model.essential;
        Component {
            first_scan: m.first_scan,
            apex_scan: m.apex_scan,
            last_scan: m.last_scan,
            model: self.model.values.clone(),
            spectrum: (self.spectrum.masses.iter().zip(&self.spectrum.values))
                .map(|(&mz, &v)| (mz.round() as u32, v.max(0) as u32))
                .collect(),
            peaks: (self.peaks.iter().enumerate())
                .map(|(i, p)| FittedPeak {
                    peak: i,
                    model_share: p.model_parameter,
                    integral: p.integral,
                })
                .collect(),
        }
    }
}

/// Components saved next to a run (`<run>.cmp`), ordered by apex like `deconvolute`'s callers use
/// them. None without a file; an error when a component lies outside the run or its apex time does
/// not match the run's `scan_times_s` (a `.cmp` of another run).
pub fn saved_components(
    run: &Path,
    scan_times_s: &[f64],
) -> Result<Option<Vec<Component>>, LibraryError> {
    let path = run.with_extension("cmp");
    if !path.is_file() {
        return Ok(None);
    }
    let compounds = read_cmp(&path)?;
    if let Some(c) = compounds.iter().find(|c| {
        let s = &c.spectrum;
        // Both come from the same raw data; 1 ms leaves room for rounding only.
        (scan_times_s.get(s.scan as usize)).is_none_or(|t| (t * 1000.0 - s.rt_ms).abs() > 1.0)
    }) {
        return Err(LibraryError::Unsupported(format!(
            "{}: component at scan {} ({} ms) is not from this run",
            path.display(),
            c.spectrum.scan,
            c.spectrum.rt_ms
        )));
    }
    let scans = scan_times_s.len();
    let mut components: Vec<_> = compounds.iter().map(CmpCompound::component).collect();
    if let Some(k) = components.iter().find(|k| {
        k.last_scan as usize >= scans
            || k.first_scan > k.apex_scan
            || k.apex_scan > k.last_scan
            || k.model.len() != (k.last_scan - k.first_scan) as usize + 1
    }) {
        return Err(LibraryError::Unsupported(format!(
            "{}: component at scan {} does not fit a run of {scans} scans",
            path.display(),
            k.apex_scan
        )));
    }
    components.sort_by_key(|k| k.apex_scan);
    Ok(Some(components))
}

const VERSION: u32 = 4;
const PEAK_BYTES: usize = 68;
/// Size of a compound whose lists and strings are all empty: a lower bound per compound.
const MIN_COMPOUND_BYTES: usize = 36 + 4 + (66 + 8 + 4 + 4) + 52;
// shortcut: whole file in RAM, capped like `.lbr`; a run has a few thousand components at most.
const MAX_CMP_BYTES: u64 = 1 << 30;

/// Reads a `.cmp` file (gzip, or an already decompressed payload).
pub fn read_cmp(path: impl AsRef<Path>) -> Result<Vec<CmpCompound>, LibraryError> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)?
        .take(MAX_CMP_BYTES + 1)
        .read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_CMP_BYTES {
        return Err(LibraryError::Unsupported(format!(
            "compound files larger than {MAX_CMP_BYTES} bytes"
        )));
    }
    parse_cmp(&bytes)
}

/// Parses `.cmp` bytes; every length and check value is validated.
pub fn parse_cmp(bytes: &[u8]) -> Result<Vec<CmpCompound>, LibraryError> {
    let unzipped;
    let data = if bytes.starts_with(&[0x1f, 0x8b]) {
        unzipped = gunzip(bytes)?;
        &unzipped[..]
    } else {
        bytes
    };
    let mut r = Reader::new(data);
    match r.u32()? {
        VERSION => {}
        // shortcut: versions 1-3 store Peak<int,float> records whose old layout is unknown;
        // decode them when a file from an older MetaboliteDetector turns up.
        v => {
            return Err(LibraryError::Unsupported(format!(
                "compound file version {v}"
            )));
        }
    }
    let compounds = r.list(MIN_COMPOUND_BYTES, read_compound)?;
    if r.remaining() != 0 {
        return r.corrupt(format!("{} bytes after the last compound", r.remaining()));
    }
    Ok(compounds)
}

fn read_essential<I>(
    r: &mut Reader,
    ion: impl FnOnce(&mut Reader) -> Result<I, LibraryError>,
) -> Result<PeakEssential<I>, LibraryError> {
    Ok(PeakEssential {
        ion: ion(r)?,
        sharpness: r.f64()?,
        integral: r.f64()?,
        quality: r.f64()?,
        signal_to_noise: r.f64()?,
        count: r.u32()?,
        first_scan: r.u32()?,
        last_scan: r.u32()?,
        model_parameter: r.f64()?,
        apex_scan: r.u32()?,
        precise_apex: r.f64()?,
    })
}

fn model_check(m: &ModelPeak) -> u32 {
    let e = &m.essential;
    (m.values.len() as u32)
        .wrapping_add(e.first_scan)
        .wrapping_add(e.last_scan)
        .wrapping_add(e.apex_scan)
}

fn compound_check(c: &CmpCompound) -> u32 {
    (c.quant_ions.len() + c.features.len() + c.name.len() + c.peaks.len()) as u32
}

fn read_compound(r: &mut Reader) -> Result<CmpCompound, LibraryError> {
    let spectrum = read_spectrum(r)?;
    let peaks = r.list(PEAK_BYTES, |r| read_essential(r, |r| r.f32()))?;
    let essential = read_essential(r, |r| r.array().map(i16::from_le_bytes))?;
    let model = ModelPeak {
        essential,
        noise_factor: r.f64()?,
        values: r.list(8, Reader::f64)?,
    };
    if r.u32()? != model_check(&model) {
        return r.corrupt("model peak check value mismatch");
    }
    let name_bytes = r.bytes()?;
    let fwhm = r.f64_list_or_null()?;
    let integrals = r.f64_list_or_null()?;
    let quant_ions = r.list(4, Reader::f32)?;
    let auto_quant_ions = r.list(4, Reader::f32)?;
    let feature_list = r.list(8, |r| Ok((r.string()?, r.string()?)))?;
    let expected = (quant_ions.len() + feature_list.len() + name_bytes.len() + peaks.len()) as u32;
    let c = CmpCompound {
        spectrum,
        peaks,
        model,
        name: text(name_bytes),
        fwhm,
        integrals,
        quant_ions,
        auto_quant_ions,
        features: feature_list.into_iter().collect(),
        drt: r.f64()?,
        time_stamp: r.i64()?,
        id: r.i64()?,
    };
    if r.u32()? != expected {
        return r.corrupt("compound check value mismatch");
    }
    Ok(c)
}

/// Writes a `.cmp` file (gzip) the old program can open.
pub fn write_cmp(path: impl AsRef<Path>, compounds: &[CmpCompound]) -> std::io::Result<()> {
    let mut gz =
        flate2::write::GzEncoder::new(std::fs::File::create(path)?, flate2::Compression::default());
    gz.write_all(&cmp_payload(compounds))?;
    gz.finish()?.sync_all()
}

/// The uncompressed `.cmp` payload. Names are written as UTF-8 (the old files carry no encoding).
pub fn cmp_payload(compounds: &[CmpCompound]) -> Vec<u8> {
    let mut w = Vec::new();
    put(&mut w, VERSION.to_le_bytes());
    count(&mut w, compounds.len());
    for c in compounds {
        write_compound(&mut w, c);
    }
    w
}

fn put<const N: usize>(w: &mut Vec<u8>, b: [u8; N]) {
    w.extend_from_slice(&b);
}

fn count(w: &mut Vec<u8>, n: usize) {
    put(w, (n as u32).to_le_bytes());
}

fn string(w: &mut Vec<u8>, s: &str) {
    count(w, s.len());
    w.extend_from_slice(s.as_bytes());
}

fn f64_list_or_null(w: &mut Vec<u8>, v: &Option<Vec<f64>>) {
    match v {
        None => put(w, (-1i32).to_le_bytes()),
        Some(v) => {
            count(w, v.len());
            v.iter().for_each(|x| put(w, x.to_le_bytes()));
        }
    }
}

fn write_essential<I: Copy>(w: &mut Vec<u8>, p: &PeakEssential<I>, ion: fn(I) -> Vec<u8>) {
    w.extend(ion(p.ion));
    for x in [p.sharpness, p.integral, p.quality, p.signal_to_noise] {
        put(w, x.to_le_bytes());
    }
    for x in [p.count, p.first_scan, p.last_scan] {
        put(w, x.to_le_bytes());
    }
    put(w, p.model_parameter.to_le_bytes());
    put(w, p.apex_scan.to_le_bytes());
    put(w, p.precise_apex.to_le_bytes());
}

fn write_compound(w: &mut Vec<u8>, c: &CmpCompound) {
    let s = &c.spectrum;
    put(w, s.scan.to_le_bytes());
    put(w, s.rt_ms.to_le_bytes());
    put(w, s.ri.to_le_bytes());
    put(w, s.bin_width.to_le_bytes());
    count(w, s.masses.len());
    s.masses.iter().for_each(|x| put(w, x.to_le_bytes()));
    count(w, s.values.len());
    s.values.iter().for_each(|x| put(w, x.to_le_bytes()));
    put(w, s.scan.wrapping_add(s.masses.len() as u32).to_le_bytes());

    count(w, c.peaks.len());
    for p in &c.peaks {
        write_essential(w, p, |i| i.to_le_bytes().to_vec());
    }
    write_essential(w, &c.model.essential, |i| i.to_le_bytes().to_vec());
    put(w, c.model.noise_factor.to_le_bytes());
    count(w, c.model.values.len());
    c.model.values.iter().for_each(|x| put(w, x.to_le_bytes()));
    put(w, model_check(&c.model).to_le_bytes());

    string(w, &c.name);
    f64_list_or_null(w, &c.fwhm);
    f64_list_or_null(w, &c.integrals);
    for set in [&c.quant_ions, &c.auto_quant_ions] {
        count(w, set.len());
        set.iter().for_each(|x| put(w, x.to_le_bytes()));
    }
    count(w, c.features.len());
    for (k, v) in &c.features {
        string(w, k);
        string(w, v);
    }
    put(w, c.drt.to_le_bytes());
    put(w, c.time_stamp.to_le_bytes());
    put(w, c.id.to_le_bytes());
    put(w, compound_check(c).to_le_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn peak<I>(ion: I, first: u32, apex: u32, last: u32) -> PeakEssential<I> {
        PeakEssential {
            ion,
            sharpness: 8.5,
            integral: 1234.5,
            quality: 0.25,
            signal_to_noise: 42.0,
            count: 1,
            first_scan: first,
            last_scan: last,
            model_parameter: 0.5,
            apex_scan: apex,
            precise_apex: apex as f64 + 0.25,
        }
    }

    fn sample() -> Vec<CmpCompound> {
        let c = CmpCompound {
            spectrum: LibrarySpectrum {
                scan: 12,
                rt_ms: 247160.0,
                ri: -1.0,
                bin_width: 1.0,
                masses: vec![73.0, 147.0],
                values: vec![900, 120],
            },
            peaks: vec![peak(73.0, 10, 12, 15), peak(147.0, 11, 12, 14)],
            model: ModelPeak {
                essential: peak(1, 10, 12, 15),
                noise_factor: 0.5,
                values: vec![1.0, 4.0, 9.0, 4.0, 2.0, 1.0],
            },
            name: "Glycine (3TMS) ä".into(),
            fwhm: None,
            integrals: Some(vec![1.5, 2.5]),
            quant_ions: vec![73.0],
            auto_quant_ions: vec![],
            features: [("k".to_string(), "v".to_string())].into(),
            drt: 0.0,
            time_stamp: 1_600_000_000,
            id: 7,
        };
        vec![
            c.clone(),
            CmpCompound {
                name: String::new(),
                peaks: vec![],
                ..c
            },
        ]
    }

    #[test]
    fn round_trip_and_component() {
        let file = std::env::temp_dir().join(format!("md-core-cmp-{}.cmp", std::process::id()));
        write_cmp(&file, &sample()).unwrap();
        let back = read_cmp(&file).unwrap();
        std::fs::remove_file(&file).ok();
        assert_eq!(back, sample());
        assert_eq!(cmp_payload(&back), cmp_payload(&sample()));

        let k = back[0].component();
        assert_eq!((k.first_scan, k.apex_scan, k.last_scan), (10, 12, 15));
        assert_eq!(k.spectrum, vec![(73, 900), (147, 120)]);
        assert_eq!((k.peaks[1].peak, k.peaks[1].integral), (1, 1234.5));
    }

    #[test]
    fn saved_components_next_to_a_run() {
        let dir = std::env::temp_dir().join(format!("md-core-saved-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let run = dir.join("run.cdf");
        // sample(): apex scan 12 at 247.16 s.
        let times: Vec<f64> = (0..20).map(|i| 247.16 + (i - 12) as f64 * 0.2).collect();
        assert!(saved_components(&run, &times).unwrap().is_none());
        write_cmp(dir.join("run.cmp"), &sample()).unwrap();
        let k = saved_components(&run, &times).unwrap().unwrap();
        assert_eq!(k.len(), 2);
        assert!(
            saved_components(&run, &times[..15]).is_err(),
            "last scan 15 needs 16 scans"
        );
        let other: Vec<f64> = times.iter().map(|t| t + 0.1).collect();
        assert!(saved_components(&run, &other).is_err(), "another run");
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn broken_files_are_errors_not_panics() {
        let good = cmp_payload(&sample());
        for n in 0..good.len() {
            assert!(parse_cmp(&good[..n]).is_err(), "truncated at {n}");
        }
        let mut v3 = good.clone();
        v3[0] = 3;
        assert!(matches!(parse_cmp(&v3), Err(LibraryError::Unsupported(_))));
        let mut bad = good.clone();
        let last = bad.len() - 1;
        bad[last] ^= 1; // compound check value
        assert!(parse_cmp(&bad).is_err());
        let mut huge = good;
        huge[4..8].copy_from_slice(&u32::MAX.to_le_bytes());
        assert!(parse_cmp(&huge).is_err());
    }
}
