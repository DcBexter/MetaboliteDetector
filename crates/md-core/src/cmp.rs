//! Saved components of a run, the old `.cmp` file (`Compound<int,float>::toDisk`): read and write.
//! Layout: `docs/reverse-engineering/file-formats.md` §3. Only version 4 (MetaboliteDetector 3.3)
//! is decoded; the writer produces the same payload, so the old program opens our files.

use std::io::{Read, Write as _};
use std::path::Path;

use crate::deconv::{Component, FittedPeak};
use crate::library::{
    LibraryCompound, LibraryError, LibraryModel, LibraryPeak, Reader, gunzip, read_compound_record,
};

impl LibraryCompound {
    /// The component as the deconvolution returns it; None without a model peak (a library compound
    /// before version 4). Scan range and model come from the model peak; `FittedPeak::peak` indexes
    /// [`Self::peaks`].
    pub fn component(&self) -> Option<Component> {
        let model = self.model.as_ref()?;
        let m = &model.peak;
        Some(Component {
            first_scan: m.first_scan,
            apex_scan: m.max_scan,
            last_scan: m.last_scan,
            model: model.values.clone(),
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
        })
    }
}

/// Components and the ion mass of each of their peaks (`FittedPeak::peak` indexes it).
pub type SavedComponents = (Vec<Component>, Vec<u32>);

/// Components saved next to a run (`<run>.cmp`), ordered by apex like `deconvolute`'s callers use
/// them, and the ion mass of every peak their `FittedPeak::peak` points to. None without a file; an
/// error when a component lies outside the run or its apex time does not match the run's
/// `scan_times_s` (a `.cmp` of another run).
pub fn saved_components(
    run: &Path,
    scan_times_s: &[f64],
) -> Result<Option<SavedComponents>, LibraryError> {
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
    // Version 4 always has a model peak, so every compound gives a component.
    let mut components: Vec<_> = compounds.iter().filter_map(|c| c.component()).collect();
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
    // One peak list for the run, like `deconvolute` gets it.
    let mut masses = Vec::new();
    for (k, c) in components.iter_mut().zip(&compounds) {
        k.peaks.iter_mut().for_each(|f| f.peak += masses.len());
        masses.extend(c.peaks.iter().map(|p| p.ion.round() as u32));
    }
    components.sort_by_key(|k| k.apex_scan);
    Ok(Some((components, masses)))
}

const VERSION: u32 = 4;
/// Size of a compound whose lists and strings are all empty: a lower bound per compound.
const MIN_COMPOUND_BYTES: usize = 36 + 4 + (66 + 8 + 4 + 4) + 52;
// shortcut: whole file in RAM, capped like `.lbr`; a run has a few thousand components at most.
const MAX_CMP_BYTES: u64 = 1 << 30;

/// Reads a `.cmp` file (gzip, or an already decompressed payload).
pub fn read_cmp(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
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

/// Parses `.cmp` bytes; every length and check value is validated. The library-only fields of the
/// compounds keep the defaults of [`LibraryCompound::new`].
pub fn parse_cmp(bytes: &[u8]) -> Result<Vec<LibraryCompound>, LibraryError> {
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
    let compounds = r.list(MIN_COMPOUND_BYTES, |r| read_compound_record(r, VERSION))?;
    if r.remaining() != 0 {
        return r.corrupt(format!("{} bytes after the last compound", r.remaining()));
    }
    Ok(compounds)
}

fn model_check(m: &LibraryModel) -> u32 {
    let p = &m.peak;
    (m.values.len() as u32)
        .wrapping_add(p.first_scan)
        .wrapping_add(p.last_scan)
        .wrapping_add(p.max_scan)
}

fn compound_check(c: &LibraryCompound) -> u32 {
    (c.quant_ions.len() + c.features.len() + c.name.len() + c.peaks.len()) as u32
}

/// Writes a `.cmp` file (gzip) the old program can open.
pub fn write_cmp(path: impl AsRef<Path>, compounds: &[LibraryCompound]) -> std::io::Result<()> {
    let mut gz =
        flate2::write::GzEncoder::new(std::fs::File::create(path)?, flate2::Compression::default());
    gz.write_all(&cmp_payload(compounds))?;
    gz.finish()?.sync_all()
}

/// The uncompressed `.cmp` payload: the compound records only (library-only fields are not part of
/// a `.cmp`). Names are written as UTF-8 (the old files carry no encoding). The old format always
/// has a model peak; a compound without one gets an empty model (all fields 0).
pub fn cmp_payload(compounds: &[LibraryCompound]) -> Vec<u8> {
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

fn write_peak(w: &mut Vec<u8>, p: &LibraryPeak, short_ion: bool) {
    if short_ion {
        put(w, (p.ion as i16).to_le_bytes());
    } else {
        put(w, p.ion.to_le_bytes());
    }
    for x in [p.sharpness, p.integral, p.quality, p.signal_to_noise] {
        put(w, x.to_le_bytes());
    }
    for x in [p.count, p.first_scan, p.last_scan] {
        put(w, x.to_le_bytes());
    }
    put(w, p.model_parameter.to_le_bytes());
    put(w, p.max_scan.to_le_bytes());
    put(w, p.interpolated_max_scan.to_le_bytes());
}

fn write_compound(w: &mut Vec<u8>, c: &LibraryCompound) {
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
        write_peak(w, p, false);
    }
    let empty;
    let model = match &c.model {
        Some(m) => m,
        None => {
            empty = LibraryModel {
                peak: zero_peak(),
                noise_factor: 0.0,
                values: Vec::new(),
            };
            &empty
        }
    };
    write_peak(w, &model.peak, true);
    put(w, model.noise_factor.to_le_bytes());
    count(w, model.values.len());
    model.values.iter().for_each(|x| put(w, x.to_le_bytes()));
    put(w, model_check(model).to_le_bytes());

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

fn zero_peak() -> LibraryPeak {
    LibraryPeak {
        ion: 0.0,
        sharpness: 0.0,
        integral: 0.0,
        quality: 0.0,
        signal_to_noise: 0.0,
        count: 0,
        first_scan: 0,
        last_scan: 0,
        model_parameter: 0.0,
        max_scan: 0,
        interpolated_max_scan: 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::library::LibrarySpectrum;

    fn peak(ion: f32, first: u32, apex: u32, last: u32) -> LibraryPeak {
        LibraryPeak {
            ion,
            sharpness: 8.5,
            integral: 1234.5,
            quality: 0.25,
            signal_to_noise: 42.0,
            count: 1,
            first_scan: first,
            last_scan: last,
            model_parameter: 0.5,
            max_scan: apex,
            interpolated_max_scan: apex as f64 + 0.25,
        }
    }

    fn sample() -> Vec<LibraryCompound> {
        let c = LibraryCompound {
            spectrum: LibrarySpectrum {
                scan: 12,
                rt_ms: 247160.0,
                ri: -1.0,
                bin_width: 1.0,
                masses: vec![73.0, 147.0],
                values: vec![900, 120],
            },
            peaks: vec![peak(73.0, 10, 12, 15), peak(147.0, 11, 12, 14)],
            model: Some(LibraryModel {
                peak: peak(1.0, 10, 12, 15),
                noise_factor: 0.5,
                values: vec![1.0, 4.0, 9.0, 4.0, 2.0, 1.0],
            }),
            name: "Glycine (3TMS) ä".into(),
            fwhm: None,
            integrals: Some(vec![1.5, 2.5]),
            quant_ions: vec![73.0],
            auto_quant_ions: vec![],
            features: [("k".to_string(), "v".to_string())].into(),
            drt: 0.0,
            time_stamp: 1_600_000_000,
            id: 7,
            ..LibraryCompound::new("")
        };
        vec![
            c.clone(),
            LibraryCompound {
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

        let k = back[0].component().unwrap();
        assert_eq!((k.first_scan, k.apex_scan, k.last_scan), (10, 12, 15));
        assert_eq!(k.spectrum, vec![(73, 900), (147, 120)]);
        assert_eq!((k.peaks[1].peak, k.peaks[1].integral), (1, 1234.5));

        // A library compound without model peak is written with an empty one.
        let bare = LibraryCompound::new("x");
        assert!(bare.component().is_none());
        let back = parse_cmp(&cmp_payload(std::slice::from_ref(&bare))).unwrap();
        assert_eq!(back[0].model.as_ref().unwrap().values, Vec::<f64>::new());
        assert_eq!(back[0].name, "x");
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
        let (k, masses) = saved_components(&run, &times).unwrap().unwrap();
        assert_eq!(k.len(), 2);
        // sample(): the first compound has peaks at 73 and 147, the second none.
        assert_eq!(masses, [73, 147]);
        let fitted: Vec<_> = k
            .iter()
            .flat_map(|k| &k.peaks)
            .map(|f| masses[f.peak])
            .collect();
        assert_eq!(fitted, [73, 147]);
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
