//! Saved components of a run, the old `.cmp` file (`Compound<int,float>::toDisk`): read and write.
//! Layout: `docs/reverse-engineering/file-formats.md` §3. Only version 4 (MetaboliteDetector 3.3)
//! is decoded; the writer produces the same payload, so the old program opens our files.

use std::io::{Read, Write as _};
use std::path::Path;

use crate::deconv::{Component, FittedPeak};
use crate::library::{
    LibraryCompound, LibraryError, LibraryModel, LibraryPeak, LibrarySpectrum, Reader, gunzip,
    read_compound_record,
};
use crate::peaks::Peak;

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
            exact_spectrum: (self.spectrum.masses.iter().zip(&self.spectrum.values))
                .map(|(&mz, &v)| (mz, v.max(0) as u32))
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

/// Whether `path` names a `.cmp` file.
pub fn is_cmp(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("cmp"))
}

/// Components, the ion mass of each peak and the peaks themselves (`FittedPeak::peak` indexes both
/// lists), as `settings::analyze` returns them. A saved peak has no ion-trace row (`Peak::row` is 0).
pub type SavedComponents = (Vec<Component>, Vec<u32>, Vec<Peak>);

/// Components saved next to a run (`<run>.cmp`), ordered by apex like `deconvolute`'s callers use
/// them, with the ion mass and the peak of every `FittedPeak::peak`. None without a file; an
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
    components_from_file(&path, scan_times_s).map(Some)
}

/// The components of the `.cmp` file `path` for the run with `scan_times_s`; see [`saved_components`].
pub fn components_from_file(
    path: &Path,
    scan_times_s: &[f64],
) -> Result<SavedComponents, LibraryError> {
    let compounds = read_cmp(path)?;
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
    components_of(path, &compounds, scan_times_s.len())
}

/// The components of a `.cmp` file without its run (batch input), and scan times that are right at
/// the apex scans only (`rt_ms / 1000`, every other scan is 0): all the batch passes read of them.
pub fn standalone_components(path: &Path) -> Result<(SavedComponents, Vec<f64>), LibraryError> {
    let compounds = read_cmp(path)?;
    let scans = (compounds.iter().filter_map(|c| c.model.as_ref()))
        .map(|m| m.peak.last_scan as usize + 1)
        .max()
        .unwrap_or(0);
    // A scan number far beyond any run would allocate for nothing.
    if scans > MAX_SCANS {
        return Err(LibraryError::Unsupported(format!(
            "{}: scan number {scans} is not plausible",
            path.display()
        )));
    }
    let mut times = vec![0.0; scans];
    for c in &compounds {
        let apex = c.model.as_ref().map(|m| m.peak.max_scan as usize);
        if let Some(t) = apex.and_then(|a| times.get_mut(a)) {
            *t = c.spectrum.rt_ms / 1000.0;
        }
    }
    Ok((components_of(path, &compounds, scans)?, times))
}

const MAX_SCANS: usize = 1 << 22;

fn components_of(
    path: &Path,
    compounds: &[LibraryCompound],
    scans: usize,
) -> Result<SavedComponents, LibraryError> {
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
    let (mut masses, mut peaks) = (Vec::new(), Vec::new());
    for (k, c) in components.iter_mut().zip(compounds) {
        k.peaks.iter_mut().for_each(|f| f.peak += masses.len());
        masses.extend(c.peaks.iter().map(|p| p.ion.round() as u32));
        peaks.extend(c.peaks.iter().map(peak_of));
    }
    components.sort_by_key(|k| k.apex_scan);
    Ok((components, masses, peaks))
}

fn peak_of(p: &LibraryPeak) -> Peak {
    Peak {
        row: 0,
        first_scan: p.first_scan,
        apex_scan: p.max_scan,
        last_scan: p.last_scan,
        precise_apex: Some(p.interpolated_max_scan as f32),
        signal_to_noise: p.signal_to_noise as f32,
        quality: p.quality as f32,
        sharpness: p.sharpness as f32,
    }
}

/// The compounds the old program saves for a run's `comps` (`masses[p]` and `peaks[p]` as in
/// [`SavedComponents`], `ris[i]` the RI of component `i`, -1 when unknown). The old map keeps one peak
/// per ion mass, the first, and writes them in ascending ion order. The name stays empty: a
/// component has none.
// shortcut: the model peak carries scan range and values only; its other fields and the noise
// factor are 0 (not recovered from the old program), and `count` is 1 for every ion peak.
pub fn compounds_of_components(
    comps: &[Component],
    masses: &[u32],
    peaks: &[Peak],
    scan_times_s: &[f64],
    ris: &[f64],
) -> Vec<LibraryCompound> {
    let stamp = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_secs() as i64);
    comps
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let mut ions: Vec<LibraryPeak> = Vec::new();
            for f in &c.peaks {
                let (Some(p), Some(&mass)) = (peaks.get(f.peak), masses.get(f.peak)) else {
                    continue;
                };
                if ions.iter().all(|x| x.ion != mass as f32) {
                    ions.push(LibraryPeak {
                        ion: mass as f32,
                        sharpness: f64::from(p.sharpness),
                        integral: f.integral,
                        quality: f64::from(p.quality),
                        signal_to_noise: f64::from(p.signal_to_noise),
                        count: 1,
                        first_scan: p.first_scan,
                        last_scan: p.last_scan,
                        model_parameter: f.model_share,
                        max_scan: p.apex_scan,
                        interpolated_max_scan: p
                            .precise_apex
                            .map_or(f64::from(p.apex_scan), f64::from),
                    });
                }
            }
            ions.sort_by(|a, b| a.ion.total_cmp(&b.ion));
            let mut model = zero_peak();
            model.ion = 1.0;
            model.count = 1;
            (model.first_scan, model.max_scan, model.last_scan) =
                (c.first_scan, c.apex_scan, c.last_scan);
            model.interpolated_max_scan = f64::from(c.apex_scan);
            let mut compound = LibraryCompound::new("");
            compound.spectrum = LibrarySpectrum {
                scan: c.apex_scan,
                rt_ms: scan_times_s
                    .get(c.apex_scan as usize)
                    .map_or(0.0, |t| t * 1000.0),
                ri: ris.get(i).copied().unwrap_or(-1.0),
                bin_width: 1.0,
                masses: c.spectrum.iter().map(|p| p.0 as f32).collect(),
                values: (c.spectrum.iter())
                    .map(|p| p.1.min(i32::MAX as u32) as i32)
                    .collect(),
            };
            compound.peaks = ions;
            compound.model = Some(LibraryModel {
                peak: model,
                noise_factor: 0.0,
                values: c.model.clone(),
            });
            compound.time_stamp = stamp;
            compound.id = i as i64;
            compound
        })
        .collect()
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

/// Writes a `.cmp` file (gzip) the old program can open. It is written beside `path` first and then
/// renamed, so a failed write never leaves a broken file in its place.
pub fn write_cmp(path: impl AsRef<Path>, compounds: &[LibraryCompound]) -> std::io::Result<()> {
    let path = path.as_ref();
    let mut part = path.as_os_str().to_owned();
    part.push(".part");
    let write = || {
        let file = std::fs::File::create(&part)?;
        let mut gz = flate2::write::GzEncoder::new(file, flate2::Compression::default());
        gz.write_all(&cmp_payload(compounds))?;
        gz.finish()?.sync_all()?;
        std::fs::rename(&part, path)
    };
    write().inspect_err(|_| {
        std::fs::remove_file(&part).ok();
    })
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
        let (k, masses, peaks) = saved_components(&run, &times).unwrap().unwrap();
        assert_eq!(
            (peaks.len(), peaks[1].signal_to_noise, peaks[1].apex_scan),
            (2, 42.0, 12)
        );
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
    fn components_survive_saving_with_their_peak_values() {
        use crate::deconv::FittedPeak;
        use crate::detail::component_detail;

        let pk = |apex, quality, sharpness, snr| Peak {
            row: 5,
            first_scan: apex - 2,
            apex_scan: apex,
            last_scan: apex + 2,
            precise_apex: Some(apex as f32 + 0.5),
            signal_to_noise: snr,
            quality,
            sharpness,
        };
        // Mass 147 sits twice (§4.5): the first, sharper one counts.
        let (masses, peaks) = (
            vec![147, 73, 147, 207],
            vec![
                pk(10, 0.5, 9.0, 40.0),
                pk(10, 1.5, 7.0, 20.0),
                pk(11, 2.5, 3.0, 10.0),
                pk(10, 0.0, 5.0, 30.0),
            ],
        );
        let fit = |peak, integral| FittedPeak {
            peak,
            model_share: 0.25,
            integral,
        };
        let comps = vec![
            Component {
                first_scan: 8,
                apex_scan: 10,
                last_scan: 12,
                model: vec![1.0, 5.0, 9.0, 5.0, 1.0],
                spectrum: vec![(73, 500), (147, 999), (207, 40)],
                peaks: vec![fit(0, 900.0), fit(1, 400.0), fit(2, 100.0), fit(3, 30.0)],
                exact_spectrum: Vec::new(),
            },
            Component {
                first_scan: 3,
                apex_scan: 4,
                last_scan: 5,
                model: vec![1.0, 2.0, 1.0],
                spectrum: vec![(50, 10)],
                peaks: vec![],
                exact_spectrum: Vec::new(),
            },
        ];
        let times: Vec<f64> = (0..20).map(|i| 100.0 + f64::from(i) * 0.5).collect();
        let compounds = compounds_of_components(&comps, &masses, &peaks, &times, &[1500.0]);
        assert_eq!(compounds[0].spectrum.ri, 1500.0);
        assert_eq!(compounds[1].spectrum.ri, -1.0, "no RI given");
        assert_eq!(compounds[0].spectrum.rt_ms, 105_000.0);
        let ions: Vec<f32> = compounds[0].peaks.iter().map(|p| p.ion).collect();
        assert_eq!(ions, [73.0, 147.0, 207.0], "one peak per mass, ascending");

        let dir = std::env::temp_dir().join(format!("md-core-comps-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("run.cmp");
        write_cmp(&file, &compounds).unwrap();
        let ((back, m2, p2), standalone_times) = standalone_components(&file).unwrap();
        // Next to its run the file is checked against the run's times; alone it rebuilds the apex times.
        let (back2, ..) = saved_components(&dir.join("run.cdf"), &times)
            .unwrap()
            .unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(back, back2);
        assert_eq!(standalone_times[10], 105.0);
        // Components come back ordered by apex, so the one saved first is second now.
        assert_eq!((back[1].apex_scan, back[0].apex_scan), (10, 4));
        let (mut before, after) = (
            component_detail(&comps[0], &peaks, &masses),
            component_detail(&back[1], &p2, &m2),
        );
        // The old file keeps its ion peaks ascending by mass, not sharpest first.
        before.peaks.sort_by_key(|r| r.mass);
        assert_eq!(
            before, after,
            "panel values from the file equal the originals"
        );
        assert_eq!(after.peaks.len(), 3);
        assert_eq!(back[1].spectrum, comps[0].spectrum);
        assert_eq!(back[1].model, comps[0].model);
        assert_eq!(p2[0].precise_apex, Some(10.5));
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
