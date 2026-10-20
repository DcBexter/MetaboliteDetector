//! Compound libraries: the old binary `.lbr` format (read), AMDIS MSL text (import and export) and the
//! library text format `.mdlib` (read and write). Layouts: `docs/reverse-engineering/file-formats.md`
//! §1 (`.lbr`) and §6 (MSL); `.mdlib` in `docs/architecture.md` (Libraries).

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::Read;
use std::path::Path;

#[derive(Debug, thiserror::Error)]
pub enum LibraryError {
    #[error("cannot read library file: {0}")]
    Io(#[from] std::io::Error),
    #[error("not supported: {0}")]
    Unsupported(String),
    #[error("corrupt library at byte {offset}: {reason}")]
    Corrupt { offset: usize, reason: String },
    #[error("MSL line {line}: {reason}")]
    Msl { line: usize, reason: String },
    #[error("library line {line}: {reason}")]
    Mdlib { line: usize, reason: String },
}

/// Reference spectrum of a library compound (`Spectrum<int,float>`).
#[derive(Debug, Clone, PartialEq)]
pub struct LibrarySpectrum {
    pub scan: u32,
    /// Retention time in ms (unit inferred from the old GUI).
    pub rt_ms: f64,
    pub ri: f64,
    pub bin_width: f32,
    /// m/z values; `values[i]` is the intensity at `masses[i]`.
    pub masses: Vec<f32>,
    pub values: Vec<i32>,
}

/// One library entry (`LibraryCompound<int,float>`), fields in file order.
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryCompound {
    pub spectrum: LibrarySpectrum,
    pub name: String,
    pub fwhm: Option<Vec<f64>>,
    pub integrals: Option<Vec<f64>>,
    /// Sorted, unique (`std::set<float>`).
    pub quant_ions: Vec<f32>,
    pub auto_quant_ions: Vec<f32>,
    pub features: BTreeMap<String, String>,
    pub drt: f64,
    pub time_stamp: i64,
    pub id: i64,
    pub comb_spec: Vec<(f32, i32)>,
    pub ion_freq: Vec<(f32, u32)>,
    pub spec_count: u32,
    pub ris: Vec<f64>,
    /// Version 2 peak records (`PeakEssential<float>`, file-formats.md §1.3), kept verbatim.
    pub blocks: Vec<[u8; BLOCK_BYTES]>,
    pub sources: Vec<String>,
    /// Peaks of the compound (version 4 and later); sorted by ion in `.lbr` files.
    pub peaks: Vec<LibraryPeak>,
    /// Model peak (version 4 and later).
    pub model: Option<LibraryModel>,
    /// Fragment `(begin, end, formula)`, formula trimmed (version 3 and later).
    pub formulas: Vec<(u32, u32, String)>,
    /// Fragment `(begin, end, mid_size)` (version 5).
    pub mid_sizes: Vec<(u32, u32, u32)>,
}

/// One peak (`PeakEssential`), fields in file order (file-formats.md §1.6).
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryPeak {
    pub ion: f32,
    pub sharpness: f64,
    pub integral: f64,
    pub quality: f64,
    pub signal_to_noise: f64,
    pub count: u32,
    pub first_scan: u32,
    pub last_scan: u32,
    pub model_parameter: f64,
    pub max_scan: u32,
    pub interpolated_max_scan: f64,
}

/// Model peak of a compound (`Peak<double,short>`); its ion is a 16-bit integer in the file.
#[derive(Debug, Clone, PartialEq)]
pub struct LibraryModel {
    pub peak: LibraryPeak,
    pub noise_factor: f64,
    pub values: Vec<f64>,
}

impl LibraryCompound {
    /// A compound without spectrum, with the defaults the old MSL importer left (sample library).
    pub fn new(name: &str) -> Self {
        Self {
            spectrum: LibrarySpectrum {
                scan: 0,
                // -1 is "unknown", what the old importers set when a record has no RT or RI.
                rt_ms: -1.0,
                ri: -1.0,
                bin_width: 1.0,
                masses: Vec::new(),
                values: Vec::new(),
            },
            name: name.to_string(),
            fwhm: None,
            integrals: None,
            quant_ions: Vec::new(),
            auto_quant_ions: Vec::new(),
            features: BTreeMap::new(),
            drt: 0.0,
            time_stamp: 0,
            id: 0,
            comb_spec: Vec::new(),
            ion_freq: Vec::new(),
            spec_count: 1,
            ris: Vec::new(),
            blocks: Vec::new(),
            sources: Vec::new(),
            peaks: Vec::new(),
            model: None,
            formulas: Vec::new(),
            mid_sizes: Vec::new(),
        }
    }
}

const BLOCK_BYTES: usize = 40;
/// `PeakEssential<float>` in versions 4 and 5.
const PEAK_BYTES: usize = 68;
/// Size of a record whose lists and strings are all empty: a lower bound per compound.
const MIN_COMPOUND_BYTES: usize = 36 + 56 + 24;
// shortcut: whole library in RAM, capped so a gzip bomb cannot exhaust memory; raise it when real
// libraries come close.
const MAX_LBR_BYTES: u64 = 1 << 30;

/// Reads a `.lbr` file (gzip, or an already decompressed payload).
pub fn read_lbr(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    parse_lbr(&std::fs::read(path)?)
}

/// Parses `.lbr` bytes of versions 2 to 5 (the old program writes 5); every length and check value is
/// validated.
pub fn parse_lbr(bytes: &[u8]) -> Result<Vec<LibraryCompound>, LibraryError> {
    let unzipped;
    let data = if bytes.starts_with(&[0x1f, 0x8b]) {
        unzipped = gunzip(bytes)?;
        &unzipped[..]
    } else {
        bytes
    };
    let mut r = Reader::new(data);
    let version = r.u32()?;
    if !(2..=5).contains(&version) {
        return Err(LibraryError::Unsupported(format!(
            "library version {version}"
        )));
    }
    let compounds = r.list(MIN_COMPOUND_BYTES, |r| read_compound(r, version))?;
    if r.remaining() != 0 {
        return r.corrupt(format!("{} bytes after the last compound", r.remaining()));
    }
    Ok(compounds)
}

pub(crate) fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, LibraryError> {
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes)
        .take(MAX_LBR_BYTES + 1)
        .read_to_end(&mut out)
        .map_err(|e| LibraryError::Corrupt {
            offset: 0,
            reason: format!("gzip: {e}"),
        })?;
    if out.len() as u64 > MAX_LBR_BYTES {
        return Err(LibraryError::Unsupported(format!(
            "libraries larger than {MAX_LBR_BYTES} bytes uncompressed"
        )));
    }
    Ok(out)
}

/// `LibraryCompound::fromStream`; the version gates follow the old reader (file-formats.md §1.2).
fn read_compound(r: &mut Reader, version: u32) -> Result<LibraryCompound, LibraryError> {
    let mut c = read_compound_record(r, version)?;
    c.comb_spec = r.list(8, |r| Ok((r.f32()?, r.i32()?)))?;
    c.ion_freq = r.list(8, |r| Ok((r.f32()?, r.u32()?)))?;
    c.spec_count = r.u32()?;
    c.ris = r.list(8, Reader::f64)?;
    if version <= 3 {
        c.blocks = r.list(BLOCK_BYTES, Reader::array)?;
    }
    c.sources = r.list(4, Reader::string)?;
    // The old reader trims each formula (`trimFormulas`, C `isspace`, which includes VT).
    if version >= 3 {
        c.formulas = r.list(12, |r| {
            Ok((
                r.u32()?,
                r.u32()?,
                r.string()?
                    .trim_matches(|c: char| c.is_ascii_whitespace() || c == '\x0b')
                    .to_string(),
            ))
        })?;
    }
    if version >= 5 {
        c.mid_sizes = r.list(12, |r| Ok((r.u32()?, r.u32()?, r.u32()?)))?;
    }
    Ok(c)
}

/// `Compound<int,float>::fromStream` (file-formats.md §3.2), the part a library compound shares with
/// a `.cmp` compound; the library-only fields keep the defaults of [`LibraryCompound::new`].
pub(crate) fn read_compound_record(
    r: &mut Reader,
    version: u32,
) -> Result<LibraryCompound, LibraryError> {
    let scan = r.u32()?;
    let rt_ms = r.f64()?;
    let ri = r.f64()?;
    let bin_width = r.f32()?;
    let masses = r.list(4, Reader::f32)?;
    let values = r.list(4, Reader::i32)?;
    if r.u32()? != scan.wrapping_add(masses.len() as u32) {
        return r.corrupt("spectrum check value mismatch");
    }
    if masses.len() != values.len() {
        return r.corrupt("spectrum has different numbers of masses and values");
    }
    let spectrum = LibrarySpectrum {
        scan,
        rt_ms,
        ri,
        bin_width,
        masses,
        values,
    };

    let n_peaks = r.u32()?;
    if n_peaks != 0 && version < 4 {
        return Err(LibraryError::Unsupported(
            "version 2 and 3 compounds with peak records (layout unknown)".into(),
        ));
    }
    let n = r.count(n_peaks as usize, PEAK_BYTES)?;
    let peaks = (0..n)
        .map(|_| read_peak(r, false))
        .collect::<Result<_, _>>()?;
    let model = if version >= 4 {
        Some(read_model(r)?)
    } else {
        None
    };
    let name_bytes = r.bytes()?;
    let fwhm = r.f64_list_or_null()?;
    let integrals = r.f64_list_or_null()?;
    let quant_ions = r.list(4, Reader::f32)?;
    let auto_quant_ions = r.list(4, Reader::f32)?;
    let feature_list = r.list(8, |r| Ok((r.string()?, r.string()?)))?;
    let drt = r.f64()?;
    let time_stamp = r.i64()?;
    let id = r.i64()?;
    let expected = (quant_ions.len() + feature_list.len() + name_bytes.len()) as u32;
    if r.u32()? != expected.wrapping_add(n_peaks) {
        return r.corrupt("compound check value mismatch");
    }

    Ok(LibraryCompound {
        spectrum,
        name: text(name_bytes),
        fwhm,
        integrals,
        quant_ions,
        auto_quant_ions,
        features: feature_list.into_iter().collect(),
        drt,
        time_stamp,
        id,
        peaks,
        model,
        ..LibraryCompound::new("")
    })
}

/// `PeakEssential<float>` (68 bytes) or, with `short_ion`, `PeakEssential<short>` (66 bytes).
pub(crate) fn read_peak(r: &mut Reader, short_ion: bool) -> Result<LibraryPeak, LibraryError> {
    Ok(LibraryPeak {
        ion: if short_ion {
            f32::from(i16::from_le_bytes(r.array()?))
        } else {
            r.f32()?
        },
        sharpness: r.f64()?,
        integral: r.f64()?,
        quality: r.f64()?,
        signal_to_noise: r.f64()?,
        count: r.u32()?,
        first_scan: r.u32()?,
        last_scan: r.u32()?,
        model_parameter: r.f64()?,
        max_scan: r.u32()?,
        interpolated_max_scan: r.f64()?,
    })
}

/// `Peak<double,short>`: peak, noise factor, values, check = first + last + max scan + value count.
pub(crate) fn read_model(r: &mut Reader) -> Result<LibraryModel, LibraryError> {
    let peak = read_peak(r, true)?;
    let noise_factor = r.f64()?;
    let values = r.list(8, Reader::f64)?;
    let expected = (peak.first_scan)
        .wrapping_add(peak.last_scan)
        .wrapping_add(peak.max_scan)
        .wrapping_add(values.len() as u32);
    if r.u32()? != expected {
        return r.corrupt("model peak check value mismatch");
    }
    Ok(LibraryModel {
        peak,
        noise_factor,
        values,
    })
}

/// Little-endian reader that fails instead of panicking on short or inconsistent input.
pub(crate) struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    pub(crate) fn new(data: &'a [u8]) -> Self {
        Self { data, pos: 0 }
    }

    pub(crate) fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    pub(crate) fn corrupt<T>(&self, reason: impl Into<String>) -> Result<T, LibraryError> {
        Err(LibraryError::Corrupt {
            offset: self.pos,
            reason: reason.into(),
        })
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8], LibraryError> {
        if n > self.remaining() {
            return self.corrupt(format!("truncated, {n} more bytes expected"));
        }
        self.pos += n;
        Ok(&self.data[self.pos - n..self.pos])
    }

    pub(crate) fn array<const N: usize>(&mut self) -> Result<[u8; N], LibraryError> {
        let mut a = [0; N];
        a.copy_from_slice(self.take(N)?);
        Ok(a)
    }

    pub(crate) fn u32(&mut self) -> Result<u32, LibraryError> {
        self.array().map(u32::from_le_bytes)
    }
    pub(crate) fn i32(&mut self) -> Result<i32, LibraryError> {
        self.array().map(i32::from_le_bytes)
    }
    pub(crate) fn i64(&mut self) -> Result<i64, LibraryError> {
        self.array().map(i64::from_le_bytes)
    }
    pub(crate) fn f32(&mut self) -> Result<f32, LibraryError> {
        self.array().map(f32::from_le_bytes)
    }
    pub(crate) fn f64(&mut self) -> Result<f64, LibraryError> {
        self.array().map(f64::from_le_bytes)
    }

    /// Reads a count of elements of at least `min_bytes` each; a count that cannot fit in the
    /// remaining bytes is rejected before anything is allocated.
    pub(crate) fn count(&mut self, n: usize, min_bytes: usize) -> Result<usize, LibraryError> {
        if n.checked_mul(min_bytes)
            .is_none_or(|b| b > self.remaining())
        {
            return self.corrupt(format!("count {n} exceeds the remaining data"));
        }
        Ok(n)
    }

    pub(crate) fn list<T>(
        &mut self,
        min_bytes: usize,
        mut item: impl FnMut(&mut Self) -> Result<T, LibraryError>,
    ) -> Result<Vec<T>, LibraryError> {
        let n = self.u32()? as usize;
        let n = self.count(n, min_bytes)?;
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(item(self)?);
        }
        Ok(out)
    }

    /// i32 count (negative = null pointer) followed by f64 values.
    pub(crate) fn f64_list_or_null(&mut self) -> Result<Option<Vec<f64>>, LibraryError> {
        let Ok(n) = usize::try_from(self.i32()?) else {
            return Ok(None);
        };
        let n = self.count(n, 8)?;
        (0..n)
            .map(|_| self.f64())
            .collect::<Result<_, _>>()
            .map(Some)
    }

    pub(crate) fn bytes(&mut self) -> Result<&'a [u8], LibraryError> {
        let n = self.u32()? as usize;
        self.take(n)
    }

    pub(crate) fn string(&mut self) -> Result<String, LibraryError> {
        self.bytes().map(text)
    }
}

/// The old files carry no encoding: UTF-8 when valid, otherwise Latin-1 (every byte is one char).
fn text(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect())
}

/// Opens a library by extension: `.lbr`, `.mdlib`, NIST text (`.sdf`, see [`is_nist`]), else AMDIS
/// `.msl` or NIST MSP text.
pub fn open(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    let path = path.as_ref();
    let ext = path.extension().unwrap_or_default().to_ascii_lowercase();
    match ext.to_str() {
        Some("lbr") => read_lbr(path),
        Some("mdlib") => read_mdlib(path),
        Some("sdf") => import_nist(path),
        Some("txt") if is_targetsearch(path) => import_targetsearch(path),
        _ => import_msl(path),
    }
}

/// True for the NIST text the old `SQLLibrary::importNIST` reads; the old search over it used the
/// fragment pre-filter (`PreparedLibrary::with_fragment_filter`).
pub fn is_nist(path: impl AsRef<Path>) -> bool {
    path.as_ref()
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("sdf"))
}

/// Imports NIST text (see [`parse_nist`]); each compound records the file as its source.
pub fn import_nist(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    let path = path.as_ref();
    let mut compounds = parse_nist(&text(&std::fs::read(path)?))?;
    for c in &mut compounds {
        c.sources
            .push(format!("Initial import from{}", path.display()));
    }
    Ok(compounds)
}

/// Parses the NIST text of the old `gcms::NISTParser` (file-formats.md §6.3): records end with a
/// `$$$$` line; a field starts with a `> <FIELD NAME>` line and its value is the following lines,
/// trimmed and joined with `\n`. Lines before the first field (an SDF mol block) are ignored. Like the
/// old `SQLLibrary::importNIST`, only peaks with m/z > 69 and intensity > 0 are kept and records without
/// such a peak are skipped. A file without any kept record is an error (e.g. a chemistry SDF file).
pub fn parse_nist(text: &str) -> Result<Vec<LibraryCompound>, LibraryError> {
    let mut out = Vec::new();
    let mut fields: Vec<(&str, String)> = Vec::new();
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    for line in text.lines().chain(["$$$$"]) {
        let line = line.trim_matches([' ', '\n', '\r', '\t']);
        if line == "$$$$" {
            out.extend(nist_record(&fields));
            fields.clear();
        } else if let Some(name) = line
            .strip_prefix('>')
            .and_then(|r| r.trim_start().strip_prefix('<'))
            .and_then(|r| r.strip_suffix('>'))
        {
            fields.push((name, String::new()));
        } else if let Some((_, v)) = fields.last_mut().filter(|_| !line.is_empty()) {
            if !v.is_empty() {
                v.push('\n');
            }
            v.push_str(line);
        }
    }
    if out.is_empty() {
        return Err(LibraryError::Unsupported(
            "no NIST record with a peak above m/z 69".into(),
        ));
    }
    Ok(out)
}

fn nist_record(fields: &[(&str, String)]) -> Option<LibraryCompound> {
    let get = |k: &str| fields.iter().find(|f| f.0 == k).map(|f| f.1.as_str());
    let mut c = LibraryCompound::new(get("NAME").unwrap_or_default());
    let s = &mut c.spectrum;
    s.ri = get("ESTIMATED KOVATS RI")
        .and_then(|v| v.parse().ok())
        .unwrap_or(-1.0);
    s.rt_ms = -1.0;
    // The old importer reads m/z with strtol (truncates) and the intensity with strtod.
    for line in get("MASS SPECTRAL PEAKS").unwrap_or_default().lines() {
        let mut it = line.split_whitespace().map(str::parse::<f64>);
        if let (Some(Ok(mz)), Some(Ok(v))) = (it.next(), it.next())
            && mz.trunc() > 69.0
            && mz < 1e6
            && v > 0.0
        {
            s.masses.push(mz.trunc() as f32);
            s.values.push(v.min(i32::MAX as f64) as i32);
        }
    }
    if s.masses.is_empty() {
        return None;
    }
    for (k, v) in fields {
        if !matches!(
            *k,
            "NAME" | "ESTIMATED KOVATS RI" | "MASS SPECTRAL PEAKS" | "NUM PEAKS"
        ) {
            c.features.insert(k.to_string(), v.clone());
        }
    }
    copy_to_combined(&mut c);
    Some(c)
}

/// Imports an AMDIS `.msl` file; each compound records the file as its source like the old app.
pub fn import_msl(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    let path = path.as_ref();
    let mut compounds = parse_msl(&text(&std::fs::read(path)?))?;
    for c in &mut compounds {
        c.sources
            .push(format!("Initial import from{}", path.display()));
    }
    Ok(compounds)
}

/// Imports TargetSearch's tab-separated `library.txt` (see [`parse_targetsearch`]).
pub fn import_targetsearch(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    let path = path.as_ref();
    let mut compounds = parse_targetsearch(&text(&std::fs::read(path)?))?;
    for c in &mut compounds {
        c.sources
            .push(format!("Initial import from{}", path.display()));
    }
    Ok(compounds)
}

const TARGETSEARCH_HEADER: &str = "Name\tRI\tWin_1\tSEL_MASS\tSPECTRUM";

/// True for a `.txt` file with TargetSearch's header line.
pub fn is_targetsearch(path: impl AsRef<Path>) -> bool {
    use std::io::Read;
    let path = path.as_ref();
    let mut head = Vec::new();
    // Only the first bytes: this also runs on dropped run files.
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("txt"))
        && std::fs::File::open(path)
            .and_then(|f| f.take(64).read_to_end(&mut head))
            .is_ok_and(|_| {
                text(&head)
                    .trim_start_matches('\u{feff}')
                    .starts_with(TARGETSEARCH_HEADER)
            })
}

/// Parses TargetSearch's `library.txt`: the header line [`TARGETSEARCH_HEADER`], then one line per compound
/// with those tab-separated columns: `RI` is on the FAME scale of the file and kept as is, `Win_1` is kept
/// as feature `WIN_1`, `SEL_MASS` (`;`-separated) gives the quantification ions, `SPECTRUM` is `mz:intensity`
/// pairs. It is rewritten as MSL and parsed by [`parse_msl`], so errors in values report that MSL's line.
pub fn parse_targetsearch(text: &str) -> Result<Vec<LibraryCompound>, LibraryError> {
    let text = text.trim_start_matches('\u{feff}');
    if text.lines().next().map(str::trim_end) != Some(TARGETSEARCH_HEADER) {
        return Err(LibraryError::Msl {
            line: 1,
            reason: format!("expected the header `{TARGETSEARCH_HEADER}`"),
        });
    }
    let mut msl = String::new();
    for (i, line) in text.lines().enumerate().skip(1) {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').map(str::trim).collect();
        let [name, ri, win, sel, spectrum] = cols[..] else {
            return Err(LibraryError::Msl {
                line: i + 1,
                reason: "expected 5 tab-separated columns".into(),
            });
        };
        let pairs: Vec<_> = spectrum.split_whitespace().collect();
        msl.push_str(&format!(
            "NAME: {name}\nRI: {ri}\nWIN_1: {win}\nQUANTIFICATION: {}\nNUM PEAKS: {}\n",
            sel.replace(';', " "),
            pairs.len()
        ));
        for p in pairs {
            msl.push_str(&format!("({})", p.replace(':', " ")));
        }
        msl.push_str("\n\n");
    }
    parse_msl(&msl)
}

/// Parses AMDIS MSL text, and NIST MSP text (peak lines `m/z intensity` separated by `;` or spaces,
/// quoted annotations ignored). `RT:` is in minutes; `RI:` and the MSP spellings `Retention_index:` and
/// `RetentionIndex:` set the retention index; keys without a field of their own become features,
/// repeated keys (MSP `Synon:`) joined with `\n`. Fields MSL has no key for keep the defaults the old
/// importer left (sample library).
pub fn parse_msl(text: &str) -> Result<Vec<LibraryCompound>, LibraryError> {
    let mut out = Vec::new();
    // Record in progress, its NAME line and its NUM PEAKS value.
    let mut cur: Option<(LibraryCompound, usize, Option<usize>)> = None;
    for (i, line) in text.lines().enumerate() {
        let (line, n) = (line.trim(), i + 1);
        let err = |reason: &str| LibraryError::Msl {
            line: n,
            reason: reason.into(),
        };
        if line.is_empty() {
            continue;
        }
        if line.starts_with('(') {
            let Some((c, _, Some(_))) = cur.as_mut() else {
                return Err(err("peak list before NAME: and NUM PEAKS:"));
            };
            parse_peaks(line, &mut c.spectrum).map_err(err)?;
            continue;
        }
        if line.starts_with(|c: char| c.is_ascii_digit())
            && let Some((c, _, Some(_))) = cur.as_mut()
        {
            parse_msp_peaks(line, &mut c.spectrum).map_err(err)?;
            continue;
        }
        let Some((key, value)) = line.split_once(':') else {
            return Err(err("expected `KEY: value` or `(m/z intensity)`"));
        };
        let (key, value) = (key.trim().to_ascii_uppercase(), value.trim());
        if key == "NAME" {
            out.extend(finish_msl(cur.take())?);
            cur = Some((LibraryCompound::new(value), n, None));
            continue;
        }
        let Some((c, _, peaks)) = cur.as_mut() else {
            return Err(err("expected NAME: before other keys"));
        };
        let num = |v: &str| v.parse::<f64>().map_err(|_| err("invalid number"));
        match key.as_str() {
            "RI" => c.spectrum.ri = num(value)?,
            // MSP spellings the old program keeps as features; an explicit `RI:` or an unreadable
            // number wins over them.
            "RETENTION_INDEX" | "RETENTIONINDEX" | "RETENTION INDEX"
                if c.spectrum.ri < 0.0 && value.parse::<f64>().is_ok_and(f64::is_finite) =>
            {
                c.spectrum.ri = num(value)?;
            }
            "RT" => c.spectrum.rt_ms = num(value)? * 60_000.0,
            "RESOLUTION" => c.spectrum.bin_width = num(value)? as f32,
            "NUM PEAKS" => {
                *peaks = Some(value.parse().map_err(|_| err("invalid peak count"))?);
            }
            "QUANTIFICATION" => {
                let mut ions = value
                    .split_whitespace()
                    .map(|v| v.parse::<f32>().map_err(|_| err("invalid ion")))
                    .collect::<Result<Vec<_>, _>>()?;
                ions.sort_by(f32::total_cmp);
                ions.dedup();
                c.quant_ions = ions;
            }
            _ => {
                c.features
                    .entry(key)
                    .and_modify(|v| {
                        v.push('\n');
                        v.push_str(value);
                    })
                    .or_insert_with(|| value.to_string());
            }
        }
    }
    out.extend(finish_msl(cur)?);
    Ok(out)
}

/// Appends the `(m/z intensity)` pairs of one peak line to the spectrum.
fn parse_peaks(line: &str, s: &mut LibrarySpectrum) -> Result<(), &'static str> {
    for chunk in line.split('(').skip(1) {
        let (pair, rest) = chunk.split_once(')').ok_or("missing `)`")?;
        if !rest.trim().is_empty() {
            return Err("text between peak pairs");
        }
        let mut it = pair.split_whitespace();
        let (Some(mz), Some(v), None) = (it.next(), it.next(), it.next()) else {
            return Err("expected `(m/z intensity)`");
        };
        push_peak(s, mz, v)?;
    }
    Ok(())
}

/// Appends the MSP pairs of one peak line: `m/z intensity` pairs separated by `;` or spaces, each pair
/// optionally followed by a quoted annotation (`71 100 "C5H11+"`).
fn parse_msp_peaks(line: &str, s: &mut LibrarySpectrum) -> Result<(), &'static str> {
    if !line.matches('"').count().is_multiple_of(2) {
        return Err("unclosed `\"` annotation");
    }
    // Every second `"`-separated piece is an annotation.
    let numbers: Vec<_> = line
        .split('"')
        .step_by(2)
        .flat_map(|p| p.split(|c: char| c == ';' || c.is_whitespace()))
        .filter(|t| !t.is_empty())
        .collect();
    if !numbers.len().is_multiple_of(2) {
        return Err("expected `m/z intensity` pairs");
    }
    for pair in numbers.chunks(2) {
        push_peak(s, pair[0], pair[1])?;
    }
    Ok(())
}

fn push_peak(s: &mut LibrarySpectrum, mz: &str, v: &str) -> Result<(), &'static str> {
    let mz: f32 = mz.parse().map_err(|_| "invalid m/z")?;
    let v: f64 = v.parse().map_err(|_| "invalid intensity")?;
    if !mz.is_finite() || !(0.0..=i32::MAX as f64).contains(&v) {
        return Err("m/z or intensity out of range");
    }
    s.masses.push(mz);
    s.values.push(v.round() as i32);
    Ok(())
}

fn finish_msl(
    cur: Option<(LibraryCompound, usize, Option<usize>)>,
) -> Result<Option<LibraryCompound>, LibraryError> {
    let Some((mut c, line, peaks)) = cur else {
        return Ok(None);
    };
    if peaks != Some(c.spectrum.masses.len()) {
        return Err(LibraryError::Msl {
            line,
            reason: format!(
                "NUM PEAKS {peaks:?} but {} pairs found",
                c.spectrum.masses.len()
            ),
        });
    }
    copy_to_combined(&mut c);
    Ok(Some(c))
}

/// The old MSL importer copies the spectrum into the combined spectrum (sample library).
fn copy_to_combined(c: &mut LibraryCompound) {
    let s = &c.spectrum;
    c.comb_spec = s
        .masses
        .iter()
        .copied()
        .zip(s.values.iter().copied())
        .collect();
    c.ion_freq = s.masses.iter().map(|&m| (m, 1)).collect();
}

/// Writes compounds as AMDIS MSL text that [`parse_msl`] reads back with identical spectra.
pub fn write_msl(compounds: &[LibraryCompound]) -> String {
    // A line break inside a value would start a new record on import.
    let one_line = |s: &str| s.replace(['\r', '\n'], " ");
    let mut s = String::new();
    for c in compounds {
        let sp = &c.spectrum;
        let _ = writeln!(s, "NAME: {}", one_line(&c.name));
        let _ = writeln!(s, "RI: {}", sp.ri);
        let _ = writeln!(s, "RT: {}", sp.rt_ms / 60_000.0);
        if !c.quant_ions.is_empty() {
            let ions: Vec<_> = c.quant_ions.iter().map(f32::to_string).collect();
            let _ = writeln!(s, "QUANTIFICATION: {}", ions.join(" "));
        }
        let _ = writeln!(s, "RESOLUTION: {}", sp.bin_width);
        for (k, v) in &c.features {
            let _ = writeln!(s, "{}: {}", one_line(k), one_line(v));
        }
        let _ = writeln!(s, "NUM PEAKS: {}", sp.masses.len());
        let pairs: Vec<_> = sp
            .masses
            .iter()
            .zip(&sp.values)
            .map(|(m, v)| format!("({m} {v})"))
            .collect();
        for line in pairs.chunks(5) {
            let _ = writeln!(s, "{}", line.join(" "));
        }
        s.push('\n');
    }
    s
}

/// The three export types of the old "Export Library" dialog (`LibrarySearch::exportMSL`,
/// docs/reverse-engineering/file-formats.md section 6.2).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MslType {
    /// "MSL Nominal": m/z as whole numbers.
    Nominal,
    /// "MSL Profile": m/z with six decimals.
    Profile,
    /// "MSL FWHM": m/z, FWHM and intensity per peak, no `RESOLUTION:` line.
    Fwhm,
}

/// `Compound::getNormalizedIntensities`: the base peak becomes 1000, the rest is rounded, and a peak
/// never drops to 0.
fn normalized_intensities(values: &[i32]) -> Vec<i32> {
    let base = f64::from(values.iter().copied().max().unwrap_or(0));
    // shortcut: the old program's result for a base peak of 0 or less is not known, keep the values.
    if base <= 0.0 {
        return values.to_vec();
    }
    values
        .iter()
        .map(|&v| match (f64::from(v) / base * 1000.0).round() as i32 {
            0 => 1,
            n => n,
        })
        .collect()
}

/// `Spectrum::getBinnedSpectrum` of the old program with bin width 1: the peaks of one whole mass
/// become one peak with the largest intensity. A bin starts at `trunc((mass + 0.5) / 1)`, so x.5 rounds
/// up, and the masses must be sorted.
fn bin_to_whole_masses(masses: &[f32], values: &[i32]) -> (Vec<f32>, Vec<i32>) {
    let bin = |m: f32| ((f64::from(m) + 0.5) as i32) as f32;
    let (mut out_m, mut out_v) = (Vec::new(), Vec::new());
    let (Some(&first), Some(&first_v)) = (masses.first(), values.first()) else {
        return (out_m, out_v);
    };
    let (mut bin_mass, mut best) = (bin(first), first_v);
    for (i, (&m, &v)) in masses.iter().zip(values).enumerate() {
        if f64::from(m) >= f64::from(bin_mass) + 0.5 {
            out_m.push(bin_mass);
            out_v.push(best);
            (bin_mass, best) = (bin(m), v);
        } else if v > best {
            best = v;
        }
        if i == masses.len() - 1 {
            out_m.push(bin_mass);
            out_v.push(best);
        }
    }
    (out_m, out_v)
}

/// `Compound::getCentroidCompound` of the old program for a compound without FWHM list. Points whose
/// neighbours are at most 0.01 apart form one profile peak (rising to the apex, then falling to the
/// valley). Its mass is the intensity-weighted mean of the points at or above half the group maximum,
/// its intensity their sum and its FWHM their count; a lone point stays as it is with FWHM 1.
/// Returns masses, intensities and FWHMs.
fn centroid(masses: &[f32], values: &[i32]) -> (Vec<f32>, Vec<i32>, Vec<f64>) {
    let n = masses.len().min(values.len());
    let (mut out_m, mut out_v, mut out_w) = (Vec::new(), Vec::new(), Vec::new());
    let mut i = 0;
    while i < n {
        let (j, max) = group_end(&masses[..n], &values[..n], i);
        let end = if j >= n { j - 1 } else { j };
        let peak = if i == end {
            Some((masses[i], values[i], 1.0))
        } else {
            weighted_peak(&masses[i..=end], &values[i..=end], max)
        };
        if let Some((m, v, w)) = peak {
            out_m.push(m);
            out_v.push(v);
            out_w.push(w);
        }
        i = j + 1;
    }
    (out_m, out_v, out_w)
}

/// Where the profile peak starting at `start` ends (`n` when it runs to the end) and its largest value.
fn group_end(masses: &[f32], values: &[i32], start: usize) -> (usize, f64) {
    let n = masses.len();
    // The old code reads past the end as "the same as the last point".
    let at = |k: usize| values[k.min(n - 1)];
    let mass_at = |k: usize| masses[k.min(n - 1)];
    let (mut j, mut max, mut rising) = (start, f64::MIN, true);
    while j < n {
        let (v0, v1, v2) = (at(j), at(j + 1), at(j + 2));
        max = max.max(f64::from(v0));
        if f64::from(mass_at(j + 1) - mass_at(j)) > 0.01 {
            break;
        }
        // Climb to the apex (two falling steps in a row end the climb), then fall to the valley
        // (a clear rise ends the peak).
        if rising && v1 < v0 && v2 < v1 {
            rising = false;
        } else if !rising && v1 > v0 && v2 > v1 {
            break;
        } else {
            j += 1;
        }
    }
    (j, max)
}

/// Mass, intensity and FWHM of several points; the old code drops a peak without weight.
fn weighted_peak(ms: &[f32], vs: &[i32], max: f64) -> Option<(f32, i32, f64)> {
    let half = max / 2.0;
    let (mut weighted, mut sum, mut width) = (0.0_f64, 0.0_f64, 0.0_f64);
    for (&m, &v) in ms.iter().zip(vs).filter(|&(_, &v)| f64::from(v) >= half) {
        weighted += f64::from(m * v as f32);
        sum += f64::from(v);
        width += 1.0;
    }
    (sum != 0.0 && weighted != 0.0).then(|| ((weighted / sum) as f32, sum as i32, width))
}

/// Writes compounds as the old program does for `kind`, byte for byte. Unlike [`write_msl`] this is for
/// other programs (AMDIS): it rounds RT, RI and quantification ions to four decimals, so reading it
/// back is not exact, and a `;` inside a feature value starts a new `KEY:` line.
pub fn write_msl_as(compounds: &[LibraryCompound], kind: MslType) -> String {
    // A line break inside a value would start a new record on import.
    let one_line = |s: &str| s.replace(['\r', '\n'], " ");
    let mut s = String::new();
    for c in compounds {
        let sp = &c.spectrum;
        let _ = write!(
            s,
            "NAME:{}\nRI:{:.4}\nRT:{:.4}\n",
            one_line(&c.name),
            sp.ri,
            // The old program writes `RT:-1.0000` for an unknown RT.
            if sp.rt_ms < 0.0 {
                -1.0
            } else {
                sp.rt_ms / 60_000.0
            }
        );
        for (k, v) in &c.features {
            let (k, v) = (one_line(k), one_line(v));
            let _ = writeln!(s, "{k}: {}", v.replace(';', &format!("\n{k}: ")));
        }
        s.push_str("QUANTIFICATION:");
        for q in &c.quant_ions {
            let _ = write!(s, "{q:.4}\t");
        }
        if kind != MslType::Fwhm {
            // Binning to whole masses sets the bin width of the exported compound to 1.
            let width = if kind == MslType::Nominal {
                1.0
            } else {
                sp.bin_width
            };
            let _ = write!(s, "\nRESOLUTION:{width:.4}");
        }
        // The old program keeps the peaks in a map, so they always come out sorted by mass.
        let mut order: Vec<usize> = (0..sp.masses.len().min(sp.values.len())).collect();
        order.sort_by(|&x, &y| sp.masses[x].total_cmp(&sp.masses[y]));
        let sorted_m: Vec<f32> = order.iter().map(|&i| sp.masses[i]).collect();
        let sorted_v: Vec<i32> = order.iter().map(|&i| sp.values[i]).collect();
        let (masses, values, fwhm) = match kind {
            MslType::Profile => (sorted_m, sorted_v, Vec::new()),
            // The old program centroids first unless the compound already has a FWHM list.
            _ => match &c.fwhm {
                Some(f) => {
                    let f = order
                        .iter()
                        .map(|&i| f.get(i).copied().unwrap_or(1.0))
                        .collect();
                    (sorted_m, sorted_v, f)
                }
                None => centroid(&sorted_m, &sorted_v),
            },
        };
        let (masses, values) = if kind == MslType::Nominal {
            bin_to_whole_masses(&masses, &values)
        } else {
            (masses, values)
        };
        let values = normalized_intensities(&values);
        let _ = write!(s, "\nNUM PEAKS:{}", masses.len());
        for (i, (m, v)) in masses.iter().zip(&values).enumerate() {
            if i % 5 == 0 {
                s.push('\n');
            }
            let _ = match kind {
                MslType::Nominal => write!(s, "({m:>4.0} {v:>4}) "),
                MslType::Profile => write!(s, "({m:>10.6} {v:>4}) "),
                MslType::Fwhm => {
                    let w = fwhm.get(i).copied().unwrap_or(1.0);
                    write!(s, "({m:>10.6} {w:>5.2}{v:>4}) ")
                }
            };
        }
        s.push_str("\n\n");
    }
    // The old program ends the file with one more blank line, also for an empty library.
    s.push_str("\n\n");
    s
}

/// First line of a `.mdlib` file; the number is the format version.
const MDLIB_HEADER: &str = "MetaboliteDetector library 2";

/// Reads a `.mdlib` library file.
pub fn read_mdlib(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    let path = path.as_ref();
    if std::fs::metadata(path)?.len() > MAX_LBR_BYTES {
        return Err(LibraryError::Unsupported(format!(
            "libraries larger than {MAX_LBR_BYTES} bytes"
        )));
    }
    parse_mdlib(&std::fs::read_to_string(path)?)
}

/// Saves `compounds` as `.mdlib`. Writes `<path>.tmp` first and renames it over `path`, so a failed
/// save never leaves half a library behind.
pub fn save_mdlib(
    path: impl AsRef<Path>,
    compounds: &[LibraryCompound],
) -> Result<(), LibraryError> {
    let path = path.as_ref();
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".tmp");
    // Synced before the rename, so a crash cannot leave an empty file under the real name.
    let result = std::fs::File::create(&tmp)
        .and_then(|mut f| {
            std::io::Write::write_all(&mut f, write_mdlib(compounds).as_bytes())?;
            f.sync_all()
        })
        .and_then(|()| std::fs::rename(&tmp, path));
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    Ok(result?)
}

/// Writes compounds as `.mdlib` text; [`parse_mdlib`] reads every field back bit for bit, except
/// that a NaN loses its sign and payload.
pub fn write_mdlib(compounds: &[LibraryCompound]) -> String {
    fn join<T: std::fmt::Display>(items: impl IntoIterator<Item = T>) -> String {
        let items: Vec<_> = items.into_iter().map(|x| x.to_string()).collect();
        items.join(" ")
    }
    fn pairs<A: std::fmt::Display, B: std::fmt::Display>(p: &[(A, B)]) -> String {
        join(p.iter().map(|(a, b)| format!("{a}:{b}")))
    }
    let opt = |v: &Option<Vec<f64>>| v.as_ref().map_or("null".into(), join);
    let mut s = format!("{MDLIB_HEADER}\n");
    for c in compounds {
        let sp = &c.spectrum;
        let peaks: Vec<_> = sp.masses.iter().zip(&sp.values).collect();
        let mut lines = vec![
            ("compound", String::new()),
            ("name", escape(&c.name)),
            ("scan", sp.scan.to_string()),
            ("rt_ms", sp.rt_ms.to_string()),
            ("ri", sp.ri.to_string()),
            ("bin_width", sp.bin_width.to_string()),
            ("peaks", pairs(&peaks)),
            ("fwhm", opt(&c.fwhm)),
            ("integrals", opt(&c.integrals)),
            ("quant_ions", join(&c.quant_ions)),
            ("auto_quant_ions", join(&c.auto_quant_ions)),
        ];
        for (k, v) in &c.features {
            lines.push(("feature", format!("{}\t{}", escape(k), escape(v))));
        }
        lines.extend([
            ("drt", c.drt.to_string()),
            ("time_stamp", c.time_stamp.to_string()),
            ("id", c.id.to_string()),
            ("comb_spec", pairs(&c.comb_spec)),
            ("ion_freq", pairs(&c.ion_freq)),
            ("spec_count", c.spec_count.to_string()),
            ("ris", join(&c.ris)),
        ]);
        for b in &c.blocks {
            lines.push(("block", b.iter().map(|x| format!("{x:02x}")).collect()));
        }
        for src in &c.sources {
            lines.push(("source", escape(src)));
        }
        for p in &c.peaks {
            lines.push(("peak", peak_fields(p)));
        }
        if let Some(m) = &c.model {
            let values = join(&m.values);
            let sep = if values.is_empty() { "" } else { " " };
            let fields = format!("{} {}{sep}{values}", peak_fields(&m.peak), m.noise_factor);
            lines.push(("model", fields));
        }
        for (begin, end, formula) in &c.formulas {
            lines.push(("formula", format!("{begin} {end} {}", escape(formula))));
        }
        for (begin, end, size) in &c.mid_sizes {
            lines.push(("mid_size", format!("{begin} {end} {size}")));
        }
        lines.push(("end", String::new()));
        for (key, value) in lines {
            s.push_str(key);
            if !value.is_empty() {
                s.push(' ');
                s.push_str(&value);
            }
            s.push('\n');
        }
    }
    s
}

/// The 11 peak fields in file order, space separated.
fn peak_fields(p: &LibraryPeak) -> String {
    format!(
        "{} {} {} {} {} {} {} {} {} {} {}",
        p.ion,
        p.sharpness,
        p.integral,
        p.quality,
        p.signal_to_noise,
        p.count,
        p.first_scan,
        p.last_scan,
        p.model_parameter,
        p.max_scan,
        p.interpolated_max_scan
    )
}

/// Parses `.mdlib` text. Strict: unknown keys, repeated keys and a newer version are errors.
pub fn parse_mdlib(text: &str) -> Result<Vec<LibraryCompound>, LibraryError> {
    // Editors on Windows may put a byte order mark in front.
    let text = text.strip_prefix('\u{feff}').unwrap_or(text);
    let mut lines = (text.lines().enumerate())
        .map(|(i, l)| (i + 1, l))
        .filter(|(_, l)| !l.trim().is_empty() && !l.starts_with('#'));
    let err = |line: usize, reason: String| LibraryError::Mdlib { line, reason };
    const PREFIX: &str = "MetaboliteDetector library ";
    match lines.next() {
        // Version 2 only added keys (peak, model, formula, mid_size), so version 1 reads unchanged.
        Some((_, l)) if [MDLIB_HEADER, "MetaboliteDetector library 1"].contains(&l.trim_end()) => {}
        Some((_, l)) if l.starts_with(PREFIX) => {
            return Err(LibraryError::Unsupported(format!(
                "library text format version {}",
                &l[PREFIX.len()..]
            )));
        }
        other => {
            let line = other.map_or(1, |(n, _)| n);
            return Err(err(line, format!("expected `{MDLIB_HEADER}`")));
        }
    }
    let mut out = Vec::new();
    // Compound in progress, the line it started on, and the single-use keys it has seen.
    let mut cur: Option<(LibraryCompound, usize, Vec<&str>)> = None;
    for (n, line) in lines {
        let (key, value) = line.split_once(' ').unwrap_or((line, ""));
        match (key, cur.as_mut()) {
            ("compound", None) => cur = Some((LibraryCompound::new(""), n, Vec::new())),
            ("end", Some((_, start, seen))) => {
                if !seen.contains(&"name") {
                    return Err(err(*start, "compound without `name`".into()));
                }
                out.extend(cur.take().map(|(c, ..)| c));
            }
            (_, None) => return Err(err(n, "expected `compound`".into())),
            ("compound", Some(_)) => return Err(err(n, "`compound` before `end`".into())),
            (key, Some((c, _, seen))) => {
                if !matches!(
                    key,
                    "feature" | "block" | "source" | "peak" | "formula" | "mid_size"
                ) {
                    if seen.contains(&key) {
                        return Err(err(n, format!("`{key}` given twice")));
                    }
                    seen.push(key);
                }
                set_mdlib_field(c, key, value).map_err(|r| err(n, r))?;
            }
        }
    }
    if let Some((_, start, _)) = cur {
        return Err(err(start, "compound without `end`".into()));
    }
    Ok(out)
}

fn num<T: std::str::FromStr>(v: &str) -> Result<T, String> {
    v.parse().map_err(|_| format!("invalid number `{v}`"))
}

fn set_mdlib_field(c: &mut LibraryCompound, key: &str, v: &str) -> Result<(), String> {
    fn list<T: std::str::FromStr>(v: &str) -> Result<Vec<T>, String> {
        v.split_whitespace().map(num).collect()
    }
    fn pairs<A: std::str::FromStr, B: std::str::FromStr>(v: &str) -> Result<Vec<(A, B)>, String> {
        (v.split_whitespace())
            .map(|p| {
                let (a, b) =
                    (p.split_once(':')).ok_or(format!("expected `m/z:value`, got `{p}`"))?;
                Ok((num(a)?, num(b)?))
            })
            .collect()
    }
    let opt = |v: &str| match v {
        "null" => Ok(None),
        v => list(v).map(Some),
    };
    let s = &mut c.spectrum;
    match key {
        "name" => c.name = unescape(v)?,
        "scan" => s.scan = num(v)?,
        "rt_ms" => s.rt_ms = num(v)?,
        "ri" => s.ri = num(v)?,
        "bin_width" => s.bin_width = num(v)?,
        "peaks" => (s.masses, s.values) = pairs::<f32, i32>(v)?.into_iter().unzip(),
        "fwhm" => c.fwhm = opt(v)?,
        "integrals" => c.integrals = opt(v)?,
        "quant_ions" => {
            c.quant_ions = list(v)?;
            c.quant_ions.sort_by(f32::total_cmp);
            c.quant_ions.dedup();
        }
        "auto_quant_ions" => c.auto_quant_ions = list(v)?,
        "feature" => {
            let (k, value) = v.split_once('\t').ok_or("expected `KEY<tab>value`")?;
            // An empty key is valid: MSL import makes one from a line like `: x`.
            let k = unescape(k)?;
            if c.features.contains_key(&k) {
                return Err(format!("feature `{k}` given twice"));
            }
            c.features.insert(k, unescape(value)?);
        }
        "drt" => c.drt = num(v)?,
        "time_stamp" => c.time_stamp = num(v)?,
        "id" => c.id = num(v)?,
        "comb_spec" => c.comb_spec = pairs(v)?,
        "ion_freq" => c.ion_freq = pairs(v)?,
        "spec_count" => c.spec_count = num(v)?,
        "ris" => c.ris = list(v)?,
        "block" => {
            if v.len() != 2 * BLOCK_BYTES || !v.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Err(format!("block must be {} hex digits", 2 * BLOCK_BYTES));
            }
            let mut b = [0; BLOCK_BYTES];
            for (i, x) in b.iter_mut().enumerate() {
                *x = u8::from_str_radix(&v[2 * i..2 * i + 2], 16).map_err(|e| e.to_string())?;
            }
            c.blocks.push(b);
        }
        "source" => c.sources.push(unescape(v)?),
        "peak" => {
            let f: Vec<&str> = v.split_whitespace().collect();
            if f.len() != 11 {
                return Err("peak must have 11 fields".into());
            }
            c.peaks.push(parse_peak(&f)?);
        }
        "model" => {
            let f: Vec<&str> = v.split_whitespace().collect();
            if f.len() < 12 {
                return Err("model must have at least 12 fields".into());
            }
            c.model = Some(LibraryModel {
                peak: parse_peak(&f[..11])?,
                noise_factor: num(f[11])?,
                values: f[12..].iter().map(|x| num(x)).collect::<Result<_, _>>()?,
            });
        }
        "formula" => {
            let mut f = v.splitn(3, ' ');
            let (Some(b), Some(e)) = (f.next(), f.next()) else {
                return Err("expected `begin end formula`".into());
            };
            // The formula may be empty.
            c.formulas
                .push((num(b)?, num(e)?, unescape(f.next().unwrap_or(""))?));
        }
        "mid_size" => match list::<u32>(v)?[..] {
            [b, e, size] => c.mid_sizes.push((b, e, size)),
            _ => return Err("expected `begin end size`".into()),
        },
        _ => return Err(format!("unknown key `{key}`")),
    }
    Ok(())
}

fn parse_peak(f: &[&str]) -> Result<LibraryPeak, String> {
    Ok(LibraryPeak {
        ion: num(f[0])?,
        sharpness: num(f[1])?,
        integral: num(f[2])?,
        quality: num(f[3])?,
        signal_to_noise: num(f[4])?,
        count: num(f[5])?,
        first_scan: num(f[6])?,
        last_scan: num(f[7])?,
        model_parameter: num(f[8])?,
        max_scan: num(f[9])?,
        interpolated_max_scan: num(f[10])?,
    })
}

/// Keeps text values on one line and frees the tab to separate a feature's key from its value.
fn escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c => out.push(c),
        }
    }
    out
}

fn unescape(s: &str) -> Result<String, String> {
    let mut out = String::with_capacity(s.len());
    let mut it = s.chars();
    while let Some(ch) = it.next() {
        out.push(match ch {
            '\\' => match it.next() {
                Some('\\') => '\\',
                Some('n') => '\n',
                Some('r') => '\r',
                Some('t') => '\t',
                Some(other) => return Err(format!("unknown escape `\\{other}`")),
                None => return Err("`\\` at the end of a value".into()),
            },
            c => c,
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    #[test]
    fn targetsearch_library_txt_opens() {
        let dir = std::env::temp_dir().join(format!("md-core-ts-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let p = dir.join("library.txt");
        std::fs::write(
            &p,
            "Name\tRI\tWin_1\tSEL_MASS\tSPECTRUM\r\nGlycine (2TMS)\t228554\t4000\t86;102\t86:26 102:999\r\n",
        )
        .unwrap();
        let libs = super::open(&p).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(libs.len(), 1);
        let c = &libs[0];
        assert_eq!(c.name, "Glycine (2TMS)");
        assert_eq!(c.spectrum.ri, 228554.0);
        assert_eq!(c.quant_ions, [86.0, 102.0]);
        assert_eq!(c.spectrum.masses, [86.0, 102.0]);
        assert_eq!(c.features["WIN_1"], "4000");
        assert!(super::parse_targetsearch("Name\tRI\nx\t1\n").is_err());
        let swapped = "Name\tRI\tSEL_MASS\tWin_1\tSPECTRUM\nx\t1\t86\t4000\t86:1\n";
        assert!(super::parse_targetsearch(swapped).is_err());
        let bom = "\u{feff}Name\tRI\tWin_1\tSEL_MASS\tSPECTRUM\nx\t1\t4000\t86\t86:1\n";
        assert_eq!(super::parse_targetsearch(bom).unwrap().len(), 1);
    }

    #[test]
    fn open_picks_the_reader_by_extension() {
        let dir = std::env::temp_dir().join(format!("md-core-open-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let msl = dir.join("lib.MSL");
        std::fs::write(&msl, "NAME: A\nNUM PEAKS: 1\n(100 1)\n").unwrap();
        assert_eq!(super::open(&msl).unwrap()[0].name, "A");
        // The same bytes are no .lbr file.
        let lbr = dir.join("lib.LBR");
        std::fs::copy(&msl, &lbr).unwrap();
        assert!(super::open(&lbr).is_err());
        let mdlib = dir.join("lib.MDLIB");
        super::save_mdlib(&mdlib, &super::open(&msl).unwrap()).unwrap();
        assert_eq!(super::open(&mdlib).unwrap(), super::open(&msl).unwrap());
        std::fs::remove_dir_all(&dir).unwrap();
    }

    use super::*;

    /// Version 2 payload with one compound, written like the old `toStream`.
    fn sample_payload() -> Vec<u8> {
        payload(2)
    }

    /// One compound as the old writer of `version` lays it out: from version 4 on with one peak and
    /// a model peak instead of the version 2 block, version 3 adds formulas and version 5 mid sizes.
    fn payload(version: u32) -> Vec<u8> {
        let mut b = Vec::new();
        let u32 = |b: &mut Vec<u8>, v: u32| b.extend(v.to_le_bytes());
        let string = |b: &mut Vec<u8>, s: &str| {
            b.extend((s.len() as u32).to_le_bytes());
            b.extend(s.as_bytes());
        };
        u32(&mut b, version);
        u32(&mut b, 1);
        // Spectrum: scan 0, rt, ri, bin width, 2 masses, 2 values, check.
        u32(&mut b, 0);
        b.extend(436080.0f64.to_le_bytes());
        b.extend(1369.05f64.to_le_bytes());
        b.extend(1.0f32.to_le_bytes());
        u32(&mut b, 2);
        b.extend(73.0f32.to_le_bytes());
        b.extend(204.1433f32.to_le_bytes());
        u32(&mut b, 2);
        b.extend(1000i32.to_le_bytes());
        b.extend(250i32.to_le_bytes());
        u32(&mut b, 2);
        // Compound: peaks and model (v4+), name, null fwhm/integrals, 1 quant ion, no auto ions,
        // 1 feature.
        let peak = |b: &mut Vec<u8>, ion: &[u8], first: u32| {
            b.extend(ion);
            for v in [14.1, 7558399.3, 0.006, 433.0] {
                b.extend(f64::to_le_bytes(v));
            }
            for v in [1, first, first + 194] {
                u32(b, v);
            }
            b.extend(0.89f64.to_le_bytes());
            u32(b, first + 54);
            b.extend(2222.7f64.to_le_bytes());
        };
        let n_peaks = u32::from(version >= 4);
        u32(&mut b, n_peaks);
        if version >= 4 {
            peak(&mut b, &87.0f32.to_le_bytes(), 2169);
            peak(&mut b, &73i16.to_le_bytes(), 10);
            b.extend(0.5f64.to_le_bytes());
            u32(&mut b, 2);
            b.extend(1.5f64.to_le_bytes());
            b.extend(2.5f64.to_le_bytes());
            u32(&mut b, 10 + 204 + 64 + 2);
        }
        string(&mut b, "Serine 3TMS");
        b.extend((-1i32).to_le_bytes());
        b.extend((-1i32).to_le_bytes());
        u32(&mut b, 1);
        b.extend(204.1433f32.to_le_bytes());
        u32(&mut b, 0);
        u32(&mut b, 1);
        string(&mut b, "LIB_SCORE");
        string(&mut b, "0.98");
        b.extend(0.0f64.to_le_bytes());
        b.extend(1227546843i64.to_le_bytes());
        b.extend(2197790076i64.to_le_bytes());
        u32(&mut b, 1 + 1 + 11 + n_peaks);
        // comb_spec, ion_freq, spec_count, ris, blocks (v2, v3), sources, formulas, mid sizes.
        u32(&mut b, 1);
        b.extend(73.0f32.to_le_bytes());
        b.extend(1000i32.to_le_bytes());
        u32(&mut b, 1);
        b.extend(73.0f32.to_le_bytes());
        u32(&mut b, 1);
        u32(&mut b, 1);
        u32(&mut b, 1);
        b.extend(1369.0f64.to_le_bytes());
        if version <= 3 {
            u32(&mut b, 1);
            b.extend([7; BLOCK_BYTES]);
        }
        u32(&mut b, 1);
        string(&mut b, "Initial import from/home/kh/Desktop/test.msl");
        if version >= 3 {
            u32(&mut b, 1);
            u32(&mut b, 73);
            u32(&mut b, 147);
            string(&mut b, "\x0b C3H9Si ");
        }
        if version >= 5 {
            u32(&mut b, 1);
            for v in [73, 147, 2] {
                u32(&mut b, v);
            }
        }
        b
    }

    #[test]
    fn reads_versions_3_to_5() {
        let v2 = parse_lbr(&sample_payload()).unwrap().remove(0);
        let v3 = parse_lbr(&payload(3)).unwrap().remove(0);
        assert_eq!(v3.formulas, [(73, 147, "C3H9Si".to_string())]);
        assert_eq!(
            LibraryCompound {
                formulas: Vec::new(),
                ..v3
            },
            v2
        );

        let v5 = parse_lbr(&payload(5)).unwrap().remove(0);
        let [p] = &v5.peaks[..] else { panic!() };
        assert_eq!(
            (p.ion, p.integral, p.signal_to_noise),
            (87.0, 7558399.3, 433.0)
        );
        assert_eq!((p.first_scan, p.last_scan, p.max_scan), (2169, 2363, 2223));
        assert_eq!((p.model_parameter, p.interpolated_max_scan), (0.89, 2222.7));
        let m = v5.model.as_ref().unwrap();
        assert_eq!(
            (m.peak.ion, m.peak.first_scan, m.noise_factor),
            (73.0, 10, 0.5)
        );
        assert_eq!(m.values, [1.5, 2.5]);
        assert_eq!(v5.mid_sizes, [(73, 147, 2)]);
        assert!(v5.blocks.is_empty());
        assert_eq!(
            (&v5.name, &v5.sources, &v5.formulas),
            (&v2.name, &v2.sources, &v3.formulas)
        );

        let v4 = parse_lbr(&payload(4)).unwrap().remove(0);
        assert_eq!(
            v4,
            LibraryCompound {
                mid_sizes: Vec::new(),
                ..v5
            }
        );
    }

    #[test]
    fn version_5_checks_counts_and_check_values() {
        let payload = payload(5);
        for len in 0..payload.len() {
            assert!(parse_lbr(&payload[..len]).is_err(), "prefix {len} accepted");
        }
        let set = |at: usize, v: u32| {
            let mut p = payload.clone();
            p[at..at + 4].copy_from_slice(&v.to_le_bytes());
            parse_lbr(&p).unwrap_err().to_string()
        };
        // Offsets: spectrum ends at 60, the peak count, 68-byte peak, 66-byte model peak.
        assert!(set(60, u32::MAX).contains("exceeds"), "peak count");
        let model = 64 + 68;
        assert!(
            set(model + 74, u32::MAX).contains("exceeds"),
            "model values"
        );
        assert!(set(model + 94, 0).contains("model peak check"));
        // The compound check value counts the peaks.
        // After the model: name, fwhm, integrals, 2 ion sets, 1 feature, drt, time stamp, id.
        let check = model + 98 + 15 + 4 + 4 + 8 + 4 + 4 + 21 + 24;
        assert!(set(check, 13).contains("compound check"));
    }

    #[test]
    fn reads_version_2_plain_and_gzip() {
        let payload = sample_payload();
        let mut gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
        std::io::Write::write_all(&mut gz, &payload).unwrap();
        let gz = gz.finish().unwrap();
        let lib = parse_lbr(&payload).unwrap();
        assert_eq!(parse_lbr(&gz).unwrap(), lib);
        let [c] = &lib[..] else { panic!() };
        assert_eq!(c.name, "Serine 3TMS");
        assert_eq!(c.spectrum.masses, [73.0, 204.1433]);
        assert_eq!(c.spectrum.values, [1000, 250]);
        assert_eq!((c.spectrum.ri, c.spectrum.rt_ms), (1369.05, 436080.0));
        assert_eq!(c.features["LIB_SCORE"], "0.98");
        assert_eq!(
            (c.fwhm.as_ref(), c.quant_ions.as_slice()),
            (None, &[204.1433][..])
        );
        assert_eq!(
            (c.time_stamp, c.id, c.spec_count),
            (1227546843, 2197790076, 1)
        );
        assert_eq!(
            (c.ion_freq.as_slice(), c.ris.as_slice()),
            (&[(73.0, 1)][..], &[1369.0][..])
        );
        assert_eq!(c.blocks, [[7; BLOCK_BYTES]]);
        assert_eq!(c.sources, ["Initial import from/home/kh/Desktop/test.msl"]);
    }

    #[test]
    fn rejects_truncated_and_manipulated_input_without_panicking() {
        let payload = sample_payload();
        for len in 0..payload.len() {
            assert!(parse_lbr(&payload[..len]).is_err(), "prefix {len} accepted");
        }
        let mut extra = payload.clone();
        extra.push(0);
        assert!(parse_lbr(&extra).is_err());
        // Every byte flipped: may parse or fail, must never panic or allocate wildly.
        for i in 0..payload.len() {
            let mut p = payload.clone();
            p[i] ^= 0xff;
            let _ = parse_lbr(&p);
        }
        let set = |at: usize, v: u32| {
            let mut p = payload.clone();
            p[at..at + 4].copy_from_slice(&v.to_le_bytes());
            parse_lbr(&p).unwrap_err().to_string()
        };
        assert!(set(4, u32::MAX).contains("exceeds"), "compound count");
        assert!(set(32, u32::MAX).contains("exceeds"), "mass count");
        assert!(set(56, 3).contains("check value"), "spectrum check");
        assert!(set(60, 1).contains("peak records"));
        assert!(set(0, 6).contains("version 6"));
        assert!(
            parse_lbr(&[0x1f, 0x8b, 0, 0])
                .unwrap_err()
                .to_string()
                .contains("gzip")
        );
    }

    const MSL: &str = "NAME: Glycine 2TMS\r\nRI: 1129.256\nRT: 5.663\nQUANTIFICATION: 174 102.1144\n\
        COMMENT: from AMDIS: test\nNUM PEAKS: 6\n( 45  99) ( 73 733) (74 10)\n(102 1000) (147 462.4)\n\
        (174.05 1)\n\nNAME: Empty\nNUM PEAKS: 0\n";

    #[test]
    fn msl_import_and_round_trip() {
        let lib = parse_msl(MSL).unwrap();
        let [g, e] = &lib[..] else { panic!() };
        assert_eq!(g.name, "Glycine 2TMS");
        assert_eq!(g.spectrum.ri, 1129.256);
        assert!((g.spectrum.rt_ms - 339_780.0).abs() < 1e-6);
        assert_eq!(g.quant_ions, [102.1144, 174.0]);
        assert_eq!(g.features["COMMENT"], "from AMDIS: test");
        assert_eq!(g.spectrum.masses, [45.0, 73.0, 74.0, 102.0, 147.0, 174.05]);
        assert_eq!(g.spectrum.values, [99, 733, 10, 1000, 462, 1]);
        assert_eq!(g.comb_spec[1], (73.0, 733));
        assert!(e.spectrum.masses.is_empty());

        let mut again = parse_msl(&write_msl(&lib)).unwrap();
        for (a, b) in again.iter_mut().zip(&lib) {
            assert!((a.spectrum.rt_ms - b.spectrum.rt_ms).abs() < 1e-6);
            a.spectrum.rt_ms = b.spectrum.rt_ms;
        }
        assert_eq!(again, lib);
    }

    /// MSP files (MoNA, MS-DIAL) write the retention index as `Retention_index:`, which the old program
    /// only keeps as a feature. Synthetic records, no library data.
    #[test]
    fn msp_retention_index_spellings_set_the_ri() {
        let text = "Name: a
Retention_index: 1817.1
Num Peaks: 1
73 999

Name: b
RETENTIONINDEX: 1234.5
RI: 1500
Num Peaks: 1
73 999

Name: c
RI: 1500
Retention_index: 1700
Num Peaks: 1
73 999

Name: d
Retention_index: n/a
Num Peaks: 1
73 999

Name: e
Retention_index: inf
Num Peaks: 1
73 999
";
        let lib = parse_msl(text).unwrap();
        let ri: Vec<_> = lib.iter().map(|c| c.spectrum.ri).collect();
        // `RI:` wins over a spelling, whatever the order; a spelling with no number stays a feature.
        assert_eq!(ri, [1817.1, 1500.0, 1500.0, -1.0, -1.0]);
        assert_eq!(lib[4].features["RETENTION_INDEX"], "inf");
        assert!(!lib[0].features.contains_key("RETENTION_INDEX"));
        assert_eq!(lib[2].features["RETENTION_INDEX"], "1700");
        assert_eq!(lib[3].features["RETENTION_INDEX"], "n/a");
    }

    /// The old `importMSL` leaves RT and RI at -1 (unknown) when a record has none; 0 stays 0.
    /// Checked with the old binary (docs/reverse-engineering/tools/golden/make-msl-export.sh).
    #[test]
    fn missing_rt_and_ri_are_unknown() {
        let text = "NAME: a
NUM PEAKS: 1
(73 1000)

NAME: b
RI: 0
RT: 0
NUM PEAKS: 1
(73 1000)
";
        let lib = parse_msl(text).unwrap();
        assert_eq!((lib[0].spectrum.ri, lib[0].spectrum.rt_ms), (-1.0, -1.0));
        assert_eq!((lib[1].spectrum.ri, lib[1].spectrum.rt_ms), (0.0, 0.0));
        let new = LibraryCompound::new("n");
        assert_eq!((new.spectrum.ri, new.spectrum.rt_ms), (-1.0, -1.0));
        // The export writes the unknown values the way the old program does.
        let out = write_msl_as(&lib[..1], MslType::Profile);
        assert!(
            out.contains(
                "RI:-1.0000
RT:-1.0000
"
            ),
            "{out}"
        );
    }

    #[test]
    fn msl_rejects_broken_records() {
        let line = |text: &str| match parse_msl(text) {
            Err(LibraryError::Msl { line, .. }) => line,
            r => panic!("{text:?} gave {r:?}"),
        };
        assert_eq!(line("NAME: a\nNUM PEAKS: 2\n(1 2)\n"), 1);
        assert_eq!(line("(1 2)\n"), 1);
        assert_eq!(line("RI: 5\n"), 1);
        assert_eq!(line("NAME: a\nRI: x\n"), 2);
        assert_eq!(line("NAME: a\nNUM PEAKS: 1\n(1 2\n"), 3);
        assert_eq!(line("NAME: a\nNUM PEAKS: 1\n(1 -2)\n"), 3);
        assert_eq!(line("NAME: a\nNUM PEAKS: 1\n(1 2 3)\n"), 3);
        assert_eq!(line("NAME: a\nno colon\n"), 2);
        assert_eq!(line("NAME: a\nNUM PEAKS: 99999999999999999999\n"), 2);
    }

    #[test]
    fn msl_export_keeps_one_record_per_compound() {
        let mut c = parse_msl("NAME: x\nNUM PEAKS: 0\n").unwrap().remove(0);
        c.name = "two\nNAME: lines".into();
        assert_eq!(parse_msl(&write_msl(&[c])).unwrap().len(), 1);
    }

    #[test]
    fn old_msl_export_types_differ_in_peak_layout() {
        let mut c = parse_msl("NAME: x\nRI: 1200.5\nRT: 2\nQUANTIFICATION: 73\nRESOLUTION: 1\nNUM PEAKS: 2\n(73 999) (147.5 1)\n")
            .unwrap()
            .remove(0);
        c.features.insert("K".into(), "a;b".into());
        let head = "NAME:x\nRI:1200.5000\nRT:2.0000\nK: a\nK: b\nQUANTIFICATION:73.0000\t";
        let end = "\n\n\n\n";
        let out = |c: &LibraryCompound, k| write_msl_as(std::slice::from_ref(c), k);
        assert_eq!(
            out(&c, MslType::Nominal),
            format!("{head}\nRESOLUTION:1.0000\nNUM PEAKS:2\n(  73 1000) ( 148    1) {end}")
        );
        assert_eq!(
            out(&c, MslType::Profile),
            format!(
                "{head}\nRESOLUTION:1.0000\nNUM PEAKS:2\n( 73.000000 1000) (147.500000    1) {end}"
            )
        );
        c.fwhm = Some(vec![0.5]);
        assert_eq!(
            out(&c, MslType::Fwhm),
            format!("{head}\nNUM PEAKS:2\n( 73.000000  0.501000) (147.500000  1.00   1) {end}")
        );
        assert_eq!(write_msl_as(&[], MslType::Nominal), "\n\n");
    }

    /// The old program keeps peaks sorted by mass, whatever order the file has.
    #[test]
    fn old_msl_export_sorts_the_peaks_by_mass() {
        let c = parse_msl(
            "NAME: u
NUM PEAKS: 3
(100 10) (50 1000) (73 20)
",
        )
        .unwrap();
        for (kind, peaks) in [
            (MslType::Nominal, "(  50 1000) (  73   20) ( 100   10)"),
            (
                MslType::Profile,
                "( 50.000000 1000) ( 73.000000   20) (100.000000   10)",
            ),
        ] {
            assert!(write_msl_as(&c, kind).contains(peaks), "{kind:?}");
        }
    }

    /// The expected peaks are what the old program wrote for the same spectra
    /// (docs/reverse-engineering/tools/golden/make-msl-export.sh with a `.msl` input).
    #[test]
    fn old_msl_export_bins_and_centroids_like_the_old_program() {
        let text = "NAME: f\nNUM PEAKS: 9\n(41.1 80) (73.2 500) (73.4 300) (74 20) (99.9 10)\n\
                    (100.004 200) (100.008 100) (100.5 1000) (147.5 50)\n\n\
                    NAME: c\nNUM PEAKS: 4\n(56.995 111) (57 1000) (57.004 333) (58 11)\n";
        let lib = parse_msl(text).unwrap();
        let peaks = |c: &LibraryCompound, k| {
            let out = write_msl_as(std::slice::from_ref(c), k);
            out.split("NUM PEAKS:")
                .nth(1)
                .unwrap()
                .split_once('\n')
                .unwrap()
                .1
                .replace('\n', "|")
        };
        assert_eq!(
            peaks(&lib[0], MslType::Nominal),
            "(  41   80) (  73  500) (  74   20) ( 100  300) ( 101 1000) |( 148   50) ||||"
        );
        assert_eq!(
            peaks(&lib[0], MslType::Fwhm),
            "( 41.099998  1.00  80) ( 73.199997  1.00 500) ( 73.400002  1.00 300) ( 74.000000  1.00  20) ( 99.900002  1.00  10) \
             |(100.005333  2.00 300) (100.500000  1.001000) (147.500000  1.00  50) ||||"
        );
        assert!(peaks(&lib[0], MslType::Profile).contains("(100.003998  200) (100.008003  100)"));
        assert_eq!(
            peaks(&lib[1], MslType::Nominal),
            "(  57 1000) (  58   11) ||||"
        );
        assert_eq!(
            peaks(&lib[1], MslType::Fwhm),
            "( 57.000000  1.001000) ( 58.000000  1.00  11) ||||"
        );
    }

    #[test]
    fn mdlib_keeps_every_field_bit_for_bit() {
        let mut lib = parse_lbr(&sample_payload()).unwrap();
        let mut odd = lib[0].clone();
        odd.name = "a\\b\nc\td\r e ".into();
        odd.features.insert("K\tE\\Y".into(), "v\nw".into());
        odd.features
            .insert(String::new(), "from `: x` in MSL".into());
        (odd.fwhm, odd.integrals) = (Some(vec![]), Some(vec![1.5, -0.0, f64::INFINITY]));
        (odd.spectrum.ri, odd.drt, odd.spectrum.bin_width) = (-1.0, f64::NAN, 0.1);
        odd.sources.clear();
        lib.push(odd);
        let mut v5 = parse_lbr(&payload(5)).unwrap().remove(0);
        v5.formulas.push((1, 2, String::new()));
        v5.formulas.push((3, 4, "a b	c".into()));
        lib.push(v5);
        let mut no_values = parse_lbr(&payload(5)).unwrap().remove(0);
        no_values.model.as_mut().unwrap().values.clear();
        lib.push(no_values);
        lib.push(LibraryCompound::new("bare"));
        let text = write_mdlib(&lib);
        let again = parse_mdlib(&text).unwrap();
        // Text equality also covers NaN and -0.0, which `==` cannot tell apart.
        assert_eq!(write_mdlib(&again), text);
        assert_eq!(again[0], lib[0]);
        assert_eq!(
            (&again[1].name, &again[1].features),
            (&lib[1].name, &lib[1].features)
        );
        assert!(again[1].drt.is_nan());
        assert_eq!(
            again[1].integrals.as_ref().unwrap()[1].to_bits(),
            (-0.0f64).to_bits()
        );
        assert_eq!(again[2..], lib[2..]);
    }

    #[test]
    fn mdlib_edited_library_reopens_with_its_values() {
        let mut lib = parse_msl(MSL).unwrap();
        lib.remove(1);
        let g = &mut lib[0];
        g.name = "Glycine 3TMS".into();
        g.spectrum.ri = 1305.5;
        g.quant_ions = vec![174.0, 248.0];
        g.features.insert("FRAG".into(), "(174/2)(248/2)".into());
        g.features
            .insert("COMMENT".into(), "edited\nin the app".into());
        lib.push(LibraryCompound::new("New compound"));
        let path = std::env::temp_dir().join(format!("md-core-{}.mdlib", std::process::id()));
        save_mdlib(&path, &lib).unwrap();
        let again = read_mdlib(&path);
        let mut tmp = path.clone().into_os_string();
        tmp.push(".tmp");
        let tmp_left = Path::new(&tmp).exists();
        std::fs::remove_file(&path).unwrap();
        assert_eq!(again.unwrap(), lib);
        assert!(!tmp_left);
    }

    #[test]
    #[allow(
        clippy::cognitive_complexity,
        reason = "a test of many cases: each assert counts, and Sonar skips tests"
    )]
    fn mdlib_rejects_broken_files() {
        let h = "MetaboliteDetector library 1\n";
        let line = |body: &str| match parse_mdlib(&format!("{h}{body}")) {
            Err(LibraryError::Mdlib { line, .. }) => line,
            r => panic!("{body:?} gave {r:?}"),
        };
        assert_eq!(line("compound\nname a\nbogus 1\nend\n"), 4);
        assert_eq!(line("compound\nname a\nname b\nend\n"), 4);
        assert_eq!(line("compound\nri 1\nend\n"), 2);
        assert_eq!(line("compound\nname a\n"), 2);
        assert_eq!(line("name a\n"), 2);
        assert_eq!(line("compound\nname a\ncompound\n"), 4);
        assert_eq!(line("compound\nname a\\x\nend\n"), 3);
        assert_eq!(line("compound\nname a\\\nend\n"), 3);
        assert_eq!(line("compound\nname a\nri x\nend\n"), 4);
        assert_eq!(line("compound\nname a\npeaks 73:1 74\nend\n"), 4);
        assert_eq!(line("compound\nname a\nfeature K\nend\n"), 4);
        assert_eq!(
            line("compound\nname a\nfeature K\tv\nfeature K\tw\nend\n"),
            5
        );
        assert_eq!(line("compound\nname a\nblock 00\nend\n"), 4);
        assert_eq!(line("compound\nname a\npeak 1 2\nend\n"), 4);
        assert_eq!(line("compound\nname a\nmodel 1 2 3\nend\n"), 4);
        assert_eq!(line("compound\nname a\nmid_size 1 2\nend\n"), 4);
        assert_eq!(line("compound\nname a\nformula x 2 N\nend\n"), 4);
        assert_eq!(
            line(&format!(
                "compound\nname a\nblock +{}\nend\n",
                "0".repeat(79)
            )),
            4
        );
        assert!(matches!(
            parse_mdlib("NAME: x\n"),
            Err(LibraryError::Mdlib { line: 1, .. })
        ));
        assert!(matches!(parse_mdlib(""), Err(LibraryError::Mdlib { .. })));
        assert!(matches!(
            parse_mdlib("MetaboliteDetector library 3\n"),
            Err(LibraryError::Unsupported(_))
        ));
        // Comments, blank lines and CRLF are fine; missing keys take a new compound's values.
        let lib = parse_mdlib(&format!(
            "\u{feff}# x\r\n{h}\r\ncompound\r\nname a\r\nquant_ions 9 7 9\r\nend\r\n"
        ))
        .unwrap();
        let mut want = LibraryCompound::new("a");
        want.quant_ions = vec![7.0, 9.0];
        assert_eq!(lib, [want]);
    }

    #[test]
    fn nist_text_fields_and_import_rules() {
        let text = "mol block line\r\n  > <NAME> \r\n Two\r\n\r\n  lines \n> <SYNONYMS>\nsyn a\nsyn b\n\
            > <MASS SPECTRAL PEAKS>\n69 100\n70 0\n71.9 5\n80 7.5\nbroken\n$$$$\n\
            > <NAME>\nno peaks\n> <MASS SPECTRAL PEAKS>\n50 1\n$$$$\n> <NAME>\nlast\n\
            > <MASS SPECTRAL PEAKS>\n100 1\n";
        let lib = parse_nist(text).unwrap();
        let [a, b] = &lib[..] else { panic!("{lib:?}") };
        assert_eq!(a.name, "Two\nlines");
        assert_eq!(a.features["SYNONYMS"], "syn a\nsyn b");
        assert_eq!(a.spectrum.masses, [71.0, 80.0]);
        assert_eq!(a.spectrum.values, [5, 7]);
        assert_eq!(a.comb_spec, [(71.0, 5), (80.0, 7)]);
        assert_eq!((a.spectrum.ri, a.spectrum.rt_ms), (-1.0, -1.0));
        // The last record needs no `$$$$`.
        assert_eq!(b.name, "last");
        let bom = parse_nist("\u{feff}> <NAME>\nx\n> <MASS SPECTRAL PEAKS>\n71 1\n").unwrap();
        assert_eq!(bom[0].name, "x");
        assert!(matches!(
            parse_nist("> <NAME>\nx\n$$$$\n"),
            Err(LibraryError::Unsupported(_))
        ));
    }

    #[test]
    fn msp_peak_lines_and_repeated_keys() {
        let msp = "Name: Synthetic msp\nSynon: one\nSynon: two\nCAS#: 0-00-0\nNum Peaks: 5\n\
            71 100; 72 12;\n91\t80 \"C7H7+; x\"\n92 15\n105 60.4\n\nName: Next\nNum Peaks: 0\n";
        let lib = parse_msl(msp).unwrap();
        let [c, _] = &lib[..] else { panic!() };
        assert_eq!(c.spectrum.masses, [71.0, 72.0, 91.0, 92.0, 105.0]);
        assert_eq!(c.spectrum.values, [100, 12, 80, 15, 60]);
        assert_eq!(c.features["SYNON"], "one\ntwo");
        assert_eq!(c.features["CAS#"], "0-00-0");
        let line = |text: &str| match parse_msl(text) {
            Err(LibraryError::Msl { line, .. }) => line,
            r => panic!("{text:?} gave {r:?}"),
        };
        assert_eq!(line("Name: a\nNum Peaks: 1\n71\n"), 3);
        assert_eq!(line("Name: a\nNum Peaks: 1\n71 1 \"open\n"), 3);
        // Peak lines before NUM PEAKS are no peaks.
        assert_eq!(line("Name: a\n71 1\n"), 2);
    }
}
