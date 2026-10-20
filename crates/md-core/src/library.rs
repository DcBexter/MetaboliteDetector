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
    /// 40-byte records of unknown meaning (file-formats.md §1.3), kept verbatim.
    pub blocks: Vec<[u8; BLOCK_BYTES]>,
    pub sources: Vec<String>,
}

impl LibraryCompound {
    /// A compound without spectrum, with the defaults the old MSL importer left (sample library).
    pub fn new(name: &str) -> Self {
        Self {
            spectrum: LibrarySpectrum {
                scan: 0,
                rt_ms: 0.0,
                ri: 0.0,
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
        }
    }
}

const BLOCK_BYTES: usize = 40;
/// Size of a record whose lists and strings are all empty: a lower bound per compound.
const MIN_COMPOUND_BYTES: usize = 36 + 56 + 24;
// shortcut: whole library in RAM, capped so a gzip bomb cannot exhaust memory; raise it when real
// libraries come close.
const MAX_LBR_BYTES: u64 = 1 << 30;

/// Reads a `.lbr` file (gzip, or an already decompressed payload).
pub fn read_lbr(path: impl AsRef<Path>) -> Result<Vec<LibraryCompound>, LibraryError> {
    parse_lbr(&std::fs::read(path)?)
}

/// Parses `.lbr` bytes. Only version 2 is decoded; every length and check value is validated.
pub fn parse_lbr(bytes: &[u8]) -> Result<Vec<LibraryCompound>, LibraryError> {
    let unzipped;
    let data = if bytes.starts_with(&[0x1f, 0x8b]) {
        unzipped = gunzip(bytes)?;
        &unzipped[..]
    } else {
        bytes
    };
    let mut r = Reader { data, pos: 0 };
    match r.u32()? {
        2 => {}
        5 => {
            return Err(LibraryError::Unsupported(
                "library version 5: its layout is unknown until a sample file is available".into(),
            ));
        }
        v => return Err(LibraryError::Unsupported(format!("library version {v}"))),
    }
    let compounds = r.list(MIN_COMPOUND_BYTES, read_compound)?;
    if r.remaining() != 0 {
        return r.corrupt(format!("{} bytes after the last compound", r.remaining()));
    }
    Ok(compounds)
}

fn gunzip(bytes: &[u8]) -> Result<Vec<u8>, LibraryError> {
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

fn read_compound(r: &mut Reader) -> Result<LibraryCompound, LibraryError> {
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

    if r.u32()? != 0 {
        return Err(LibraryError::Unsupported(
            "compounds with peak records (layout unknown)".into(),
        ));
    }
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
    if r.u32()? != expected {
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
        comb_spec: r.list(8, |r| Ok((r.f32()?, r.i32()?)))?,
        ion_freq: r.list(8, |r| Ok((r.f32()?, r.u32()?)))?,
        spec_count: r.u32()?,
        ris: r.list(8, Reader::f64)?,
        blocks: r.list(BLOCK_BYTES, Reader::array)?,
        sources: r.list(4, Reader::string)?,
    })
}

/// Little-endian reader that fails instead of panicking on short or inconsistent input.
struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }

    fn corrupt<T>(&self, reason: impl Into<String>) -> Result<T, LibraryError> {
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

    fn array<const N: usize>(&mut self) -> Result<[u8; N], LibraryError> {
        let mut a = [0; N];
        a.copy_from_slice(self.take(N)?);
        Ok(a)
    }

    fn u32(&mut self) -> Result<u32, LibraryError> {
        self.array().map(u32::from_le_bytes)
    }
    fn i32(&mut self) -> Result<i32, LibraryError> {
        self.array().map(i32::from_le_bytes)
    }
    fn i64(&mut self) -> Result<i64, LibraryError> {
        self.array().map(i64::from_le_bytes)
    }
    fn f32(&mut self) -> Result<f32, LibraryError> {
        self.array().map(f32::from_le_bytes)
    }
    fn f64(&mut self) -> Result<f64, LibraryError> {
        self.array().map(f64::from_le_bytes)
    }

    /// Reads a count of elements of at least `min_bytes` each; a count that cannot fit in the
    /// remaining bytes is rejected before anything is allocated.
    fn count(&mut self, n: usize, min_bytes: usize) -> Result<usize, LibraryError> {
        if n.checked_mul(min_bytes)
            .is_none_or(|b| b > self.remaining())
        {
            return self.corrupt(format!("count {n} exceeds the remaining data"));
        }
        Ok(n)
    }

    fn list<T>(
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
    fn f64_list_or_null(&mut self) -> Result<Option<Vec<f64>>, LibraryError> {
        let Ok(n) = usize::try_from(self.i32()?) else {
            return Ok(None);
        };
        let n = self.count(n, 8)?;
        (0..n)
            .map(|_| self.f64())
            .collect::<Result<_, _>>()
            .map(Some)
    }

    fn bytes(&mut self) -> Result<&'a [u8], LibraryError> {
        let n = self.u32()? as usize;
        self.take(n)
    }

    fn string(&mut self) -> Result<String, LibraryError> {
        self.bytes().map(text)
    }
}

/// The old files carry no encoding: UTF-8 when valid, otherwise Latin-1 (every byte is one char).
fn text(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).unwrap_or_else(|_| bytes.iter().map(|&b| b as char).collect())
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

/// Parses AMDIS MSL text. `RT:` is in minutes; keys without a field of their own become features.
/// Fields MSL has no key for keep the defaults the old importer left (sample library).
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
                c.features.insert(key, value.to_string());
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
        let mz: f32 = mz.parse().map_err(|_| "invalid m/z")?;
        let v: f64 = v.parse().map_err(|_| "invalid intensity")?;
        if !mz.is_finite() || !(0.0..=i32::MAX as f64).contains(&v) {
            return Err("m/z or intensity out of range");
        }
        s.masses.push(mz);
        s.values.push(v.round() as i32);
    }
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
    // The old importer copies the spectrum into the combined spectrum (sample library).
    let s = &c.spectrum;
    c.comb_spec = s
        .masses
        .iter()
        .copied()
        .zip(s.values.iter().copied())
        .collect();
    c.ion_freq = s.masses.iter().map(|&m| (m, 1)).collect();
    Ok(Some(c))
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

/// First line of a `.mdlib` file; the number is the format version.
const MDLIB_HEADER: &str = "MetaboliteDetector library 1";

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
        Some((_, l)) if l.trim_end() == MDLIB_HEADER => {}
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
                if !matches!(key, "feature" | "block" | "source") {
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

fn set_mdlib_field(c: &mut LibraryCompound, key: &str, v: &str) -> Result<(), String> {
    fn num<T: std::str::FromStr>(v: &str) -> Result<T, String> {
        v.parse().map_err(|_| format!("invalid number `{v}`"))
    }
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
        _ => return Err(format!("unknown key `{key}`")),
    }
    Ok(())
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
    use super::*;

    /// Version 2 payload with one compound, written like the old `toStream`.
    fn sample_payload() -> Vec<u8> {
        let mut b = Vec::new();
        let u32 = |b: &mut Vec<u8>, v: u32| b.extend(v.to_le_bytes());
        let string = |b: &mut Vec<u8>, s: &str| {
            b.extend((s.len() as u32).to_le_bytes());
            b.extend(s.as_bytes());
        };
        u32(&mut b, 2);
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
        // Compound: no peaks, name, null fwhm/integrals, 1 quant ion, no auto ions, 1 feature.
        u32(&mut b, 0);
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
        u32(&mut b, 1 + 1 + 11);
        // comb_spec, ion_freq, spec_count, ris, blocks, sources.
        u32(&mut b, 1);
        b.extend(73.0f32.to_le_bytes());
        b.extend(1000i32.to_le_bytes());
        u32(&mut b, 1);
        b.extend(73.0f32.to_le_bytes());
        u32(&mut b, 1);
        u32(&mut b, 1);
        u32(&mut b, 1);
        b.extend(1369.0f64.to_le_bytes());
        u32(&mut b, 1);
        b.extend([7; BLOCK_BYTES]);
        u32(&mut b, 1);
        string(&mut b, "Initial import from/home/kh/Desktop/test.msl");
        b
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
        assert!(set(0, 5).contains("version 5"));
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
        assert_eq!(again[2], lib[2]);
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
            parse_mdlib("MetaboliteDetector library 2\n"),
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
}
