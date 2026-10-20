//! Saved state (S4, #253). eframe keeps everything in one file, `app.ron` in the config dir: egui's
//! memory (panel sizes, tool window positions), the window geometry and, under [`KEY`], our own
//! `KEY=value` lines: the analysis settings under the old QSettings names, the shown panels, the
//! libraries and the recent files. Unknown keys are kept as they are, so a newer version's values
//! survive; `mode` (S7) and `language` (S8) are reserved for those. The first start imports the old
//! 3.3 settings once (`~/.config/TU-Braunschweig/MetaboliteDetector.conf`, on Windows the registry).
//! Every value is checked; a bad one keeps its default, and the unreadable file is copied aside.

use std::io::Read;
use std::path::{Path, PathBuf};

use eframe::egui;
use md_core::peaks::Detector;
use md_core::settings::{GcmsSettings, MAX_LEN};

use crate::App;

/// Storage key of our lines in `app.ron`.
pub const KEY: &str = "metabolite-detector";
const RECENT: usize = 10;

/// What the app keeps between starts besides settings and panels.
#[derive(Default)]
pub struct Saved {
    /// Last opened runs, newest first.
    recent: Vec<PathBuf>,
    /// Libraries opened at start; kept while they load, so a quick exit does not forget them.
    library: Option<PathBuf>,
    nist: Option<PathBuf>,
    /// Keys this version does not know, written back unchanged.
    extra: Vec<(String, String)>,
}

impl Saved {
    pub fn remember(&mut self, path: PathBuf) {
        self.recent.retain(|p| p != &path);
        self.recent.insert(0, path);
        self.recent.truncate(RECENT);
    }
}

/// (section, key, value) with the value already unescaped.
type Entry = (String, String, String);

fn env(k: &str) -> Option<PathBuf> {
    std::env::var_os(k)
        .map(PathBuf::from)
        .filter(|p| p.is_absolute())
}

/// `$XDG_CONFIG_HOME`, else `~/.config`.
fn xdg() -> Option<PathBuf> {
    env("XDG_CONFIG_HOME").or_else(|| Some(env("HOME")?.join(".config")))
}

/// `%APPDATA%\MetaboliteDetector`, else `$XDG_CONFIG_HOME` or `~/.config` + `/metabolite-detector`.
pub fn config_dir() -> Option<PathBuf> {
    if cfg!(windows) {
        Some(env("APPDATA")?.join("MetaboliteDetector"))
    } else {
        Some(xdg()?.join("metabolite-detector"))
    }
}

/// A field of [`GcmsSettings`] as read from and written to a line.
enum Field<'a> {
    F64(&'a mut f64),
    F32(&'a mut f32),
    I32(&'a mut i32),
    Bool(&'a mut bool),
    Ranges(&'a mut Vec<(f32, f32)>),
    Detector(&'a mut Detector),
}

type Access = fn(&mut GcmsSettings) -> Field<'_>;

/// Every setting under its old key; reading and writing both go through this one table.
const FIELDS: [(&str, Access); 34] = [
    ("SF_PEAK_THRESHOLD_BEGIN", |s| {
        Field::F64(&mut s.sf_peak_threshold_begin)
    }),
    ("SF_PEAK_THRESHOLD_END", |s| {
        Field::F64(&mut s.sf_peak_threshold_end)
    }),
    ("SF_MAX_BASE_DEV", |s| Field::F64(&mut s.sf_max_base_dev)),
    ("SF_FWHM", |s| Field::F64(&mut s.sf_fwhm)),
    ("SF_CENTROID_THRESHOLD_BEGIN", |s| {
        Field::F64(&mut s.sf_centroid_threshold_begin)
    }),
    ("SF_CENTROID_THRESHOLD_END", |s| {
        Field::F64(&mut s.sf_centroid_threshold_end)
    }),
    ("SF_CENTROID_MAX_BASE_DEV", |s| {
        Field::F64(&mut s.sf_centroid_max_base_dev)
    }),
    ("AN_PEAK_THRESHOLD_BEGIN", |s| {
        Field::F64(&mut s.an_peak_threshold_begin)
    }),
    ("AN_PEAK_THRESHOLD_END", |s| {
        Field::F64(&mut s.an_peak_threshold_end)
    }),
    ("AN_MAX_BASE_DEV", |s| Field::F64(&mut s.an_max_base_dev)),
    ("AN_MIN_PEAK_HEIGHT", |s| {
        Field::F64(&mut s.an_min_peak_height)
    }),
    ("AN_DECONVOLUTION_WIDTH", |s| {
        Field::F64(&mut s.an_deconvolution_width)
    }),
    ("AN_NUMBER_BINS", |s| Field::I32(&mut s.an_number_bins)),
    ("AN_MIN_PEAK_NUMBER", |s| {
        Field::I32(&mut s.an_min_peak_number)
    }),
    ("AN_MIN_HEIGHT_BP", |s| Field::F64(&mut s.an_min_height_bp)),
    ("AN_MODELED_VALUES_FOR_INTEGRATION", |s| {
        Field::Bool(&mut s.an_modeled_values_for_integration)
    }),
    ("LS_RI_DIFF", |s| Field::F64(&mut s.ls_ri_diff)),
    ("LS_PURE_FACTOR", |s| Field::F64(&mut s.ls_pure_factor)),
    ("LS_IMPURE_FACTOR", |s| Field::F64(&mut s.ls_impure_factor)),
    ("LS_SCALE_FACTOR_ION", |s| {
        Field::F64(&mut s.ls_scale_factor_ion)
    }),
    ("LS_SCALE_FACTOR_VALUE", |s| {
        Field::F64(&mut s.ls_scale_factor_value)
    }),
    ("LS_SCALE_FACTOR_NO_MATCH", |s| {
        Field::F64(&mut s.ls_scale_factor_no_match)
    }),
    ("LS_RESOLUTION", |s| Field::F32(&mut s.ls_resolution)),
    ("LS_USE_RI", |s| Field::Bool(&mut s.ls_use_ri)),
    ("LS_USE_SCALED_LIB", |s| {
        Field::Bool(&mut s.ls_use_scaled_lib)
    }),
    ("LS_USE_SPEC", |s| Field::Bool(&mut s.ls_use_spec)),
    ("LS_MASS_FILTER", |s| Field::Ranges(&mut s.ls_mass_filter)),
    ("QUANT_ION_NUMBER", |s| Field::I32(&mut s.quant_ion_number)),
    ("QUANT_MIN_DIST_IONS", |s| {
        Field::F64(&mut s.quant_min_dist_ions)
    }),
    ("QUANT_QUALITY_CUTOFF", |s| {
        Field::F64(&mut s.quant_quality_cutoff)
    }),
    ("QUANT_EXCL_IONS", |s| Field::Ranges(&mut s.quant_excl_ions)),
    ("SQL_MAX_NUMBER_FRAGMENTS", |s| {
        Field::I32(&mut s.sql_max_number_fragments)
    }),
    ("SQL_MINIMAL_NUMBER_EQUAL_FRAGMENTS", |s| {
        Field::I32(&mut s.sql_minimal_number_equal_fragments)
    }),
    ("PEAK_DETECTOR", |s| Field::Detector(&mut s.peak_detector)),
];

fn put<T: std::str::FromStr>(f: &mut T, v: &str) -> bool {
    v.parse().map(|x| *f = x).is_ok()
}

impl Field<'_> {
    fn parse(self, v: &str) -> bool {
        match self {
            Field::F64(f) => put(f, v),
            Field::F32(f) => put(f, v) || variant_f32(v).map(|x| *f = x).is_some(),
            Field::I32(f) => put(f, v),
            Field::Bool(f) => put(f, v),
            Field::Ranges(f) => ranges(v).map(|r| *f = r).is_some(),
            Field::Detector(f) => detector(v).map(|d| *f = d).is_some(),
        }
    }

    fn text(self) -> String {
        match self {
            Field::F64(f) => f.to_string(),
            Field::F32(f) => f.to_string(),
            Field::I32(f) => f.to_string(),
            Field::Bool(f) => f.to_string(),
            Field::Ranges(f) => (f.iter().map(|(a, b)| format!("{a}:{b}")))
                .collect::<Vec<_>>()
                .join(";"),
            Field::Detector(Detector::Default) => "default".into(),
            Field::Detector(Detector::Matched { width, threshold }) => {
                format!("matched {width} {threshold}")
            }
        }
    }
}

/// `73:75;147:148`; the old app's `@Invalid()` is its empty list.
fn ranges(v: &str) -> Option<Vec<(f32, f32)>> {
    if v.is_empty() || v == "@Invalid()" {
        return Some(Vec::new());
    }
    // shortcut: a non-empty old list (a QVariantList, never seen) is not decoded and keeps the default.
    (v.split(';'))
        .map(|r| {
            let (a, b) = r.split_once(':')?;
            Some((a.parse().ok()?, b.parse().ok()?))
        })
        .collect()
}

fn detector(v: &str) -> Option<Detector> {
    match v.split_whitespace().collect::<Vec<_>>()[..] {
        ["default"] => Some(Detector::Default),
        ["matched", w, t] => Some(Detector::Matched {
            width: w.parse().ok()?,
            threshold: t.parse().ok()?,
        }),
        _ => None,
    }
}

/// The old `LS_RESOLUTION=@Variant(\0\0\0\x87>\x99\x99\x9a)`: a QDataStream QVariant of type 135
/// (float), big-endian, here already unescaped to one char per byte.
fn variant_f32(v: &str) -> Option<f32> {
    let inner = v.strip_prefix("@Variant(")?.strip_suffix(')')?;
    let bytes: Vec<u8> = inner
        .chars()
        .map(|c| u8::try_from(u32::from(c)).ok())
        .collect::<Option<_>>()?;
    match bytes[..] {
        [0, 0, 0, 135, a, b, c, d] => Some(f32::from_be_bytes([a, b, c, d])),
        _ => None,
    }
}

/// Keys the old app named differently from its settings file (Inferred: same meaning and default).
fn alias(key: &str) -> &str {
    match key {
        "AN_BP_HEIGHT" => "AN_MIN_HEIGHT_BP",
        "MASS_FILTER" => "LS_MASS_FILTER",
        k => k,
    }
}

/// Our own lines: `key=value`, no sections, no escaping.
fn own_entries(text: &str) -> Vec<Entry> {
    (text.lines())
        .filter_map(|l| l.split_once('='))
        .map(|(k, v)| (String::new(), k.into(), v.into()))
        .collect()
}

/// An old QSettings INI file: `[section]` lines and `key=value` with Qt's quoting and escapes.
fn ini_entries(text: &str) -> Vec<Entry> {
    let mut section = String::new();
    let mut out = Vec::new();
    for line in text.lines().map(str::trim) {
        if let Some(s) = line.strip_prefix('[').and_then(|l| l.strip_suffix(']')) {
            section = s.into();
        } else if let Some((k, v)) = line.split_once('=') {
            out.push((section.clone(), k.trim().into(), unescape(v.trim())));
        }
    }
    out
}

/// Qt's INI value escapes: `"…"` quoting, `\xHH…` hex, `\0…` octal and the C escapes.
fn unescape(v: &str) -> String {
    fn number(it: &mut std::iter::Peekable<std::str::Chars<'_>>, radix: u32, mut n: u32) -> char {
        while let Some(d) = it.peek().and_then(|c| c.to_digit(radix)) {
            n = n.saturating_mul(radix).saturating_add(d);
            it.next();
        }
        char::from_u32(n).unwrap_or(char::REPLACEMENT_CHARACTER)
    }
    let mut out = String::new();
    let mut it = v.chars().peekable();
    while let Some(c) = it.next() {
        match c {
            '"' => {}
            '\\' => match it.next() {
                Some('x') => out.push(number(&mut it, 16, 0)),
                Some(d @ '0'..='7') => out.push(number(&mut it, 8, d.to_digit(8).unwrap_or(0))),
                Some('a') => out.push('\x07'),
                Some('b') => out.push('\x08'),
                Some('f') => out.push('\x0c'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some('t') => out.push('\t'),
                Some('v') => out.push('\x0b'),
                Some(c) => out.push(c),
                None => {}
            },
            c => out.push(c),
        }
    }
    out
}

const REG_ROOT: &str = r"HKEY_CURRENT_USER\Software\TU-Braunschweig\MetaboliteDetector";

/// `reg query <root> /s` output. Qt writes ints as REG_DWORD, strings with a NUL (the `@Variant`)
/// as REG_BINARY UTF-16, everything else as REG_SZ, all unescaped.
#[cfg_attr(not(windows), allow(dead_code))]
fn reg_entries(text: &str) -> Vec<Entry> {
    let mut section = String::new();
    let mut out = Vec::new();
    for line in text.lines() {
        if let Some(s) = line.strip_prefix(REG_ROOT) {
            section = s.trim_start_matches('\\').into();
            continue;
        }
        let Some(rest) = line.strip_prefix("    ") else {
            continue;
        };
        let mut parts = rest.splitn(3, "    ");
        let (Some(name), Some(kind)) = (parts.next(), parts.next()) else {
            continue;
        };
        let data = parts.next().unwrap_or("");
        let hex = |d: &str| d.trim_start_matches("0x").to_owned();
        let value = match kind {
            "REG_SZ" => data.into(),
            "REG_DWORD" => u32::from_str_radix(&hex(data), 16).map_or(String::new(), |n| {
                i32::from_ne_bytes(n.to_ne_bytes()).to_string()
            }),
            "REG_QWORD" => u64::from_str_radix(&hex(data), 16).map_or(String::new(), |n| {
                i64::from_ne_bytes(n.to_ne_bytes()).to_string()
            }),
            "REG_BINARY" => {
                let bytes: Vec<u8> = (0..data.len() / 2)
                    .filter_map(|i| u8::from_str_radix(data.get(2 * i..2 * i + 2)?, 16).ok())
                    .collect();
                let units: Vec<u16> = (bytes.as_chunks::<2>().0.iter())
                    .map(|b| u16::from_le_bytes(*b))
                    .collect();
                String::from_utf16_lossy(&units)
            }
            _ => continue,
        };
        out.push((section.clone(), name.into(), value));
    }
    out
}

/// The old app's settings and where they came from, if it ever ran here.
fn old_entries() -> Option<(String, Vec<Entry>)> {
    registry().or_else(|| {
        let path = xdg()?.join("TU-Braunschweig/MetaboliteDetector.conf");
        let text = read_capped(&path)?;
        Some((path.display().to_string(), ini_entries(&text)))
    })
}

#[cfg(windows)]
fn registry() -> Option<(String, Vec<Entry>)> {
    use std::os::windows::process::CommandExt;
    // CREATE_NO_WINDOW: no console flashes up next to the app.
    let out = std::process::Command::new("reg")
        .args(["query", REG_ROOT, "/s"])
        .creation_flags(0x0800_0000)
        .output()
        .ok()?;
    // shortcut: reg.exe prints in the OEM code page, so a non-ASCII path arrives garbled and is
    // skipped at start as missing; read the registry directly if that matters.
    let text = String::from_utf8_lossy(&out.stdout);
    (out.status.success()).then(|| (format!("the registry ({REG_ROOT})"), reg_entries(&text)))
}

#[cfg(not(windows))]
fn registry() -> Option<(String, Vec<Entry>)> {
    None
}

fn read_capped(path: &Path) -> Option<String> {
    let mut bytes = Vec::new();
    std::fs::File::open(path)
        .and_then(|f| f.take(MAX_LEN).read_to_end(&mut bytes))
        .ok()?;
    Some(String::from_utf8_lossy(&bytes).into_owned())
}

/// The old `size=@Size(1591 800)`.
fn old_size(v: &str) -> Option<egui::Vec2> {
    let inner = v.strip_prefix("@Size(")?.strip_suffix(')')?;
    let (w, h) = inner.split_once(' ')?;
    let ok = |x: f32| (200.0..=20_000.0).contains(&x);
    let (w, h) = (
        w.parse().ok().filter(|x| ok(*x))?,
        h.parse().ok().filter(|x| ok(*x))?,
    );
    Some(egui::vec2(w, h))
}

fn path(v: &str) -> Option<PathBuf> {
    (!v.is_empty()).then(|| v.into())
}

/// Applies saved lines to a fresh app. Returns the keys whose values were bad (they keep the default).
fn apply(app: &mut App, ctx: &egui::Context, entries: Vec<Entry>) -> Vec<String> {
    let mut s = app.settings.applied.clone();
    let mut bad = Vec::new();
    let mut recent = Vec::new();
    for (section, key, v) in entries {
        let key = alias(&key).to_owned();
        let ok = if let Some((_, field)) = FIELDS.iter().find(|(k, _)| *k == key) {
            field(&mut s).parse(&v)
        } else if let Some(name) = key.strip_prefix("panel.")
            && let Some((_, on)) = app.panels.saved().into_iter().find(|(k, _)| *k == name)
        {
            put(on, &v)
        } else {
            match (section.as_str(), key.as_str()) {
                (_, "library" | "rc_lib_file") => {
                    app.saved.library = path(&v);
                    true
                }
                (_, "nist_library" | "rc_lib_nist_file") => {
                    app.saved.nist = path(&v);
                    true
                }
                (_, "recent" | "rc_open_file") => {
                    recent.extend(path(&v));
                    true
                }
                ("GUI", "size") => old_size(&v)
                    .map(|size| ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(size)))
                    .is_some(),
                // The old app's other settings stay for later versions (e.g. LS_CUTOFF for S5);
                // its window state blobs do not.
                ("" | "Algorithm", _) => {
                    app.saved.extra.push((key.clone(), v));
                    true
                }
                _ => true,
            }
        };
        if !ok {
            bad.push(key);
        }
    }
    for p in recent.into_iter().rev() {
        app.saved.remember(p);
    }
    if let Err(e) = s.validate() {
        bad.push(e);
        s = GcmsSettings::default();
    }
    app.settings.restore(s);
    bad
}

/// The lines saved under [`KEY`].
pub fn text(app: &App) -> String {
    let mut s = app.settings.applied.clone();
    let mut lines: Vec<(String, String)> = (FIELDS.iter())
        .map(|(k, field)| ((*k).to_owned(), field(&mut s).text()))
        .collect();
    let mut panels = app.panels;
    lines.extend((panels.saved()).map(|(k, on)| (format!("panel.{k}"), on.to_string())));
    let library = (app.library.as_ref()).map(|l| l.path.clone());
    let nist = (app.nist.lib.as_ref()).map(|l| l.path.clone());
    let paths = [
        ("library", library.or_else(|| app.saved.library.clone())),
        ("nist_library", nist.or_else(|| app.saved.nist.clone())),
    ];
    let paths = paths.into_iter().filter_map(|(k, p)| Some((k, p?)));
    let recent = app.saved.recent.iter().map(|p| ("recent", p.clone()));
    // shortcut: a path that is not UTF-8 is not saved; store it as bytes if that ever matters.
    lines.extend(
        (paths.chain(recent)).filter_map(|(k, p)| Some((k.to_owned(), p.to_str()?.to_owned()))),
    );
    lines.extend(app.saved.extra.iter().cloned());
    // A line break would split a value into a second, wrong line.
    (lines.into_iter())
        .filter(|(_, v)| !v.contains(['\n', '\r']))
        .map(|(k, v)| format!("{k}={v}\n"))
        .collect()
}

/// Copies an unreadable `app.ron` aside before eframe overwrites it; returns what to tell the user.
fn backup(file: Option<&Path>) -> String {
    let Some(file) = file.filter(|f| f.exists()) else {
        return String::new();
    };
    let secs = (std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH))
        .map_or(0, |d| d.as_secs());
    let copy = file.with_extension(format!("ron.{secs}.bak"));
    match std::fs::copy(file, &copy) {
        Ok(_) => format!(" The old file is kept as {}.", copy.display()),
        Err(e) => format!(" Copying it aside failed: {e}."),
    }
}

/// At start: the saved state, else on the first start the old app's settings; then the libraries.
pub fn restore(
    app: &mut App,
    ctx: &egui::Context,
    storage: Option<&dyn eframe::Storage>,
    file: Option<&Path>,
) {
    let saved = storage.and_then(|s| s.get_string(KEY));
    let note = match saved {
        Some(text) => {
            let bad = apply(app, ctx, own_entries(&text));
            (!bad.is_empty()).then(|| {
                Err(format!(
                    "Saved values could not be read and use their defaults: {}.{}",
                    bad.join(", "),
                    backup(file)
                ))
            })
        }
        // eframe drops a file it cannot parse and would write over it.
        None if file.is_some_and(Path::exists) => Some(Err(format!(
            "The saved settings could not be read; the defaults are used.{}",
            backup(file)
        ))),
        None => old_entries().map(|(source, entries)| {
            let bad = apply(app, ctx, entries);
            let done = format!("Imported the settings of MetaboliteDetector 3.3 from {source}.");
            if bad.is_empty() {
                Ok(done)
            } else {
                Err(format!(
                    "{done} These keep their defaults: {}.",
                    bad.join(", ")
                ))
            }
        }),
    };
    if let Some(note) = note {
        app.settings.show_note(note);
    }
    // A library that is gone is forgotten without a message; one that fails to load says so.
    if let Some(p) = app.saved.library.clone().filter(|p| p.exists()) {
        app.open_library(ctx, p);
    } else {
        app.saved.library = None;
    }
    if let Some(p) = app.saved.nist.clone().filter(|p| p.exists()) {
        app.open_nist(ctx, p);
    } else {
        app.saved.nist = None;
    }
}

impl App {
    /// File → Open recent: the last opened runs, newest first.
    pub fn recent_menu(&mut self, ui: &mut egui::Ui) {
        let empty = self.saved.recent.is_empty();
        ui.add_enabled_ui(!empty, |ui| {
            ui.menu_button("Open recent", |ui| {
                for p in self.saved.recent.clone() {
                    let name = p.file_name().unwrap_or(p.as_os_str()).to_string_lossy();
                    if ui
                        .button(name)
                        .on_hover_text(p.display().to_string())
                        .clicked()
                    {
                        self.open(ui.ctx(), p);
                    }
                }
            })
            .response
            .on_disabled_hover_text("No runs opened yet");
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply_text(text: &str, entries: fn(&str) -> Vec<Entry>) -> (App, Vec<String>) {
        let mut app = App::default();
        let bad = apply(&mut app, &egui::Context::default(), entries(text));
        (app, bad)
    }

    /// Lines in the form the old 3.3 wrote them (`[Algorithm]` of a real file, `[GUI]` shortened).
    const OLD: &str = r#"[Algorithm]
AN_BP_HEIGHT=0
AN_DECONVOLUTION_WIDTH=1
AN_PEAK_THRESHOLD_END=-5
LS_CUTOFF=0.85
LS_RESOLUTION=@Variant(\0\0\0\x87>\x99\x99\x9a)
LS_RI_DIFF=10
LS_USE_SCALED_LIB=true
MASS_FILTER=@Invalid()
QUANT_EXCL_IONS=@Invalid()
QUANT_ION_NUMBER=3
rc_lib_file=
rc_open_file=/data/run 1.bin

[GUI]
mw_state="@ByteArray(\0\0\0\xff\0\0\0\0\xfd)"
size=@Size(1591 800)
"#;

    #[test]
    fn old_conf_values_win_over_defaults() {
        let (app, bad) = apply_text(OLD, ini_entries);
        assert!(bad.is_empty(), "{bad:?}");
        let s = &app.settings.applied;
        assert_eq!(
            s.ls_ri_diff, 10.0,
            "a saved LS_RI_DIFF beats the default 20"
        );
        assert_eq!(s.ls_resolution, 0.3);
        assert_eq!(s.an_peak_threshold_end, -5.0);
        assert_eq!(app.saved.recent, [PathBuf::from("/data/run 1.bin")]);
        assert_eq!(app.saved.library, None);
        // Kept for later, the GUI blob is not.
        assert_eq!(app.saved.extra, [("LS_CUTOFF".into(), "0.85".into())]);
    }

    /// Reads a real old `MetaboliteDetector.conf` copied into `example-data/old-results/`.
    #[test]
    fn real_old_conf() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../example-data/old-results/MetaboliteDetector.conf");
        let Some(text) = read_capped(&path) else {
            eprintln!("skipped: no {}", path.display());
            return;
        };
        let (app, bad) = apply_text(&text, ini_entries);
        assert!(bad.is_empty(), "{bad:?}");
        assert_eq!(app.settings.applied.ls_ri_diff, 10.0);
        assert_eq!(app.settings.applied.ls_resolution, 0.3);
        assert_eq!(app.saved.recent.len(), 1);
    }

    #[test]
    fn registry_values_read_like_the_ini() {
        // LS_RESOLUTION: "@Variant(\0\0\0\x87>\x80\0\0)" (0.25) as UTF-16LE.
        let reg = format!(
            "\r\n{REG_ROOT}\\Algorithm\r\n    LS_RI_DIFF    REG_SZ    10\r\n    QUANT_ION_NUMBER    REG_DWORD    0x4\r\n    \
             LS_RESOLUTION    REG_BINARY    4000560061007200690061006E00740028000000000000008700\
             3E008000000000002900\r\n    rc_save_file    REG_SZ    \r\n"
        );
        let (app, bad) = apply_text(&reg, reg_entries);
        assert!(bad.is_empty(), "{bad:?}");
        assert_eq!(app.settings.applied.ls_ri_diff, 10.0);
        assert_eq!(app.settings.applied.quant_ion_number, 4);
        assert_eq!(app.settings.applied.ls_resolution, 0.25);
        assert_eq!(app.saved.extra, [("rc_save_file".into(), String::new())]);
    }

    #[test]
    fn own_lines_round_trip_and_keep_unknown_keys() {
        let mut app = App::default();
        let s = &mut app.settings.applied;
        (s.ls_ri_diff, s.ls_use_ri) = (12.5, false);
        s.ls_mass_filter = vec![(73.0, 75.0), (147.0, 148.5)];
        s.peak_detector = Detector::Matched {
            width: 7,
            threshold: -3,
        };
        app.panels.info = false;
        app.saved.remember("C:\\runs\\a=b.cdf".into());
        app.saved.extra.push(("language".into(), "de".into()));
        let text = text(&app);
        let (back, bad) = apply_text(&text, own_entries);
        assert!(bad.is_empty(), "{bad:?}");
        assert_eq!(back.settings.applied, app.settings.applied);
        assert_eq!(back.panels, app.panels);
        assert_eq!(back.saved.recent, app.saved.recent);
        assert_eq!(back.saved.extra, app.saved.extra);
    }

    #[test]
    fn bad_values_keep_defaults() {
        let (app, bad) = apply_text(
            "LS_RI_DIFF=abc\nAN_NUMBER_BINS=0\npanel.runs=yes\n",
            own_entries,
        );
        assert_eq!(
            bad,
            [
                "LS_RI_DIFF",
                "panel.runs",
                "AN_NUMBER_BINS: must be between 1 and 100"
            ]
        );
        assert_eq!(app.settings.applied, GcmsSettings::default());
        assert!(app.panels.runs);
    }

    /// A file eframe could not read is copied aside and the user is told; the defaults apply.
    #[test]
    fn unreadable_file_is_kept_as_backup() {
        let dir = std::env::temp_dir().join(format!("md-persist-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("app.ron");
        std::fs::write(&file, "not ron {").unwrap();
        let mut app = App::default();
        restore(&mut app, &egui::Context::default(), None, Some(&file));
        let backups: Vec<_> = (std::fs::read_dir(&dir).unwrap())
            .filter_map(|e| e.ok()?.file_name().into_string().ok())
            .filter(|n| n.ends_with(".bak"))
            .collect();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(backups.len(), 1, "{backups:?}");
        assert!(app.settings.open);
        assert_eq!(app.settings.applied, GcmsSettings::default());
    }

    #[test]
    fn unescape_qt() {
        assert_eq!(unescape(r#""a,b" \x41\x42\0\aP"#), "a,b AB\0\x07P");
    }
}
