//! Saved state (S4, #253). eframe keeps everything in one file, `app.ron` in the config dir: egui's
//! memory (panel sizes, tool window positions), the window geometry and, under [`KEY`], our own
//! `KEY=value` lines: the analysis settings under the old QSettings names, the shown panels, the
//! libraries and the recent files. Unknown keys are kept as they are, so a newer version's values
//! survive. `language` is `en` or `de` (S8), saved once the user picks one other than the system's. `mode` is the guided or expert mode (S7); without it
//! the start asks for it. The first start imports the old
//! 3.3 settings once (`~/.config/TU-Braunschweig/MetaboliteDetector.conf`, on Windows the registry).
//! Every value is checked; a bad one keeps its default. A saved deconvolution width 1 becomes Auto once
//! (#311, marker `auto_width_once`).
//!
//! eframe writes `app.ron` itself, not atomically, so [`guard`] runs before it reads the file and
//! keeps the last good copy as `app.ron.bak` (the only backup, rewritten via temp file + rename).
//! A damaged file is replaced by that copy (#286, #288). Two open windows: the last one to save wins
//! for settings, panels and libraries; the recent files are merged with the file on disk (#287).

use std::collections::HashMap;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use eframe::egui;
use md_core::settings::{MAX_LEN, keys};

use crate::App;
use crate::lang::{Msg, msg, t};

/// Storage key of our lines in `app.ron`.
pub const KEY: &str = "metabolite-detector";
const RECENT: usize = 10;
/// Marker line: the saved width 1 was switched to Auto once (#311).
const WIDTH_ONCE: &str = "auto_width_once";

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
    /// The `app.ron` eframe writes, to merge with at save time.
    file: Option<PathBuf>,
    /// The one-time width 1 → Auto switch ran (#311); saved as [`WIDTH_ONCE`].
    width_once: bool,
    /// A language was saved; without one the app follows the system language until the user picks another.
    language: bool,
}

impl Saved {
    /// Adds the recent files another window saved since we started (only those that still exist,
    /// so a file we dropped does not come back); ours stay first.
    pub fn merge_disk(&mut self) {
        let Some(state) = self.file.as_deref().and_then(read_state) else {
            return;
        };
        let disk = own_entries(state.get(KEY).map_or("", String::as_str));
        let mut new: Vec<PathBuf> = (disk.into_iter())
            .filter_map(|(_, k, v)| (k == "recent").then(|| path(&v)).flatten())
            .filter(|p| p.exists() && !self.recent.contains(p))
            .collect();
        self.recent.append(&mut new);
        self.recent.truncate(RECENT);
    }

    pub fn remember(&mut self, path: PathBuf) {
        // Absolute, so a command-line `a.cdf` and the same file from the dialog are one entry.
        let path = std::path::absolute(&path).unwrap_or(path);
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

fn put<T: std::str::FromStr>(f: &mut T, v: &str) -> bool {
    v.parse().map(|x| *f = x).is_ok()
}

/// Keys the old app named differently from its settings file (Inferred: same meaning and default).
fn alias(key: &str) -> &str {
    match key {
        keys::LEGACY_AN_BP_HEIGHT => keys::AN_MIN_HEIGHT_BP,
        keys::LEGACY_MASS_FILTER => keys::LS_MASS_FILTER,
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
fn old_entries() -> Option<(Msg, Vec<Entry>)> {
    registry().or_else(|| {
        let path = xdg()?.join("TU-Braunschweig/MetaboliteDetector.conf");
        let text = read_capped(&path)?;
        Some((path.display().to_string().into(), ini_entries(&text)))
    })
}

#[cfg(windows)]
fn registry() -> Option<(Msg, Vec<Entry>)> {
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
    (out.status.success()).then(|| {
        (
            msg!("the registry ({REG_ROOT})", "der Registry ({REG_ROOT})"),
            reg_entries(&text),
        )
    })
}

#[cfg(not(windows))]
fn registry() -> Option<(Msg, Vec<Entry>)> {
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
    let mut s = app.settings.applied().clone();
    let mut bad = Vec::new();
    let mut recent = Vec::new();
    for (section, key, v) in entries {
        let key = alias(&key).to_owned();
        let ok = if let Some(ok) = s.set(&key, &v) {
            ok
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
                ("", "mode") => {
                    match crate::guided::Mode::from_key(&v) {
                        Some(m) => (app.guided.mode, app.guided.ask) = (m, false),
                        // Maybe from a newer version: kept until the user answers, so not "unreadable".
                        None => app.saved.extra.push((key.clone(), v)),
                    }
                    true
                }
                ("", "language") => {
                    match crate::lang::Lang::from_key(&v) {
                        Some(l) => {
                            crate::lang::set(l);
                            app.saved.language = true;
                        }
                        // Maybe from a newer version: kept, the system language is used meanwhile.
                        None => app.saved.extra.push((key.clone(), v)),
                    }
                    true
                }
                ("", WIDTH_ONCE) => {
                    app.saved.width_once = true;
                    true
                }
                (_, "recent" | "rc_open_file") => {
                    recent.extend(path(&v));
                    true
                }
                // The old app's Appearance checkboxes, saved under [GUI].
                ("", "axes_labels") | ("GUI", "axes_label") => {
                    put(&mut app.settings.appearance.axes_labels, &v)
                }
                ("", "peak_labels") | ("GUI", "peak_label") => {
                    put(&mut app.settings.appearance.peak_labels, &v)
                }
                ("GUI", "size") => old_size(&v)
                    .map(|size| ctx.send_viewport_cmd(egui::ViewportCommand::InnerSize(size)))
                    .is_some(),
                // The old app's other settings stay for later versions; its window state blobs do not.
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
    app.settings.restore(s.sanitized(&mut bad));
    bad
}

/// #311: a saved width 1 (3.3's only value, the app's old default) becomes Auto once. The marker keeps
/// a width 1 chosen afterwards. Only start-up calls this; settings files from Load… keep their value.
fn width_once(app: &mut App) {
    let done = std::mem::replace(&mut app.saved.width_once, true);
    if !done && app.settings.applied().an_deconvolution_width == Some(1.0) {
        let mut s = app.settings.applied().clone();
        s.an_deconvolution_width = None;
        app.settings.restore(s);
    }
}

/// The lines saved under [`KEY`].
pub fn text(app: &App) -> String {
    let mut lines: Vec<(String, String)> = (app.settings.applied().lines().into_iter())
        .map(|(k, v)| (k.to_owned(), v))
        .collect();
    let mut panels = app.panels;
    lines.extend((panels.saved()).map(|(k, on)| (format!("panel.{k}"), on.to_string())));
    let a = app.settings.appearance;
    lines.extend([
        ("axes_labels".to_owned(), a.axes_labels.to_string()),
        ("peak_labels".to_owned(), a.peak_labels.to_string()),
    ]);
    // Not answered yet: the next start asks again.
    if !app.guided.ask {
        lines.push(("mode".into(), app.guided.mode.key().into()));
    }
    // Only a choice is saved, so the default keeps following the system language.
    let lang = crate::lang::get();
    let chosen = app.saved.language || lang != crate::lang::Lang::system();
    if chosen {
        lines.push(("language".into(), lang.key().into()));
    }
    if app.saved.width_once {
        lines.push((WIDTH_ONCE.into(), "true".into()));
    }
    let library = (app.library.as_ref()).map(|l| l.path.clone());
    let nist = (app.nist.lib.as_ref()).map(|l| l.path.clone());
    let paths = [
        ("library", library.or_else(|| app.saved.library.clone())),
        ("nist_library", nist.or_else(|| app.saved.nist.clone())),
    ];
    let paths = paths.into_iter().filter_map(|(k, p)| Some((k, p?)));
    let recent = app.saved.recent.iter().map(|p| ("recent", p.clone()));
    // shortcut: a path that is not UTF-8 is not saved; store it as bytes if that ever matters.
    // Absolute, so a path given on the command line still points to the same file from elsewhere.
    lines.extend((paths.chain(recent)).filter_map(|(k, p)| {
        let p = std::path::absolute(&p).unwrap_or(p);
        Some((k.to_owned(), p.to_str()?.to_owned()))
    }));
    // An answered question replaces a kept unknown `mode`, a chosen language a kept unknown one.
    let answered = !app.guided.ask;
    lines.extend(
        (app.saved.extra.iter())
            .filter(|(k, _)| !(answered && k == "mode") && !(chosen && k == "language"))
            .cloned(),
    );
    // A line break would split a value into a second, wrong line.
    (lines.into_iter())
        .filter(|(_, v)| !v.contains(['\n', '\r']))
        .map(|(k, v)| format!("{k}={v}\n"))
        .collect()
}

/// What `app.ron` holds as eframe wrote it, or None if it is missing or does not parse.
fn read_state(file: &Path) -> Option<HashMap<String, String>> {
    parse(&std::fs::read(file).ok()?)
}

fn parse(bytes: &[u8]) -> Option<HashMap<String, String>> {
    ron::de::from_str(std::str::from_utf8(bytes).ok()?).ok()
}

/// Writes a temp file next to `file`, then renames it over: a crash leaves the old file whole.
fn write_atomic(file: &Path, bytes: &[u8]) -> std::io::Result<()> {
    let tmp = file.with_extension("ron.tmp");
    let done = std::fs::File::create(&tmp)
        .and_then(|mut f| f.write_all(bytes).and_then(|()| f.sync_all()))
        .and_then(|()| std::fs::rename(&tmp, file));
    if done.is_err() {
        std::fs::remove_file(&tmp).ok();
    }
    done
}

/// What [`guard`] tells the user.
type Note = Result<Msg, Msg>;

/// Before eframe reads `app.ron`. A good file is copied to `app.ron.bak`. A damaged one (eframe
/// drops it and would write over it) is replaced by that copy; without one it becomes the backup.
/// If neither works, eframe gets another file name, so nothing is overwritten. Returns the file
/// eframe may use and what to tell the user. A damaged file counts only if it is still damaged
/// after a short pause, so a second start does not take another window's save for damage (#305).
pub fn guard(file: Option<PathBuf>) -> (Option<PathBuf>, Option<Note>) {
    guard_after(file, || {
        std::thread::sleep(std::time::Duration::from_millis(200))
    })
}

/// [`guard`] with the pause before a damaged file is read a second time.
fn guard_after(file: Option<PathBuf>, pause: impl FnOnce()) -> (Option<PathBuf>, Option<Note>) {
    let Some(file) = file else {
        return (None, None);
    };
    if !file.exists() {
        return (Some(file), None);
    }
    let bak = file.with_extension("ron.bak");
    let lost = msg!(
        "The saved settings could not be read.",
        "Die gespeicherten Einstellungen konnten nicht gelesen werden.",
    );
    // eframe gets another name, so nothing is written over `app.ron`.
    let aside = |e: std::io::Error| {
        let other = file.with_extension("ron.new");
        let msg = msg!(
            "{lost} It could not be kept aside ({e}), so this session saves to {} instead.",
            "{lost} Sie konnten nicht beiseitegelegt werden ({e}), daher speichert diese Sitzung stattdessen in {}.",
            other.display()
        );
        (Some(other), Some(Err(msg)))
    };
    // A file we cannot read (locked by a scanner) is not a damaged one.
    let read = || std::fs::read(&file);
    let mut bytes = match read() {
        Ok(bytes) => bytes,
        Err(e) => return aside(e),
    };
    if parse(&bytes).is_none() {
        // Another window may be saving right now (eframe truncates first): look once more.
        pause();
        bytes = match read() {
            Ok(bytes) => bytes,
            // Gone: another starting window kept it aside and says so; nothing left to guard.
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return (Some(file), None),
            Err(e) => return aside(e),
        };
    }
    if parse(&bytes).is_some() {
        // A failed copy loses nothing: the file itself is fine.
        if std::fs::read(&bak).ok().as_deref() != Some(&bytes[..]) {
            write_atomic(&bak, &bytes).ok();
        }
        return (Some(file), None);
    }
    let kept = match std::fs::read(&bak).ok().filter(|b| parse(b).is_some()) {
        Some(good) => write_atomic(&file, &good).map(|()| {
            msg!(
                "The last good copy is used.",
                "Die letzte gute Kopie wird verwendet."
            )
        }),
        // Gone already: the other window got there first.
        None => std::fs::rename(&file, &bak)
            .or_else(|e| {
                (e.kind() == std::io::ErrorKind::NotFound)
                    .then_some(())
                    .ok_or(e)
            })
            .map(|()| {
                msg!(
                    "The defaults are used; the old file is kept as {} until the next start.",
                    "Die Standardwerte werden verwendet; die alte Datei bleibt bis zum nächsten Start als {} erhalten.",
                    bak.display()
                )
            }),
    };
    match kept {
        Ok(what) => (
            Some(file.clone()),
            Some(Err(msg!("{lost} {what}", "{lost} {what}"))),
        ),
        Err(e) => aside(e),
    }
}

/// The values `apply` reset, in the current language (called inside `msg!`, so in both).
fn shown(bad: &[String]) -> String {
    let bad: Vec<String> = bad.iter().map(|b| crate::errors::text(b)).collect();
    bad.join(", ")
}

/// At start: the saved state, else on the first start the old app's settings; then the libraries.
/// `guarded` is what [`guard`] said about the file.
pub fn restore(
    app: &mut App,
    ctx: &egui::Context,
    storage: Option<&dyn eframe::Storage>,
    file: Option<PathBuf>,
    guarded: Option<Note>,
) {
    app.saved.file = file;
    let saved = storage.and_then(|s| s.get_string(KEY));
    // Asked with Guided preselected unless a saved `mode` answers it.
    (app.guided.mode, app.guided.ask) = (crate::guided::Mode::Guided, true);
    let note = match saved {
        Some(text) => {
            let bad = apply(app, ctx, own_entries(&text));
            (!bad.is_empty()).then(|| {
                Err(msg!(
                    "Saved values could not be read and use their defaults: {}.",
                    "Gespeicherte Werte konnten nicht gelesen werden und nutzen ihre Standardwerte: {}.",
                    shown(&bad)
                ))
            })
        }
        // A damaged file says why on its own; the old settings are no stand-in for it.
        None if guarded.is_some() => None,
        None => old_entries().map(|(source, entries)| {
            let bad = apply(app, ctx, entries);
            let done = msg!(
                "Imported the settings of MetaboliteDetector 3.3 from {source}.",
                "Die Einstellungen von MetaboliteDetector 3.3 aus {source} übernommen."
            );
            if bad.is_empty() {
                Ok(done)
            } else {
                Err(msg!(
                    "{done} These keep their defaults: {}.",
                    "{done} Diese behalten ihre Standardwerte: {}.",
                    shown(&bad)
                ))
            }
        }),
    };
    width_once(app);
    if let Some(note) = guarded.or(note) {
        app.settings.show_note(note);
    }
    // A library that is missing now (a drive not mounted) is kept for the next start without a
    // message; one that fails to load says so.
    if let Some(p) = app.saved.library.clone().filter(|p| p.exists()) {
        app.open_library(ctx, p);
    }
    if let Some(p) = app.saved.nist.clone().filter(|p| p.exists()) {
        app.open_nist(ctx, p);
    }
}

impl App {
    /// File → Open recent: the last opened runs, newest first.
    pub fn recent_menu(&mut self, ui: &mut egui::Ui) {
        let empty = self.saved.recent.is_empty();
        ui.add_enabled_ui(!empty, |ui| {
            ui.menu_button(t("Open recent", "Zuletzt geöffnet"), |ui| {
                for p in self.saved.recent.clone() {
                    let name = p.file_name().unwrap_or(p.as_os_str()).to_string_lossy();
                    if ui
                        .button(name)
                        .on_hover_text(p.display().to_string())
                        .clicked()
                    {
                        // Loading it again puts it back on top; a file that is gone leaves the list.
                        self.saved.recent.retain(|q| q != &p);
                        self.open(ui.ctx(), p);
                    }
                }
            })
            .response
            .on_hover_text(t(
                "Open one of the last runs again",
                "Einen der letzten Läufe wieder öffnen",
            ))
            .on_disabled_hover_text(t("No runs opened yet", "Noch keine Läufe geöffnet"));
        });
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use md_core::peaks::Detector;
    use md_core::settings::GcmsSettings;

    fn apply_text(text: &str, entries: fn(&str) -> Vec<Entry>) -> (App, Vec<String>) {
        let mut app = App::default();
        let bad = apply(&mut app, &egui::Context::default(), entries(text));
        (app, bad)
    }

    /// Hand-written lines in the form the old 3.3 writes them (Qt INI quoting and escapes).
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
baseline_correction=true
rc_lib_file=
rc_open_file=/data/run 1.bin

[GUI]
axes_label=false
mw_state="@ByteArray(\0\0\0\xff\0\0\0\0\xfd)"
peak_label=true
size=@Size(1591 800)
"#;

    #[test]
    fn old_conf_values_win_over_defaults() {
        let (app, bad) = apply_text(OLD, ini_entries);
        assert!(bad.is_empty(), "{bad:?}");
        let s = &app.settings.applied();
        assert_eq!(
            s.ls_ri_diff, 10.0,
            "a saved LS_RI_DIFF beats the default 20"
        );
        assert_eq!(s.ls_resolution, 0.3);
        assert_eq!(s.an_peak_threshold_end, -5.0);
        assert_eq!(s.ls_cutoff, 0.85);
        assert!(s.baseline_correction, "the old QSettings key, default off");
        assert_eq!(
            s.an_deconvolution_width,
            Some(1.0),
            "a saved width beats the automatic one (#279)"
        );
        let shown = app.settings.appearance;
        assert_eq!((shown.axes_labels, shown.peak_labels), (false, true));
        let run = std::path::absolute("/data/run 1.bin").unwrap();
        assert_eq!(app.saved.recent, [run]);
        assert_eq!(app.saved.library, None);
        // The GUI blob is not kept.
        assert!(app.saved.extra.is_empty(), "{:?}", app.saved.extra);
    }

    /// #311: a saved width 1 becomes Auto once, whether from `app.ron` or the 3.3 import; a width
    /// chosen afterwards stays.
    #[test]
    fn saved_width_one_becomes_auto_once() {
        let width = |a: &App| a.settings.applied().an_deconvolution_width;
        let width_one = format!("{}=1\n", keys::AN_DECONVOLUTION_WIDTH);
        for (text, entries) in [
            (width_one.as_str(), own_entries as fn(&str) -> Vec<Entry>),
            (OLD, ini_entries),
        ] {
            let (mut app, _) = apply_text(text, entries);
            assert_eq!(width(&app), Some(1.0), "precondition: the file holds 1");
            width_once(&mut app);
            assert_eq!(width(&app), None);
            let saved = super::text(&app);
            assert!(saved.contains("\nauto_width_once=true\n"), "{saved}");
            // The user picks 1 again; the next start keeps it.
            let (mut next, _) = apply_text(&saved.replace("=auto\n", "=1\n"), own_entries);
            assert_eq!(width(&next), Some(1.0));
            width_once(&mut next);
            assert_eq!(width(&next), Some(1.0), "the marker stops a second switch");
        }
        let (mut app, _) = apply_text(
            &format!("{}=2.5\n", keys::AN_DECONVOLUTION_WIDTH),
            own_entries,
        );
        width_once(&mut app);
        assert_eq!(width(&app), Some(2.5), "only 1 is switched");
    }

    /// The start itself runs the switch, not only the helper.
    #[test]
    fn restore_switches_the_saved_width() {
        struct Fake(String);
        impl eframe::Storage for Fake {
            fn get_string(&self, _: &str) -> Option<String> {
                Some(self.0.clone())
            }
            fn set_string(&mut self, _: &str, _: String) {}
            fn remove_string(&mut self, _: &str) {}
            fn flush(&mut self) {}
        }
        let mut app = App::default();
        let ctx = egui::Context::default();
        let fake = Fake(format!("{}=1\n", keys::AN_DECONVOLUTION_WIDTH));
        restore(&mut app, &ctx, Some(&fake), None, None);
        assert_eq!(app.settings.applied().an_deconvolution_width, None);
        assert!(text(&app).contains("\nauto_width_once=true\n"));
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
        assert_eq!(app.settings.applied().ls_ri_diff, 10.0);
        assert_eq!(app.settings.applied().ls_resolution, 0.3);
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
        assert_eq!(app.settings.applied().ls_ri_diff, 10.0);
        assert_eq!(app.settings.applied().quant_ion_number, 4);
        assert_eq!(app.settings.applied().ls_resolution, 0.25);
        assert_eq!(app.saved.extra, [("rc_save_file".into(), String::new())]);
    }

    #[test]
    fn own_lines_round_trip_and_keep_unknown_keys() {
        let mut app = App::default();
        app.settings.edit_applied(|s| {
            (s.ls_ri_diff, s.ls_use_ri, s.ls_cutoff) = (12.5, false, 0.7);
            s.ls_mass_filter = vec![(73.0, 75.0), (147.0, 148.5)];
            s.peak_detector = Detector::Matched {
                width: 7,
                threshold: -3,
            };
        });
        app.settings.appearance.peak_labels = false;
        app.panels.info = false;
        // Absolute on every platform, with a '=' the line format must keep.
        app.saved.remember(std::env::temp_dir().join("a=b.cdf"));
        app.saved.extra.push(("future_key".into(), "x".into()));
        app.guided.mode = crate::guided::Mode::Guided;
        crate::lang::set(crate::lang::Lang::De);
        app.saved.language = true;
        let text = text(&app);
        crate::lang::set(crate::lang::Lang::En);
        assert!(text.contains("\nAN_DECONVOLUTION_WIDTH=auto\n"), "{text}");
        let (mut back, bad) = apply_text(&text, own_entries);
        assert!(bad.is_empty(), "{bad:?}");
        assert_eq!(crate::lang::get(), crate::lang::Lang::De);
        crate::lang::set(crate::lang::Lang::En);
        assert_eq!(back.guided.mode, app.guided.mode);
        // An unanswered first-start question is not saved, so the next start asks again.
        back.guided.ask = true;
        assert!(!super::text(&back).contains("mode="));
        assert_eq!(back.settings.applied(), app.settings.applied());
        let set = apply_text(
            &format!("{}=2.5\n", keys::AN_DECONVOLUTION_WIDTH),
            own_entries,
        )
        .0;
        assert_eq!(set.settings.applied().an_deconvolution_width, Some(2.5));
        assert_eq!(back.settings.appearance, app.settings.appearance);
        assert_eq!(back.panels, app.panels);
        assert_eq!(back.saved.recent, app.saved.recent);
        assert_eq!(back.saved.extra, app.saved.extra);
    }

    #[test]
    fn the_system_language_is_not_saved_and_an_unknown_one_kept() {
        let system = crate::lang::Lang::system();
        crate::lang::set(system);
        assert!(!text(&App::default()).contains("language="));
        let (app, bad) = apply_text("language=fr\n", own_entries);
        assert!(bad.is_empty(), "{bad:?}");
        assert!(text(&app).contains("\nlanguage=fr\n"));
        crate::lang::set(crate::lang::Lang::En);
    }

    #[test]
    fn unknown_mode_is_asked_again_and_kept_until_answered() {
        let mut app = App::default();
        // As the real start: asked unless a saved `mode` answers it.
        app.guided.ask = true;
        let bad = apply(
            &mut app,
            &egui::Context::default(),
            own_entries("mode=wizard\n"),
        );
        assert!(bad.is_empty(), "a kept mode is not unreadable: {bad:?}");
        assert!(
            app.guided.ask,
            "an unknown mode does not answer the question"
        );
        assert!(super::text(&app).contains("mode=wizard\n"));
        (app.guided.mode, app.guided.ask) = (crate::guided::Mode::Expert, false);
        let text = super::text(&app);
        assert!(
            text.contains("mode=expert\n") && !text.contains("wizard"),
            "{text}"
        );
    }

    #[test]
    fn bad_values_keep_defaults() {
        let (app, bad) = apply_text(
            &format!(
                "{}=abc\n{}=0\n{}=0.7\npanel.runs=yes\n{}=matched 0 5\n",
                keys::LS_RI_DIFF,
                keys::AN_NUMBER_BINS,
                keys::LS_PURE_FACTOR,
                keys::PEAK_DETECTOR
            ),
            own_entries,
        );
        assert_eq!(
            bad,
            [
                keys::LS_RI_DIFF,
                "panel.runs",
                &format!("{}: must be between 1 and 100", keys::AN_NUMBER_BINS),
                "Peak width: must be between 1 and 1000"
            ]
        );
        // Only the bad fields fall back; the good one stays.
        let expected = GcmsSettings {
            ls_pure_factor: 0.7,
            ..GcmsSettings::default()
        };
        assert_eq!(*app.settings.applied(), expected);
        assert!(app.panels.runs);
    }

    /// A fresh directory with an `app.ron` holding `recent` lines for these paths.
    fn state_dir(name: &str, recent: &[&Path]) -> (PathBuf, PathBuf) {
        let dir = std::env::temp_dir().join(format!("md-persist-{name}-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        let file = dir.join("app.ron");
        write_state(&file, recent);
        (dir, file)
    }

    fn write_state(file: &Path, recent: &[&Path]) {
        let lines: String = (recent.iter())
            .map(|p| format!("recent={}\n", p.display()))
            .collect();
        let state = HashMap::from([(KEY.to_owned(), lines)]);
        std::fs::write(file, ron::to_string(&state).unwrap()).unwrap();
    }

    fn files(dir: &Path) -> Vec<String> {
        let mut names: Vec<_> = (std::fs::read_dir(dir).unwrap())
            .filter_map(|e| e.ok()?.file_name().into_string().ok())
            .collect();
        names.sort();
        names
    }

    /// A truncated write (crash, full disk) does not cost the last good state; no backup piles up.
    #[test]
    fn truncated_file_is_restored_from_the_last_good_copy() {
        let (dir, file) = state_dir("truncated", &[Path::new("/a")]);
        for _ in 0..3 {
            assert!(guard(Some(file.clone())).1.is_none());
        }
        let good = std::fs::read(&file).unwrap();
        std::fs::write(&file, &good[..good.len() / 2]).unwrap();
        let (used, note) = guard(Some(file.clone()));
        let (names, restored) = (files(&dir), std::fs::read(&file).unwrap());
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(used, Some(file));
        assert!(note.is_some_and(|n| n.is_err()));
        assert_eq!(restored, good);
        assert_eq!(names, ["app.ron", "app.ron.bak"]);
    }

    /// A file caught mid-save by another window is whole again on the second read: left alone.
    #[test]
    fn file_being_saved_is_read_again() {
        let (dir, file) = state_dir("saving", &[Path::new("/a")]);
        let good = std::fs::read(&file).unwrap();
        std::fs::write(file.with_extension("ron.bak"), "{}").unwrap();
        std::fs::write(&file, &good[..good.len() / 2]).unwrap();
        let (used, note) =
            guard_after(Some(file.clone()), || std::fs::write(&file, &good).unwrap());
        let (now, bak) = (
            std::fs::read(&file).unwrap(),
            std::fs::read(file.with_extension("ron.bak")),
        );
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!((used, note), (Some(file), None));
        assert_eq!(now, good, "the older copy did not replace the saved file");
        assert_eq!(
            bak.unwrap(),
            good,
            "the whole file is the new last good copy"
        );
        // Kept aside by another starting window during the pause: not a reason to save elsewhere.
        let (dir, file) = state_dir("gone", &[]);
        std::fs::write(&file, "not ron {").unwrap();
        let gone = guard_after(Some(file.clone()), || std::fs::remove_file(&file).unwrap());
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(gone, (Some(file), None));
    }

    /// Without a good copy the damaged file becomes the one backup; the defaults apply.
    #[test]
    fn damaged_file_without_copy_is_kept_as_backup() {
        let (dir, file) = state_dir("damaged", &[]);
        std::fs::write(&file, "not ron {").unwrap();
        let (used, note) = guard(Some(file.clone()));
        let again = guard(Some(file.clone()));
        let (names, bak) = (files(&dir), std::fs::read(file.with_extension("ron.bak")));
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(used, Some(file.clone()));
        assert!(note.is_some_and(|n| n.is_err()));
        assert_eq!(again, (Some(file), None), "nothing left to guard");
        assert_eq!(names, ["app.ron.bak"]);
        assert_eq!(bak.unwrap(), b"not ron {");
    }

    /// If the damaged file cannot be kept aside, eframe gets another name and never overwrites it.
    #[test]
    fn failed_backup_protects_the_file() {
        let (dir, file) = state_dir("failed", &[]);
        std::fs::write(&file, "not ron {").unwrap();
        // A directory in the way makes the rename fail on every platform.
        std::fs::create_dir(file.with_extension("ron.bak")).unwrap();
        let (used, note) = guard(Some(file.clone()));
        let kept = std::fs::read(&file);
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(used, Some(file.with_extension("ron.new")));
        assert!(note.is_some_and(|n| n.is_err()));
        assert_eq!(kept.unwrap(), b"not ron {");
    }

    /// A file that cannot be read is not called damaged: the good copy must not replace it.
    #[test]
    fn unreadable_file_is_left_alone() {
        let (dir, file) = state_dir("locked", &[]);
        std::fs::remove_file(&file).unwrap();
        std::fs::create_dir(&file).unwrap();
        std::fs::write(file.with_extension("ron.bak"), "{}").unwrap();
        let (used, note) = guard(Some(file.clone()));
        let still_dir = file.is_dir();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(used, Some(file.with_extension("ron.new")));
        assert!(note.is_some_and(|n| n.is_err()) && still_dir);
    }

    /// Two windows save one after the other: the settings of the last win, the recent files merge.
    #[test]
    fn second_instance_merges_recent_files() {
        let tmp = std::env::temp_dir();
        let (a, b) = (tmp.join("md-merge-a.cdf"), tmp.join("md-merge-b.cdf"));
        let gone = tmp.join("md-merge-gone.cdf");
        std::fs::write(&a, "").unwrap();
        std::fs::write(&b, "").unwrap();
        let (dir, file) = state_dir("merge", &[&a, &gone]);
        // The second window knows only b and has not seen what the first saved.
        let mut saved = Saved {
            file: Some(file),
            ..Saved::default()
        };
        saved.remember(b.clone());
        saved.merge_disk();
        std::fs::remove_dir_all(&dir).ok();
        std::fs::remove_file(&a).ok();
        std::fs::remove_file(&b).ok();
        assert_eq!(saved.recent, [b, a], "a is added; gone is dropped");
    }

    #[test]
    fn unescape_qt() {
        assert_eq!(unescape(r#""a,b" \x41\x42\0\aP"#), "a,b AB\0\x07P");
    }
}
