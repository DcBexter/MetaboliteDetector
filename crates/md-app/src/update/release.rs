//! Pure parts of the update: which release is newer, which asset fits, how the app was installed,
//! and the saved update settings.

use std::path::Path;

use semver::Version;
use serde_json::Value;
use sha2::{Digest, Sha256};

use crate::lang::{t, tf};

/// A published release newer than the running app.
#[derive(Clone, Debug, PartialEq)]
pub struct Release {
    pub version: Version,
    /// Release page on GitHub, for the notes and a manual download.
    pub page: String,
    /// (file name, download URL) of every asset.
    pub assets: Vec<(String, String)>,
}

/// The newest release above `current` in the GitHub releases JSON; drafts are skipped,
/// pre-releases count (every release so far is one).
pub fn newest(json: &str, current: &Version) -> Result<Option<Release>, String> {
    let list: Value = serde_json::from_str(json).map_err(|e| {
        tf!(
            "Unreadable release list: {e}",
            "Unlesbare Versionsliste: {e}"
        )
    })?;
    let list = list.as_array().ok_or(t(
        "Unexpected release list from GitHub",
        "Unerwartete Versionsliste von GitHub",
    ))?;
    Ok(list
        .iter()
        .filter(|r| r["draft"].as_bool() == Some(false))
        .filter_map(release)
        .filter(|r| r.version > *current)
        .max_by(|a, b| a.version.cmp(&b.version)))
}

/// One entry of the list; `None` for a tag that is not `vX.Y.Z[-pre]`.
fn release(r: &Value) -> Option<Release> {
    let tag = r["tag_name"].as_str()?;
    let version = Version::parse(tag.strip_prefix('v').unwrap_or(tag)).ok()?;
    let assets = (r["assets"].as_array()?.iter())
        .filter_map(|a| {
            let name = a["name"].as_str()?;
            Some((
                name.to_owned(),
                a["browser_download_url"].as_str()?.to_owned(),
            ))
        })
        .collect();
    let page = r["html_url"].as_str()?.to_owned();
    Some(Release {
        version,
        page,
        assets,
    })
}

/// The SHA-256 of `name` in a `sha256sum` listing (`<hex>  <name>`; `*` marks binary mode).
pub fn expected_sum<'a>(sums: &'a str, name: &str) -> Option<&'a str> {
    sums.lines().find_map(|line| {
        let (hex, file) = line.trim_end().split_once(' ')?;
        let file = file.strip_prefix([' ', '*']).unwrap_or(file);
        let valid = hex.len() == 64 && hex.bytes().all(|b| b.is_ascii_hexdigit());
        (valid && file == name).then_some(hex)
    })
}

/// Accepts `data` only if `SHA256SUMS.txt` of the same release lists `name` with its exact hash.
pub fn verify(data: &[u8], sums: &str, name: &str) -> Result<(), String> {
    let expected = expected_sum(sums, name).ok_or(tf!(
        "SHA256SUMS.txt does not list {name}",
        "SHA256SUMS.txt führt {name} nicht auf"
    ))?;
    let actual: String = Sha256::digest(data)
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect();
    if actual.eq_ignore_ascii_case(expected) {
        Ok(())
    } else {
        Err(tf!(
            "The checksum of {name} does not match; nothing was installed.",
            "Die Prüfsumme von {name} stimmt nicht; nichts wurde installiert."
        ))
    }
}

/// How the running app was installed, which decides the asset and the way it is replaced.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Kind {
    /// The MSI; `all_users` when an admin installed it under Program Files (`ALLUSERS=1`).
    Msi {
        all_users: bool,
    },
    Portable,
    Linux,
}

/// The one asset name accepted for `kind`, as the release workflow names it.
pub fn asset_name(kind: Kind, version: &Version) -> String {
    let suffix = match kind {
        Kind::Msi { .. } => "windows-x64.msi",
        Kind::Portable => "windows-x64-portable.exe",
        Kind::Linux => "linux-x64.tar.gz",
    };
    format!("MetaboliteDetector-{version}-{suffix}")
}

/// The install kind from the exe path: the MSI installs into `<LocalAppData>\Programs` (per user)
/// or `<ProgramFiles>` (per machine); any other place on Windows is the portable exe. `None` where
/// no asset is published (macOS).
pub fn kind(
    os: &str,
    exe: &str,
    local_app_data: Option<&str>,
    program_files: Option<&str>,
) -> Option<Kind> {
    // Strings, not `Path`, so the Windows rules also run in the Linux tests; NTFS ignores case.
    let norm = |p: &str| p.replace('\\', "/").trim_end_matches('/').to_lowercase();
    let dir = norm(exe).rsplit_once('/')?.0.to_owned();
    let installed = |base: Option<&str>, sub: &str| {
        base.is_some_and(|b| !b.is_empty() && dir == format!("{}/{sub}", norm(b)))
    };
    match os {
        "linux" => Some(Kind::Linux),
        "windows" if installed(local_app_data, "programs/metabolitedetector") => {
            Some(Kind::Msi { all_users: false })
        }
        "windows" if installed(program_files, "metabolitedetector") => {
            Some(Kind::Msi { all_users: true })
        }
        "windows" => Some(Kind::Portable),
        _ => None,
    }
}

/// Seconds between automatic checks.
const DAY: u64 = 24 * 60 * 60;

/// Update settings, kept in `update.txt` as `key=value` lines.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Config {
    /// Check automatically at start.
    pub auto: bool,
    /// Unix seconds of the last successful automatic or manual check.
    pub last_check: u64,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            auto: true,
            last_check: 0,
        }
    }
}

impl Config {
    /// Missing, unknown or garbled lines keep their defaults.
    pub fn parse(text: &str) -> Self {
        let mut c = Self::default();
        for line in text.lines() {
            match line.trim().split_once('=') {
                Some(("auto", "1")) => c.auto = true,
                Some(("auto", "0")) => c.auto = false,
                Some(("last_check", t)) => c.last_check = t.parse().unwrap_or(c.last_check),
                _ => {}
            }
        }
        c
    }

    pub fn text(&self) -> String {
        format!(
            "auto={}\nlast_check={}\n",
            u8::from(self.auto),
            self.last_check
        )
    }

    pub fn read(path: &Path) -> Self {
        std::fs::read_to_string(path).map_or_else(|_| Self::default(), |t| Self::parse(&t))
    }

    /// Best effort: without a writable config dir the app only checks again next start.
    pub fn write(&self, path: &Path) {
        if let Some(dir) = path.parent() {
            std::fs::create_dir_all(dir).ok();
        }
        std::fs::write(path, self.text()).ok();
    }

    /// Whether the automatic check runs now; a clock set back also counts as due.
    pub fn due(&self, now: u64) -> bool {
        self.auto && (now < self.last_check || now - self.last_check >= DAY)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &str) -> Version {
        Version::parse(s).unwrap()
    }

    const RELEASES: &str = r#"[
        {"tag_name": "v0.3.0", "draft": true, "prerelease": false,
         "html_url": "https://github.com/x/draft", "assets": []},
        {"tag_name": "v0.2.0-rc1", "draft": false, "prerelease": true,
         "html_url": "https://github.com/DcBexter/MetaboliteDetector/releases/tag/v0.2.0-rc1",
         "assets": [{"name": "SHA256SUMS.txt", "browser_download_url": "https://example/sums"}]},
        {"tag_name": "v0.1.2", "draft": false, "prerelease": true,
         "html_url": "https://github.com/x/0.1.2", "assets": []},
        {"tag_name": "nightly", "draft": false, "prerelease": true,
         "html_url": "https://github.com/x/nightly", "assets": []},
        {"tag_name": "v0.1.1", "draft": false, "prerelease": true,
         "html_url": "https://github.com/x/0.1.1", "assets": []}
    ]"#;

    #[test]
    fn newest_skips_drafts_and_takes_pre_releases() {
        let r = newest(RELEASES, &v("0.1.1")).unwrap().unwrap();
        assert_eq!(r.version, v("0.2.0-rc1"));
        assert!(r.page.ends_with("/tag/v0.2.0-rc1"));
        assert_eq!(
            r.assets,
            [(
                "SHA256SUMS.txt".to_owned(),
                "https://example/sums".to_owned()
            )]
        );
    }

    #[test]
    fn nothing_newer_is_none() {
        assert_eq!(newest(RELEASES, &v("0.2.0")).unwrap(), None);
        assert_eq!(newest("[]", &v("0.1.1")).unwrap(), None);
        assert!(newest("{\"message\": \"rate limit\"}", &v("0.1.1")).is_err());
        assert!(newest("<html>", &v("0.1.1")).is_err());
    }

    #[test]
    fn versions_order_like_semver() {
        assert!(v("0.1.2") > v("0.1.1"));
        assert!(v("0.2.0-rc1") > v("0.1.1"));
        assert!(v("0.2.0-rc1") < v("0.2.0"));
        // Equal or older is no update.
        let one = r#"[{"tag_name": "v0.1.1", "draft": false, "html_url": "u", "assets": []}]"#;
        assert_eq!(newest(one, &v("0.1.1")).unwrap(), None);
        assert_eq!(newest(one, &v("0.1.2")).unwrap(), None);
        assert_eq!(
            newest(one, &v("0.1.0")).unwrap().map(|r| r.version),
            Some(v("0.1.1"))
        );
    }

    const HASH_ABC: &str = "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad";

    #[test]
    fn sums_are_read_and_checked() {
        let sums = format!(
            "{HASH_ABC}  MetaboliteDetector-0.1.2-linux-x64.tar.gz\n\
             {}  MetaboliteDetector-0.1.2-windows-x64.msi\n\
             {HASH_ABC} *binary.exe\n\
             garbage line\n",
            "0".repeat(64)
        );
        let linux = "MetaboliteDetector-0.1.2-linux-x64.tar.gz";
        assert_eq!(expected_sum(&sums, linux), Some(HASH_ABC));
        assert_eq!(expected_sum(&sums, "binary.exe"), Some(HASH_ABC));
        assert_eq!(expected_sum(&sums, "MetaboliteDetector-0.1.2"), None);
        assert_eq!(verify(b"abc", &sums, linux), Ok(()));
        assert_eq!(verify(b"abc", &sums.to_uppercase(), "BINARY.EXE"), Ok(()));
        let msi = "MetaboliteDetector-0.1.2-windows-x64.msi";
        assert!(
            verify(b"abc", &sums, msi)
                .unwrap_err()
                .contains("does not match")
        );
        assert!(verify(b"abd", &sums, linux).is_err());
        assert!(
            verify(b"abc", &sums, "other.exe")
                .unwrap_err()
                .contains("does not list")
        );
    }

    #[test]
    fn asset_names_follow_the_release_workflow() {
        let ver = v("0.2.0-rc1");
        let msi = Kind::Msi { all_users: true };
        assert_eq!(
            asset_name(msi, &ver),
            "MetaboliteDetector-0.2.0-rc1-windows-x64.msi"
        );
        assert_eq!(
            asset_name(Kind::Portable, &v("0.1.2")),
            "MetaboliteDetector-0.1.2-windows-x64-portable.exe"
        );
        assert_eq!(
            asset_name(Kind::Linux, &v("0.1.2")),
            "MetaboliteDetector-0.1.2-linux-x64.tar.gz"
        );
    }

    #[test]
    fn install_kind_comes_from_the_path() {
        let local = Some(r"C:\Users\Ana\AppData\Local");
        let pf = Some(r"C:\Program Files\");
        let k = |exe| kind("windows", exe, local, pf);
        assert_eq!(
            k(r"C:\Users\Ana\AppData\Local\Programs\MetaboliteDetector\metabolite-detector.exe"),
            Some(Kind::Msi { all_users: false })
        );
        assert_eq!(
            k(r"c:\users\ana\appdata\local\PROGRAMS\metabolitedetector\metabolite-detector.exe"),
            Some(Kind::Msi { all_users: false })
        );
        assert_eq!(
            k(r"C:\Program Files\MetaboliteDetector\metabolite-detector.exe"),
            Some(Kind::Msi { all_users: true })
        );
        assert_eq!(
            k(r"C:\Users\Ana\Downloads\MetaboliteDetector-0.1.1-windows-x64-portable.exe"),
            Some(Kind::Portable)
        );
        // A subfolder of the install dir is not the MSI's.
        assert_eq!(
            k(r"C:\Program Files\MetaboliteDetector\old\metabolite-detector.exe"),
            Some(Kind::Portable)
        );
        assert_eq!(
            kind("windows", r"D:\md.exe", None, None),
            Some(Kind::Portable)
        );
        assert_eq!(
            kind("linux", "/opt/md/metabolite-detector", None, None),
            Some(Kind::Linux)
        );
        assert_eq!(kind("macos", "/Applications/md", None, None), None);
    }

    #[test]
    fn config_round_trips_and_garbage_is_default() {
        let dir = std::env::temp_dir().join(format!("md-update-{}", std::process::id()));
        let path = dir.join("sub").join("update.txt");
        assert_eq!(Config::read(&path), Config::default());
        let c = Config {
            auto: false,
            last_check: 1_760_000_000,
        };
        c.write(&path);
        assert_eq!(Config::read(&path), c);
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(
            Config::parse("\u{0}auto=maybe\nlast_check=-3\n=="),
            Config::default()
        );
        assert_eq!(Config::parse("auto=0").last_check, 0);
    }

    #[test]
    fn checks_at_most_daily_when_enabled() {
        let c = Config {
            auto: true,
            last_check: 1000,
        };
        assert!(!c.due(1000 + DAY - 1));
        assert!(c.due(1000 + DAY));
        assert!(c.due(10), "clock set back");
        assert!(!Config { auto: false, ..c }.due(1000 + 2 * DAY));
        assert!(Config::default().due(1_760_000_000));
    }
}
