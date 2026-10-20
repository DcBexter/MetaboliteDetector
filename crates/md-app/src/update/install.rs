//! Network and file work of the update: fetch the release list, download and verify an asset,
//! put the new binaries in place and start them.

use std::ffi::OsString;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::time::Duration;

use semver::Version;

use super::release::{self, Kind, Release};

const RELEASES_URL: &str =
    "https://api.github.com/repos/DcBexter/MetaboliteDetector/releases?per_page=20";
/// GitHub's API rejects requests without one.
const USER_AGENT: &str = concat!("MetaboliteDetector/", env!("CARGO_PKG_VERSION"));
/// Largest download accepted; the assets are 10–30 MB, ureq's default cap is 10 MB.
const MAX_DOWNLOAD: u64 = 512 * 1024 * 1024;

/// HTTPS only, also after GitHub's redirect of a download to its file host.
fn agent(timeout: Duration) -> ureq::Agent {
    ureq::Agent::config_builder()
        .https_only(true)
        .timeout_connect(Some(Duration::from_secs(15)))
        .timeout_global(Some(timeout))
        .build()
        .into()
}

fn get(agent: &ureq::Agent, url: &str, accept: &str) -> Result<Vec<u8>, String> {
    let mut response = agent
        .get(url)
        .header("User-Agent", USER_AGENT)
        .header("Accept", accept)
        .call()
        .map_err(|e| format!("Could not reach {url}: {e}"))?;
    response
        .body_mut()
        .with_config()
        .limit(MAX_DOWNLOAD)
        .read_to_vec()
        .map_err(|e| format!("Download from {url} failed: {e}"))
}

/// The newest published release above `current`, if any.
pub fn check(current: &Version) -> Result<Option<Release>, String> {
    let agent = agent(Duration::from_secs(30));
    let body = get(&agent, RELEASES_URL, "application/vnd.github+json")?;
    release::newest(&String::from_utf8_lossy(&body), current)
}

/// Downloads the asset for `kind`, checks it against the release's `SHA256SUMS.txt`, installs it
/// and starts the new version; the caller then closes the app. Nothing changes on any error before
/// the files are replaced.
pub fn install(release: &Release, kind: Kind, exe: &Path) -> Result<(), String> {
    let name = release::asset_name(kind, &release.version);
    let url = |file: &str| {
        (release.assets.iter())
            .find(|(n, u)| n == file && u.starts_with("https://"))
            .map(|(_, u)| u.clone())
            .ok_or(format!("The release has no {file}"))
    };
    let agent = agent(Duration::from_secs(600));
    let octets = "application/octet-stream";
    let sums = get(&agent, &url("SHA256SUMS.txt")?, octets)?;
    let data = get(&agent, &url(&name)?, octets)?;
    release::verify(&data, &String::from_utf8_lossy(&sums), &name)?;
    match kind {
        Kind::Msi { all_users } => install_msi(&data, &name, all_users, exe),
        Kind::Portable => replace_exe(&data, exe),
        Kind::Linux => {
            replace_from_tar_gz(&data, exe).map_err(|e| {
                let dir = exe.parent().unwrap_or(exe).display();
                format!(
                    "Could not replace the binaries in {dir}: {e}. Download the update from {}",
                    release.page
                )
            })?;
            start(exe)
        }
    }
}

/// `exe` with `suffix` appended to its file name (`app.exe` → `app.exe.old`).
fn side(exe: &Path, suffix: &str) -> PathBuf {
    let mut s = OsString::from(exe.as_os_str());
    s.push(suffix);
    s.into()
}

/// Removes the exe a portable update renamed; best effort, it may still be closing.
pub fn remove_old(exe: &Path) {
    std::fs::remove_file(side(exe, ".old")).ok();
}

fn start(exe: &Path) -> Result<(), String> {
    Command::new(exe)
        .spawn()
        .map(drop)
        .map_err(|e| format!("Installed, but {} did not start: {e}", exe.display()))
}

/// Runs the MSI from a fresh temp dir through a hidden helper that waits for this app to exit,
/// then starts the installed exe again.
fn install_msi(data: &[u8], name: &str, all_users: bool, exe: &Path) -> Result<(), String> {
    let nanos = (std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH))
        .map_or(0, |d| d.as_nanos());
    let dir = std::env::temp_dir().join(format!(
        "MetaboliteDetector-update-{}-{nanos}",
        std::process::id()
    ));
    let msi = dir.join(name);
    std::fs::create_dir(&dir)
        .and_then(|()| std::fs::write(&msi, data))
        .map_err(|e| format!("Could not save {}: {e}", msi.display()))?;
    // The ping is a ~2 s wait so msiexec does not meet the exe still running.
    let all = if all_users { " ALLUSERS=1" } else { "" };
    let script = format!(
        "/c ping -n 3 127.0.0.1 >nul & msiexec /i \"{}\" /passive /norestart{all} & start \"\" \"{}\"",
        msi.display(),
        exe.display()
    );
    hidden_cmd(&script)
}

#[cfg(windows)]
fn hidden_cmd(script: &str) -> Result<(), String> {
    use std::os::windows::process::CommandExt;
    const CREATE_NO_WINDOW: u32 = 0x0800_0000;
    // Raw, because cmd parses its own quotes and Rust's argument escaping would break them.
    Command::new("cmd")
        .raw_arg(script)
        .creation_flags(CREATE_NO_WINDOW)
        .spawn()
        .map(drop)
        .map_err(|e| format!("Could not start the installer: {e}"))
}

#[cfg(not(windows))]
fn hidden_cmd(_: &str) -> Result<(), String> {
    Err("The MSI installs only on Windows".into())
}

/// Swaps the portable exe: Windows lets a running exe be renamed but not overwritten.
fn replace_exe(data: &[u8], exe: &Path) -> Result<(), String> {
    let (new, old) = (side(exe, ".new"), side(exe, ".old"));
    let fail = |what: &str, e: io::Error| {
        std::fs::remove_file(&new).ok();
        format!("Could not {what} {}: {e}", exe.display())
    };
    std::fs::write(&new, data).map_err(|e| fail("write the update next to", e))?;
    std::fs::remove_file(&old).ok();
    std::fs::rename(exe, &old).map_err(|e| fail("rename", e))?;
    if let Err(e) = std::fs::rename(&new, exe) {
        std::fs::rename(&old, exe).ok();
        return Err(fail("replace", e));
    }
    start(exe)
}

/// Replaces `exe` and, if it exists, `md-cli` next to it with `metabolite-detector` and `md-cli`
/// from the release tar.gz (`MetaboliteDetector-X.Y.Z-linux-x64/<file>`). All are written as
/// `<file>.new` before any is renamed over its original, so a failed extraction changes nothing.
pub fn replace_from_tar_gz(archive: &[u8], exe: &Path) -> io::Result<()> {
    let cli = exe.with_file_name("md-cli");
    let mut targets = vec![("metabolite-detector", exe.to_path_buf(), false)];
    if cli.exists() {
        targets.push(("md-cli", cli, false));
    }
    let result =
        write_new(archive, &mut targets).and_then(|()| match targets.iter().find(|t| !t.2) {
            Some((name, ..)) => Err(io::Error::other(format!("the archive has no {name}"))),
            None => Ok(()),
        });
    if let Err(e) = result {
        for (_, path, _) in &targets {
            std::fs::remove_file(side(path, ".new")).ok();
        }
        return Err(e);
    }
    for (_, path, _) in &targets {
        std::fs::rename(side(path, ".new"), path)?;
    }
    Ok(())
}

/// Writes each wanted file of the archive to `<target>.new` and marks it found.
fn write_new(archive: &[u8], targets: &mut [(&str, PathBuf, bool)]) -> io::Result<()> {
    let mut tar = tar::Archive::new(flate2::read::GzDecoder::new(archive));
    for entry in tar.entries()? {
        let mut entry = entry?;
        let path = entry.path()?.into_owned();
        // Only regular files exactly one folder down; the names written are ours, never the archive's.
        let wanted = targets.iter_mut().find(|(name, ..)| {
            path.components().count() == 2 && path.file_name() == Some(name.as_ref())
        });
        let Some((_, target, found)) = wanted.filter(|_| entry.header().entry_type().is_file())
        else {
            continue;
        };
        let new = side(target, ".new");
        let mut out = std::fs::File::create(&new)?;
        io::copy(&mut entry, &mut out)?;
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            out.set_permissions(std::fs::Permissions::from_mode(0o755))?;
        }
        *found = true;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tar.gz like the release workflow's, with `files` under one top folder.
    fn archive(files: &[(&str, &[u8])]) -> Vec<u8> {
        let gz = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        let mut tar = tar::Builder::new(gz);
        for (name, data) in files {
            let mut header = tar::Header::new_gnu();
            header.set_size(data.len() as u64);
            header.set_mode(0o644);
            header.set_cksum();
            let path = format!("MetaboliteDetector-0.1.2-linux-x64/{name}");
            tar.append_data(&mut header, path, *data).unwrap();
        }
        tar.into_inner().unwrap().finish().unwrap()
    }

    fn folder(name: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("md-update-{name}-{}", std::process::id()));
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    #[test]
    fn linux_binaries_are_replaced_from_the_tar_gz() {
        let dir = folder("replace");
        let (exe, cli) = (dir.join("metabolite-detector"), dir.join("md-cli"));
        std::fs::write(&exe, "old app").unwrap();
        std::fs::write(&cli, "old cli").unwrap();
        let data = archive(&[
            ("README.md", b"readme"),
            ("metabolite-detector", b"new app"),
            ("md-cli", b"new cli"),
        ]);
        replace_from_tar_gz(&data, &exe).unwrap();
        assert_eq!(std::fs::read(&exe).unwrap(), b"new app");
        assert_eq!(std::fs::read(&cli).unwrap(), b"new cli");
        assert!(!dir.join("README.md").exists());
        assert!(!side(&exe, ".new").exists());
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let mode = std::fs::metadata(&exe).unwrap().permissions().mode();
            assert_eq!(mode & 0o777, 0o755);
        }
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn md_cli_is_only_replaced_where_it_was() {
        let dir = folder("no-cli");
        let exe = dir.join("metabolite-detector");
        std::fs::write(&exe, "old app").unwrap();
        let data = archive(&[("metabolite-detector", b"new app"), ("md-cli", b"cli")]);
        replace_from_tar_gz(&data, &exe).unwrap();
        assert_eq!(std::fs::read(&exe).unwrap(), b"new app");
        assert!(!dir.join("md-cli").exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn incomplete_archive_changes_nothing() {
        let dir = folder("incomplete");
        let (exe, cli) = (dir.join("metabolite-detector"), dir.join("md-cli"));
        std::fs::write(&exe, "old app").unwrap();
        std::fs::write(&cli, "old cli").unwrap();
        // The app is there but md-cli is not; and a garbled archive.
        let data = archive(&[("metabolite-detector", b"new app")]);
        assert!(replace_from_tar_gz(&data, &exe).is_err());
        assert!(replace_from_tar_gz(b"not a tar.gz", &exe).is_err());
        assert_eq!(std::fs::read(&exe).unwrap(), b"old app");
        assert_eq!(std::fs::read(&cli).unwrap(), b"old cli");
        assert!(!side(&exe, ".new").exists());
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn deeper_or_misplaced_entries_are_ignored() {
        let dir = folder("nested");
        let exe = dir.join("metabolite-detector");
        std::fs::write(&exe, "old app").unwrap();
        let data = archive(&[("sub/metabolite-detector", b"evil")]);
        assert!(replace_from_tar_gz(&data, &exe).is_err());
        assert_eq!(std::fs::read(&exe).unwrap(), b"old app");
        std::fs::remove_dir_all(&dir).ok();
    }
}
