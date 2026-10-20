//! Help → Check for updates: finds a newer GitHub release and installs it after a click (#235).

mod install;
mod release;

use std::path::{Path, PathBuf};
use std::sync::mpsc::Receiver;

use eframe::egui;
use semver::Version;

use crate::jobs::{poll_job, spawn};
use crate::layout::ToolWindow;
use release::{Config, Kind, Release};

const RELEASES_PAGE: &str = "https://github.com/DcBexter/MetaboliteDetector/releases";

enum State {
    Idle,
    /// `manual` checks report their errors; the automatic one stays silent.
    Checking {
        rx: Receiver<Result<Option<Release>, String>>,
        manual: bool,
    },
    UpToDate,
    Available(Release),
    Installing(Version, Receiver<Result<(), String>>),
    Failed(String),
}

pub struct Updater {
    pub open: bool,
    config: Config,
    /// `update.txt` in the user config dir; `None` until `start` (tests) or without one.
    path: Option<PathBuf>,
    state: State,
}

impl Default for Updater {
    fn default() -> Self {
        Self {
            open: false,
            config: Config::default(),
            path: None,
            state: State::Idle,
        }
    }
}

fn current() -> Version {
    Version::parse(env!("CARGO_PKG_VERSION")).expect("Cargo versions are semver")
}

fn now() -> u64 {
    (std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)).map_or(0, |d| d.as_secs())
}

/// `%APPDATA%\MetaboliteDetector`, else `$XDG_CONFIG_HOME` or `~/.config` + `/metabolite-detector`.
fn config_path() -> Option<PathBuf> {
    let env = |k| {
        std::env::var_os(k)
            .map(PathBuf::from)
            .filter(|p| p.is_absolute())
    };
    let dir = if cfg!(windows) {
        env("APPDATA")?.join("MetaboliteDetector")
    } else {
        (env("XDG_CONFIG_HOME").or_else(|| Some(env("HOME")?.join(".config"))))?
            .join("metabolite-detector")
    };
    Some(dir.join("update.txt"))
}

/// Only x64 builds are published.
fn install_kind(exe: &Path) -> Option<Kind> {
    let env = |k| std::env::var(k).ok();
    let kind = release::kind(
        std::env::consts::OS,
        &exe.to_string_lossy(),
        env("LOCALAPPDATA").as_deref(),
        env("ProgramFiles").as_deref(),
    );
    kind.filter(|_| cfg!(target_arch = "x86_64"))
}

impl Updater {
    /// At start: loads the settings, removes an exe a portable update replaced, runs the daily check.
    pub fn start(&mut self, ctx: &egui::Context) {
        self.path = config_path();
        if let Some(p) = &self.path {
            self.config = Config::read(p);
        }
        if let Ok(exe) = std::env::current_exe() {
            install::remove_old(&exe);
        }
        // Debug builds never check on their own, so `cargo run` stays offline.
        if !cfg!(debug_assertions) && self.config.due(now()) {
            self.check(ctx, false);
        }
    }

    fn check(&mut self, ctx: &egui::Context, manual: bool) {
        match &mut self.state {
            // A click during the automatic check still gets its error shown.
            State::Checking { manual: m, .. } => *m |= manual,
            State::Installing(..) => {}
            _ => {
                let rx = spawn(ctx, || install::check(&current()));
                self.state = State::Checking { rx, manual };
            }
        }
    }

    fn save(&self) {
        if let Some(p) = &self.path {
            self.config.write(p);
        }
    }

    /// Saves only the check time, re-read from the file so another window's switch survives.
    fn save_last_check(&mut self) {
        // An unreadable file keeps the in-memory switch rather than the default.
        if let Some(Ok(text)) = self.path.as_ref().map(std::fs::read_to_string) {
            self.config = Config::parse(&text);
        }
        self.config.last_check = now();
        self.save();
    }

    pub fn menu(&mut self, ui: &mut egui::Ui) {
        if ui
            .button("Check for updates…")
            .on_hover_text("Look for a newer release on GitHub")
            .clicked()
        {
            self.open = true;
            self.check(ui.ctx(), true);
        }
        if ui
            .checkbox(&mut self.config.auto, "Check for updates at start")
            .on_hover_text("At most once a day; nothing is installed without your click")
            .changed()
        {
            self.save();
        }
    }

    fn poll(&mut self, ctx: &egui::Context) {
        match &self.state {
            State::Checking { rx, manual } => {
                let manual = *manual;
                let Some(result) = poll_job(rx, "The update check") else {
                    return;
                };
                self.state = match result.and_then(|r| r) {
                    Ok(found) => {
                        self.save_last_check();
                        // A found update opens the window: the notice, which installs nothing by itself.
                        self.open |= found.is_some();
                        found.map_or(State::UpToDate, State::Available)
                    }
                    Err(e) if manual => State::Failed(e),
                    Err(_) => State::Idle,
                };
            }
            State::Installing(_, rx) => match poll_job(rx, "The update").map(|r| r.and_then(|r| r))
            {
                Some(Ok(())) => ctx.send_viewport_cmd(egui::ViewportCommand::Close),
                Some(Err(e)) => self.state = State::Failed(e),
                None => {}
            },
            _ => {}
        }
    }

    pub fn window(&mut self, ctx: &egui::Context) {
        self.poll(ctx);
        let mut open = self.open;
        egui::Window::new("Updates")
            .tool_window(ctx, 12, [340.0, 240.0])
            .open(&mut open)
            .show(ctx, |ui| self.contents(ui));
        // "Later" closes it from inside.
        self.open &= open;
    }

    fn contents(&mut self, ui: &mut egui::Ui) {
        ui.label(format!("Installed version: {}", env!("CARGO_PKG_VERSION")));
        let mut install = None;
        match &self.state {
            State::Idle => {}
            State::Checking { .. } => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label("Looking for a newer release…");
                });
            }
            State::UpToDate => {
                ui.label("This is the newest version.");
            }
            State::Available(r) => {
                ui.strong(format!("Version {} is available.", r.version));
                ui.hyperlink_to("Release notes and downloads", &r.page);
                ui.label(
                    "Installing downloads the update, verifies its checksum and restarts \
                     MetaboliteDetector. Unsaved results are lost.",
                );
                let exe = std::env::current_exe().ok();
                let kind = exe.as_deref().and_then(install_kind);
                // A debug build would overwrite cargo's own binaries.
                let blocked = if cfg!(debug_assertions) {
                    Some("This is a debug build: install the release by hand.")
                } else if kind.is_none() {
                    Some(
                        "No automatic install on this platform: download it from the release page.",
                    )
                } else {
                    None
                };
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(blocked.is_none(), egui::Button::new("Download and install"))
                        .clicked()
                        && let (Some(exe), Some(kind)) = (exe, kind)
                    {
                        install = Some((r.clone(), kind, exe));
                    }
                    if ui.button("Later").clicked() {
                        self.open = false;
                    }
                });
                if let Some(why) = blocked {
                    ui.label(why);
                }
            }
            State::Installing(version, _) => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Downloading and installing {version}…"));
                });
            }
            State::Failed(e) => {
                ui.colored_label(ui.visuals().error_fg_color, e);
                ui.hyperlink_to("Releases page", RELEASES_PAGE);
            }
        }
        if matches!(self.state, State::Idle | State::UpToDate | State::Failed(_))
            && ui.button("Check again").clicked()
        {
            self.check(ui.ctx(), true);
        }
        if let Some((release, kind, exe)) = install {
            let version = release.version.clone();
            let rx = spawn(ui.ctx(), move || install::install(&release, kind, &exe));
            self.state = State::Installing(version, rx);
        }
    }
}
