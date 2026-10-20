//! File → "Export all runs as NetCDF…": writes every open run into one folder in the background.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, mpsc};
use std::time::Duration;

use eframe::egui;
use md_core::Chromatogram;
use md_core::centroid::{CentroidSettings, ExportKind, reduce};

use crate::layout::ToolWindow;
use crate::{App, poll_job, spawn};

struct Job {
    total: usize,
    done: Arc<AtomicUsize>,
    cancel: Arc<AtomicBool>,
    rx: mpsc::Receiver<Vec<String>>,
}

#[derive(Default)]
pub struct ExportAll {
    pub open: bool,
    kind: ExportKind,
    /// Only the runs drawn in the plots.
    only_shown: bool,
    dir: Option<PathBuf>,
    job: Option<Job>,
    /// One line per failed file, or the summary of the last export.
    note: String,
    errors: Vec<String>,
}

/// Writes each run as `<name><suffix>.cdf` in `dir`, the name the single-run export proposes. A
/// repeated run name gets a counter, so no file overwrites another. Returns one line per failure.
fn export_runs(
    runs: &[(String, Arc<Chromatogram>)],
    dir: &Path,
    kind: ExportKind,
    s: &CentroidSettings,
    done: &AtomicUsize,
    cancel: &AtomicBool,
) -> Vec<String> {
    let (mut errors, mut used) = (Vec::new(), HashSet::new());
    for (name, data) in runs {
        if cancel.load(Relaxed) {
            break;
        }
        let mut file = format!("{name}{}.cdf", kind.suffix());
        for n in 2.. {
            if used.insert(file.clone()) {
                break;
            }
            file = format!("{name}{}_{n}.cdf", kind.suffix());
        }
        let path = dir.join(&file);
        let result = if kind == ExportKind::Profile {
            md_core::andi::write(&path, data)
        } else {
            md_core::andi::write_with_tic(&path, &reduce(data, kind, s), &data.tic())
        };
        if let Err(e) = result {
            errors.push(format!("{file}: {e}"));
        }
        done.fetch_add(1, Relaxed);
    }
    errors
}

impl App {
    pub fn export_all_window(&mut self, ctx: &egui::Context) {
        let e = &mut self.export_all;
        if let Some(job) = &e.job {
            match poll_job(&job.rx, "NetCDF export") {
                Some(r) => {
                    let (done, total) = (job.done.load(Relaxed), job.total);
                    e.errors = r.unwrap_or_else(|m| vec![m]);
                    e.note = match job.cancel.load(Relaxed) {
                        true => format!("Cancelled after {done} of {total} runs."),
                        false => format!("Exported {} of {total} runs.", total - e.errors.len()),
                    };
                    e.job = None;
                }
                None => ctx.request_repaint_after(Duration::from_millis(200)),
            }
        }
        let mut open = e.open;
        egui::Window::new("Export all runs as NetCDF")
            .tool_window(ctx, 13, [420.0, 320.0])
            .open(&mut open)
            .show(ctx, |ui| self.export_all_contents(ui));
        self.export_all.open = open;
    }

    fn export_all_contents(&mut self, ui: &mut egui::Ui) {
        let (runs, applied) = (&self.runs, &self.settings.applied);
        let e = &mut self.export_all;
        ui.add_enabled_ui(e.job.is_none(), |ui| {
            for kind in ExportKind::ALL {
                ui.radio_value(&mut e.kind, kind, kind.label());
            }
            ui.checkbox(&mut e.only_shown, "Only runs shown in the plots");
            ui.horizontal(|ui| {
                if ui.button("Folder…").clicked()
                    && let Some(dir) = rfd::FileDialog::new().pick_folder()
                {
                    e.dir = Some(dir);
                }
                ui.label(
                    e.dir
                        .as_ref()
                        .map_or("No folder chosen".into(), |d| d.display().to_string()),
                );
            });
        });
        match &e.job {
            Some(job) => {
                let d = job.done.load(Relaxed);
                ui.horizontal(|ui| {
                    ui.add(
                        egui::ProgressBar::new(d as f32 / job.total.max(1) as f32)
                            .text(format!("{d} / {} runs", job.total))
                            .desired_width(240.0),
                    );
                    if ui.button("Cancel").clicked() {
                        job.cancel.store(true, Relaxed);
                    }
                });
            }
            None => {
                let picked: Vec<_> = runs
                    .iter()
                    .filter(|r| r.visible || !e.only_shown)
                    .map(|r| (r.name.clone(), r.data.clone()))
                    .collect();
                let dir = e.dir.clone().filter(|_| !picked.is_empty());
                if ui
                    .add_enabled(
                        dir.is_some(),
                        egui::Button::new(format!("Export {} runs", picked.len())),
                    )
                    .clicked()
                    && let Some(dir) = dir
                {
                    let (done, cancel) = (
                        Arc::new(AtomicUsize::new(0)),
                        Arc::new(AtomicBool::new(false)),
                    );
                    let (d, c, kind, s) = (
                        done.clone(),
                        cancel.clone(),
                        e.kind,
                        CentroidSettings::from(applied),
                    );
                    let total = picked.len();
                    let rx = spawn(ui.ctx(), move || {
                        export_runs(&picked, &dir, kind, &s, &d, &c)
                    });
                    (e.note, e.errors) = (String::new(), Vec::new());
                    e.job = Some(Job {
                        total,
                        done,
                        cancel,
                        rx,
                    });
                }
            }
        }
        if !e.note.is_empty() {
            ui.label(&e.note);
        }
        for line in &e.errors {
            ui.colored_label(ui.visuals().error_fg_color, line);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_each_run_and_reports_failures_per_file() {
        let run = |mz: Vec<f32>, intensity: Vec<f32>, offsets| {
            Arc::new(Chromatogram {
                scan_times_s: vec![1.5, 2.0],
                offsets,
                mz,
                intensity,
                exact_masses: false,
            })
        };
        let c = run(
            vec![73.0, 147.5, 73.0],
            vec![10.0, 0.5, 20.0],
            vec![0, 2, 3],
        );
        // No points at all: the writer refuses it.
        let empty = run(vec![], vec![], vec![0, 0, 0]);
        let dir = std::env::temp_dir().join(format!("md-export-all-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let runs = vec![
            ("a".to_string(), c.clone()),
            ("a".to_string(), c),
            ("bad".to_string(), empty),
        ];
        let done = AtomicUsize::new(0);
        let errors = export_runs(
            &runs,
            &dir,
            ExportKind::Profile,
            &CentroidSettings::from(&md_core::settings::GcmsSettings::default()),
            &done,
            &AtomicBool::new(false),
        );
        let names: Vec<_> = (std::fs::read_dir(&dir).unwrap())
            .map(|f| f.unwrap().file_name().into_string().unwrap())
            .collect();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(done.load(Relaxed), 3);
        assert_eq!(errors.len(), 1);
        assert!(errors[0].starts_with("bad_profile.cdf"));
        assert!(
            names.contains(&"a_profile.cdf".into()) && names.contains(&"a_profile_2.cdf".into())
        );
    }
}
