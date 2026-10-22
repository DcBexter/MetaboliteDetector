//! File → "Export all runs as NetCDF…": writes every open run into one folder in the background.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Duration;

use eframe::egui;
use md_core::Chromatogram;
use md_core::andi::AndiError;
use md_core::centroid::{CentroidSettings, ExportKind, reduce};
use md_core::job::{Cancel, Progress};

use crate::lang::{Msg, msg, t, tf};
use crate::layout::ToolWindow;
use crate::{App, Job};

#[derive(Default)]
pub struct ExportAll {
    pub open: bool,
    kind: ExportKind,
    /// Only the runs drawn in the plots.
    only_shown: bool,
    dir: Option<PathBuf>,
    job: Option<Job<Vec<Msg>>>,
    /// One line per failed file, or the summary of the last export.
    note: Msg,
    errors: Vec<Msg>,
}

/// Writes each run as `<name><suffix>.cdf` in `dir`, the name the single-run export proposes. A
/// repeated run name gets a counter, so no file overwrites another. Returns one line per failure.
fn export_runs(
    runs: &[(String, Arc<Chromatogram>)],
    dir: &Path,
    kind: ExportKind,
    s: &CentroidSettings,
    done: &Progress,
    cancel: &Cancel,
) -> Vec<Msg> {
    let (mut errors, mut used) = (Vec::new(), HashSet::new());
    // Quitting mid-write leaves `<file>.cdf.tmp` behind (the writer renames it when done): clear them and say so.
    for f in std::fs::read_dir(dir).into_iter().flatten().flatten() {
        let n = f.file_name().to_string_lossy().into_owned();
        if n.to_lowercase().ends_with(".cdf.tmp") && std::fs::remove_file(f.path()).is_ok() {
            errors.push(msg!(
                "{n}: leftover of an interrupted export, removed",
                "{n}: Rest eines abgebrochenen Exports, entfernt"
            ));
        }
    }
    for (name, data) in runs {
        if cancel.is_cancelled() {
            break;
        }
        let mut file = format!("{name}{}.cdf", kind.suffix());
        for n in 2.. {
            // Lower case: NTFS and macOS treat names that differ only in case as one file.
            if used.insert(file.to_lowercase()) {
                break;
            }
            file = format!("{name}{}_{n}.cdf", kind.suffix());
        }
        let path = dir.join(&file);
        // The single-run dialog asks before replacing; here a file in the folder (maybe a source run) is kept.
        if path.exists() {
            errors.push(msg!(
                "{file}: already exists, skipped",
                "{file}: existiert bereits, übersprungen"
            ));
            done.step();
            continue;
        }
        let result = if kind == ExportKind::Profile {
            md_core::andi::write_with_tic(&path, data, &data.tic(), cancel)
        } else {
            md_core::andi::write_with_tic(&path, &reduce(data, kind, s), &data.tic(), cancel)
        };
        if let Err(AndiError::Cancelled) = result {
            break;
        }
        if let Err(e) = result {
            errors.push(format!("{file}: {e}").into());
        }
        done.step();
    }
    errors
}

impl App {
    pub fn export_all_window(&mut self, ctx: &egui::Context) {
        let e = &mut self.export_all;
        if let Some(job) = &e.job {
            match job.poll(|| t("NetCDF export", "Der NetCDF-Export").into()) {
                Some(r) => {
                    let (done, total) = (job.progress.done(), job.progress.total());
                    e.errors = r.unwrap_or_else(|m| vec![m]);
                    let failed = e.errors.len();
                    e.note = if job.cancel.is_cancelled() && done < total {
                        msg!(
                            "Cancelled: {done} of {total} runs done, {failed} failed.",
                            "Abgebrochen: {done} von {total} Läufen fertig, {failed} fehlgeschlagen."
                        )
                    } else {
                        msg!(
                            "{done} of {total} runs done, {failed} failed.",
                            "{done} von {total} Läufen fertig, {failed} fehlgeschlagen."
                        )
                    };
                    e.job = None;
                }
                None => ctx.request_repaint_after(Duration::from_millis(200)),
            }
        }
        let mut open = e.open;
        egui::Window::new(t(
            "Export all runs as NetCDF",
            "Alle Läufe als NetCDF exportieren",
        ))
        .id(egui::Id::unique("Export all runs as NetCDF"))
        .tool_window(ctx, 13, [420.0, 320.0])
        .open(&mut open)
        .show(ctx, |ui| self.export_all_contents(ui));
        self.export_all.open = open;
    }

    fn export_all_contents(&mut self, ui: &mut egui::Ui) {
        let (runs, applied) = (&self.runs, self.settings.applied());
        let e = &mut self.export_all;
        ui.add_enabled_ui(e.job.is_none(), |ui| {
            for kind in ExportKind::ALL {
                let label = match kind {
                    ExportKind::Profile => t("Profile", "Profil"),
                    ExportKind::Centroid => t("Centroid", "Zentroid"),
                    ExportKind::Nominal => t("Nominal mass", "Nominalmasse"),
                };
                ui.radio_value(&mut e.kind, kind, label);
            }
            ui.checkbox(
                &mut e.only_shown,
                t(
                    "Only runs shown in the plots",
                    "Nur die in den Grafiken gezeigten Läufe",
                ),
            );
            ui.horizontal(|ui| {
                if ui.button(t("Folder…", "Ordner…")).clicked()
                    && let Some(dir) = rfd::FileDialog::new().pick_folder()
                {
                    e.dir = Some(dir);
                }
                ui.label(
                    e.dir
                        .as_ref()
                        .map_or(t("No folder chosen", "Kein Ordner gewählt").into(), |d| {
                            d.display().to_string()
                        }),
                );
            });
        });
        match &e.job {
            Some(job) => {
                let (d, total) = (job.progress.done(), job.progress.total());
                ui.horizontal(|ui| {
                    ui.add(
                        egui::ProgressBar::new(d as f32 / total.max(1) as f32)
                            .text(tf!("{d} / {} runs", "{d} / {} Läufe", total))
                            .desired_width(240.0),
                    );
                    if ui.button(t("Cancel", "Abbrechen")).clicked() {
                        job.cancel.cancel();
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
                        egui::Button::new(tf!(
                            "Export {} runs",
                            "{} Läufe exportieren",
                            picked.len()
                        )),
                    )
                    .clicked()
                    && let Some(dir) = dir
                {
                    let (kind, s) = (e.kind, CentroidSettings::from(applied));
                    (e.note, e.errors) = (Msg::default(), Vec::new());
                    let total = picked.len();
                    let job = Job::spawn(ui.ctx(), move |done, cancel| {
                        export_runs(&picked, &dir, kind, &s, done, cancel)
                    });
                    job.progress.set_total(total);
                    e.job = Some(job);
                }
            }
        }
        if !e.note.is_empty() {
            ui.label(&*e.note);
        }
        egui::ScrollArea::vertical().show(ui, |ui| {
            for line in &e.errors {
                ui.colored_label(ui.visuals().error_fg_color, &**line);
            }
        });
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
                baseline: Default::default(),
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
        std::fs::write(dir.join("keep_profile.cdf"), b"old").unwrap();
        std::fs::write(dir.join("orphan.cdf.tmp"), b"half").unwrap();
        let runs = vec![
            ("a".to_string(), c.clone()),
            ("a".to_string(), c.clone()),
            ("bad".to_string(), empty),
            ("keep".to_string(), c.clone()),
            ("A".to_string(), c),
        ];
        let done = Progress::default();
        let errors = export_runs(
            &runs,
            &dir,
            ExportKind::Profile,
            &CentroidSettings::from(&md_core::settings::GcmsSettings::default()),
            &done,
            &Cancel::default(),
        );
        let names: Vec<_> = (std::fs::read_dir(&dir).unwrap())
            .map(|f| f.unwrap().file_name().into_string().unwrap())
            .collect();
        let kept = std::fs::read(dir.join("keep_profile.cdf")).unwrap();
        std::fs::remove_dir_all(&dir).ok();
        assert_eq!(kept, b"old");
        assert_eq!(done.done(), 5);
        assert_eq!(errors.len(), 3);
        assert!(errors[0].starts_with("orphan.cdf.tmp") && errors[0].contains("removed"));
        assert!(errors[1].starts_with("bad_profile.cdf"));
        assert!(errors[2].contains("already exists"));
        // Three runs whose names differ at most in case give three files, plus the kept one.
        assert_eq!(names.len(), 4, "{names:?}");
    }
}
