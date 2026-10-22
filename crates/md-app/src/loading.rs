//! Opening files and keeping the per-run background results (loading, SIC) up to date.

use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use eframe::egui;
use md_core::job::Cancel;

use crate::lang::{Msg, msg, t};
use crate::settings::{self, Key};
use crate::{App, Job, Loaded, Run};

/// Reads and analyses one run on a job thread. The read itself cannot be interrupted; a set `cancel`
/// stops before the analysis and inside it.
fn load_run(p: &Path, key: Key, correct: bool, cancel: &Cancel) -> Result<Loaded, Msg> {
    md_core::read_run(p)
        .map(|mut data| {
            data.set_baseline_correction(correct);
            data
        })
        .and_then(|data| match data.scan_count() {
            0 => Err("no scans".into()),
            // The read itself cannot be interrupted; the analysis after it is not started.
            _ if cancel.is_cancelled() => Err("cancelled".into()),
            _ => {
                // Results the old program saved next to the run replace the analysis
                // until the settings change. A broken side file only gives a warning.
                let mut warnings = Vec::new();
                let saved =
                    md_core::cmp::saved_components(p, &data.scan_times_s).unwrap_or_else(|e| {
                        warnings.push(e.to_string().into());
                        None
                    });
                let (components, peak_masses, peaks) = match saved {
                    Some(saved) => saved,
                    None => settings::analyze(&data, &key, cancel)?,
                };
                let saved_calibration = md_core::ri::saved_calibration(p).unwrap_or_else(|e| {
                    warnings.push(msg!("calibration: {e}", "Kalibrierung: {e}"));
                    None
                });
                let run = Run {
                    name: p.file_stem().unwrap_or_default().to_string_lossy().into(),
                    minutes: data.scan_times_s.iter().map(|t| t / 60.0).collect(),
                    tic: data.tic(),
                    data: Arc::new(data),
                    path: p.to_path_buf(),
                    sic: Vec::new(),
                    sic_for: None,
                    sic_job: None,
                    d_tic: Default::default(),
                    d_sic: Default::default(),
                    components,
                    saved_calibration,
                    peak_masses,
                    peaks,
                    analysis: settings::Analysis::new(key),
                    scan: 0,
                    component: None,
                    visible: true,
                    ident: None,
                    ident_job: None,
                    rev: 0,
                };
                Ok((run, warnings))
            }
        })
        .map_err(|e| Msg::from(format!("{}: {e}", p.display())))
}

impl App {
    /// Every file loads on its own thread, so many runs load in parallel.
    pub fn open(&mut self, ctx: &egui::Context, path: PathBuf) {
        let p = path.clone();
        let key = self.settings.key();
        let correct = self.settings.applied().baseline_correction;
        let job = Job::spawn(ctx, move |_, cancel| load_run(&p, key, correct, cancel));
        self.load_started.get_or_insert_with(Instant::now);
        self.loading.push((path, job));
    }

    pub fn poll_loading(&mut self) {
        let (count, active) = (
            self.runs.len(),
            self.runs.get(self.active).map(|r| r.path.clone()),
        );
        self.loading.retain(|(path, job)| {
            match job.poll(|| format!("{}: {}", path.display(), t("opening", "Öffnen"))) {
                Some(Ok(Ok((mut run, warnings)))) => {
                    if let Some(cal) = run.saved_calibration.take() {
                        self.ri.set_own(run.path.clone(), cal);
                    }
                    if !warnings.is_empty() {
                        self.error = Some(Msg::new(|| {
                            let warnings: Vec<&str> = warnings.iter().map(|w| &**w).collect();
                            format!("{}: {}", run.path.display(), warnings.join("; "))
                        }));
                    }
                    self.saved.remember(run.path.clone());
                    self.runs.push(run);
                }
                Some(Ok(Err(e)) | Err(e)) => self.error = Some(e),
                None => return true,
            }
            false
        });
        // Runs finish in any order; sorting by path keeps list order and colours stable.
        if self.runs.len() != count {
            self.runs.sort_by(|a, b| a.path.cmp(&b.path));
            self.active = active
                .and_then(|p| self.runs.iter().position(|r| r.path == p))
                .unwrap_or(0);
        }
        if self.loading.is_empty()
            && let Some(t) = self.load_started.take()
        {
            let took = crate::view::took(t.elapsed());
            self.load_note = match self.runs.len() {
                1 => msg!(
                    "1 run open, loaded in {took}",
                    "1 Lauf geöffnet, geladen in {took}"
                ),
                n => msg!(
                    "{n} runs open, last batch loaded in {took}",
                    "{n} Läufe geöffnet, letzter Stapel geladen in {took}"
                ),
            };
        }
    }

    /// Collects finished SICs and starts one job per run whose SIC is stale.
    pub fn update_sics(&mut self, ctx: &egui::Context) {
        let key = (self.mz, self.tolerance);
        for run in &mut self.runs {
            if let Some((k, job)) = &run.sic_job
                && let Some(result) = job.poll(|| format!("{}: SIC", run.name))
            {
                // A crashed SIC is not retried for this key; the trace stays empty.
                let sic = result.unwrap_or_else(|e| {
                    self.error = Some(e);
                    Vec::new()
                });
                (run.sic, run.sic_for) = (sic, Some(*k));
                run.d_sic.take();
                run.sic_job = None;
            }
            if run.sic_for != Some(key) && run.sic_job.as_ref().map(|j| j.0) != Some(key) {
                let data = run.data.clone();
                run.sic_job = Some((
                    key,
                    Job::spawn(ctx, move |_, c| {
                        data.sic(key.0, key.1, c).unwrap_or_default()
                    }),
                ));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cancelled_open_stops_before_the_analysis() {
        let run = md_core::Chromatogram {
            offsets: vec![0, 1, 2],
            mz: vec![57.0, 57.0],
            intensity: vec![1.0, 2.0],
            scan_times_s: vec![0.0, 1.0],
            exact_masses: false,
            baseline: Default::default(),
        };
        let path = std::env::temp_dir().join(format!("md-open-cancel-{}.cdf", std::process::id()));
        md_core::andi::write(&path, &run).unwrap();
        let key = crate::settings::Settings::default().key();
        let cancel = Cancel::default();
        cancel.cancel();
        let result = load_run(&path, key, false, &cancel);
        std::fs::remove_file(&path).ok();
        assert!(result.err().unwrap().ends_with(": cancelled"));
    }
}
