//! Opening files and keeping the per-run background results (loading, SIC) up to date.

use std::path::PathBuf;
use std::sync::{Arc, mpsc};
use std::time::Instant;

use eframe::egui;

use crate::{App, Run, poll_job, settings, spawn};

impl App {
    /// Every file loads on its own thread, so many runs load in parallel.
    pub fn open(&mut self, ctx: &egui::Context, path: PathBuf) {
        let p = path.clone();
        let key = self.settings.key();
        let rx = spawn(ctx, move || {
            md_core::read_run(&p)
                .and_then(|data| match data.scan_count() {
                    0 => Err("no scans".into()),
                    _ => {
                        let components = settings::analyze(&data, &key)?;
                        Ok(Run {
                            name: p.file_stem().unwrap_or_default().to_string_lossy().into(),
                            minutes: data.scan_times_s.iter().map(|t| t / 60.0).collect(),
                            tic: data.tic(),
                            data: Arc::new(data),
                            path: p.clone(),
                            sic: Vec::new(),
                            sic_for: None,
                            sic_job: None,
                            components,
                            analysis: settings::Analysis::new(key),
                            scan: 0,
                            component: None,
                            visible: true,
                            ident: None,
                            ident_job: None,
                        })
                    }
                })
                .map_err(|e| format!("{}: {e}", p.display()))
        });
        self.load_started.get_or_insert_with(Instant::now);
        self.loading.push((path, rx));
    }

    pub fn poll_loading(&mut self) {
        let (count, active) = (
            self.runs.len(),
            self.runs.get(self.active).map(|r| r.path.clone()),
        );
        self.loading.retain(|(_, rx)| {
            match rx.try_recv() {
                Ok(Ok(run)) => self.runs.push(run),
                Ok(Err(e)) => self.error = Some(e),
                Err(mpsc::TryRecvError::Empty) => return true,
                Err(mpsc::TryRecvError::Disconnected) => self.error = Some("loader crashed".into()),
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
            self.load_note = format!(
                "{} runs open, last batch loaded in {:.2?}",
                self.runs.len(),
                t.elapsed()
            );
        }
    }

    /// Collects finished SICs and starts one job per run whose SIC is stale.
    pub fn update_sics(&mut self, ctx: &egui::Context) {
        let key = (self.mz, self.tolerance);
        for run in &mut self.runs {
            if let Some((k, rx)) = &run.sic_job
                && let Some(result) = poll_job(rx, &format!("{}: SIC", run.name))
            {
                // A crashed SIC is not retried for this key; the trace stays empty.
                let sic = result.unwrap_or_else(|e| {
                    self.error = Some(e);
                    Vec::new()
                });
                (run.sic, run.sic_for) = (sic, Some(*k));
                run.sic_job = None;
            }
            if run.sic_for != Some(key) && run.sic_job.as_ref().map(|j| j.0) != Some(key) {
                let data = run.data.clone();
                run.sic_job = Some((key, spawn(ctx, move || data.sic(key.0, key.1))));
            }
        }
    }
}
