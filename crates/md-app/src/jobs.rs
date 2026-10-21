//! Background jobs: run work on its own thread and poll the result from the UI thread.

use std::sync::mpsc;

use eframe::egui;
use md_core::job::{Cancel, Progress};

use crate::lang::{Msg, tf};

/// A background job: `work` runs on its own thread and the UI polls the result, so the window never
/// freezes. The job reports its [`Progress`] and stops when [`Cancel`] is set; dropping the `Job`
/// (closing the window, starting the next job in its place) cancels it, and its result is discarded.
pub struct Job<T> {
    rx: mpsc::Receiver<T>,
    pub progress: Progress,
    pub cancel: Cancel,
}

impl<T: Send + 'static> Job<T> {
    /// Wakes the UI when `work` is done. `work` gets the job's progress and cancel handles.
    pub fn spawn(
        ctx: &egui::Context,
        work: impl FnOnce(&Progress, &Cancel) -> T + Send + 'static,
    ) -> Self {
        let (tx, rx) = mpsc::channel();
        let (progress, cancel) = (Progress::default(), Cancel::default());
        let (p, c) = (progress.clone(), cancel.clone());
        let ctx = ctx.clone();
        let lang = crate::lang::get();
        std::thread::spawn(move || {
            // Texts made by the job come in the UI's language.
            crate::lang::set(lang);
            // Dropping `tx` before the repaint lets the UI see a panic as `Disconnected` right away.
            if let Ok(v) = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| work(&p, &c))) {
                tx.send(v).ok();
            }
            drop(tx);
            ctx.request_repaint();
        });
        Self {
            rx,
            progress,
            cancel,
        }
    }
}

impl<T> Job<T> {
    /// `None` while the job runs, `Err` when its thread panicked (the sender is dropped unsent), so
    /// callers can report it and not wait or restart forever.
    pub fn poll(&self, what: impl Fn() -> String) -> Option<Result<T, Msg>> {
        match self.rx.try_recv() {
            Ok(v) => Some(Ok(v)),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err(Msg::new(|| {
                let what = what();
                tf!("{what} crashed", "{what} ist abgestürzt")
            }))),
        }
    }
}

impl<T> Drop for Job<T> {
    fn drop(&mut self) {
        self.cancel.cancel();
    }
}

#[cfg(test)]
mod tests {
    use std::time::{Duration, Instant};

    use super::*;

    fn wait<T: Send + 'static>(job: &Job<T>) -> Result<T, Msg> {
        let started = Instant::now();
        loop {
            if let Some(r) = job.poll(|| "test job".into()) {
                return r;
            }
            assert!(started.elapsed().as_secs() < 10, "job never finished");
            std::thread::yield_now();
        }
    }

    fn spin(c: &Cancel) {
        while !c.is_cancelled() {
            std::thread::yield_now();
        }
    }

    #[test]
    fn panicking_job_is_reported_not_lost() {
        let job = Job::spawn(&egui::Context::default(), |_, _| -> u8 { panic!("boom") });
        assert_eq!(wait(&job).unwrap_err().to_string(), "test job crashed");
    }

    #[test]
    fn job_reports_progress() {
        let job = Job::spawn(&egui::Context::default(), |p, _| {
            p.set_total(2);
            p.step();
            p.step();
        });
        wait(&job).unwrap();
        assert_eq!((job.progress.done(), job.progress.total()), (2, 2));
    }

    #[test]
    fn cancel_and_drop_stop_the_job() {
        let ctx = egui::Context::default();
        let job = Job::spawn(&ctx, |_, c| {
            spin(c);
            true
        });
        job.cancel.cancel();
        assert!(wait(&job).unwrap());
        // Dropping the job sets the flag the worker still holds.
        let (tx, rx) = mpsc::channel();
        drop(Job::spawn(&ctx, move |_, c| {
            spin(c);
            tx.send(()).ok();
        }));
        rx.recv_timeout(Duration::from_secs(10)).unwrap();
    }
}
