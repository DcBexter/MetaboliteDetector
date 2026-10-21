//! Background jobs: run work on its own thread and poll the result from the UI thread.

use std::sync::mpsc;

use eframe::egui;

use crate::lang::{Msg, tf};

/// Runs `job` on its own thread and wakes the UI when it is done, so the window never freezes.
pub fn spawn<T: Send + 'static>(
    ctx: &egui::Context,
    job: impl FnOnce() -> T + Send + 'static,
) -> mpsc::Receiver<T> {
    let (tx, rx) = mpsc::channel();
    let ctx = ctx.clone();
    let lang = crate::lang::get();
    std::thread::spawn(move || {
        // Texts made by the job come in the UI's language.
        crate::lang::set(lang);
        // Dropping `tx` before the repaint lets the UI see a panic as `Disconnected` right away.
        if let Ok(v) = std::panic::catch_unwind(std::panic::AssertUnwindSafe(job)) {
            tx.send(v).ok();
        }
        drop(tx);
        ctx.request_repaint();
    });
    rx
}

/// Polls a job from `spawn`: `None` while it runs, `Err` when its thread panicked (the sender is
/// dropped unsent), so callers can report it and not wait or restart forever.
pub fn poll_job<T>(rx: &mpsc::Receiver<T>, what: impl Fn() -> String) -> Option<Result<T, Msg>> {
    match rx.try_recv() {
        Ok(v) => Some(Ok(v)),
        Err(mpsc::TryRecvError::Empty) => None,
        Err(mpsc::TryRecvError::Disconnected) => Some(Err(Msg::new(|| {
            let what = what();
            tf!("{what} crashed", "{what} ist abgestürzt")
        }))),
    }
}

#[cfg(test)]
mod tests {
    use std::time::Instant;

    use super::*;

    #[test]
    fn panicking_job_is_reported_not_lost() {
        let ctx = egui::Context::default();
        let rx = spawn(&ctx, || -> u8 { panic!("boom") });
        let started = Instant::now();
        let result = loop {
            if let Some(r) = poll_job(&rx, || "test job".into()) {
                break r;
            }
            assert!(started.elapsed().as_secs() < 10, "job never finished");
            std::thread::yield_now();
        };
        assert_eq!(result.unwrap_err().to_string(), "test job crashed");
    }
}
