//! Cancel and progress of long jobs, shared by the app and the command line.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering::Relaxed};

/// Asks a job to stop. Cooperative: the job looks at [`Cancel::is_cancelled`] between its steps and
/// returns early (usually `None`); work in a single step runs to its end. Clones share one flag, and
/// once set it stays set.
#[derive(Clone, Default, Debug)]
pub struct Cancel(
    Arc<AtomicBool>,
    /// Tests only: polls left until the flag sets itself; 0 is off.
    #[cfg(test)]
    Arc<AtomicUsize>,
);

impl Cancel {
    pub fn cancel(&self) {
        self.0.store(true, Relaxed);
    }

    pub fn is_cancelled(&self) -> bool {
        #[cfg(test)]
        if self.1.load(Relaxed) > 0 && self.1.fetch_sub(1, Relaxed) == 1 {
            self.cancel();
        }
        self.0.load(Relaxed)
    }

    /// Tests only: cancels itself at the `polls`-th call of [`Self::is_cancelled`], so a test can stop a
    /// job inside its loop instead of before it.
    #[cfg(test)]
    pub(crate) fn after(polls: usize) -> Self {
        let c = Self::default();
        c.1.store(polls, Relaxed);
        c
    }
}

/// How far a job is: `done` of `total` steps. `done` only grows. Clones share the counters.
#[derive(Clone, Default, Debug)]
pub struct Progress(Arc<(AtomicUsize, AtomicUsize)>);

impl Progress {
    pub fn set_total(&self, total: usize) {
        self.0.1.store(total, Relaxed);
    }

    /// One more step is finished.
    pub fn step(&self) {
        self.0.0.fetch_add(1, Relaxed);
    }

    pub fn done(&self) -> usize {
        self.0.0.load(Relaxed)
    }

    /// 0 while the job has not said how many steps it has.
    pub fn total(&self) -> usize {
        self.0.1.load(Relaxed)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn clones_share_state() {
        let (c, p) = (Cancel::default(), Progress::default());
        let (c2, p2) = (c.clone(), p.clone());
        assert!(!c2.is_cancelled());
        c.cancel();
        assert!(c2.is_cancelled());
        p.set_total(3);
        p.step();
        p2.step();
        assert_eq!((p2.done(), p2.total()), (2, 3));
    }
}
