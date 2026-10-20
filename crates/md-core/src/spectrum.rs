//! Spectrum display normalisation, as in the old `SpectrumPlot::getNormalizedData`
//! (`docs/reverse-engineering/algorithms.md` §7c): every intensity is divided by one norm of the spectrum.

/// Old `Normalization` enum (`NONE`, `BP`, `SSUM`, `SUM`) with the menu names of the issue.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Normalization {
    /// Raw intensities (norm 1).
    #[default]
    None,
    /// Divide by the base peak (largest intensity); the tallest ion becomes 1.
    BasePeak,
    /// Divide by the vector length `sqrt(Σ i²)` (old `SSUM`, `getSpectrumVectorLength`).
    Unit,
    /// Divide by the total signal `Σ i` (old `SUM`, `getTotalSignal`); fractions add up to 1.
    Relative,
}

impl Normalization {
    pub const ALL: [Self; 4] = [Self::None, Self::BasePeak, Self::Unit, Self::Relative];

    pub fn label(self) -> &'static str {
        match self {
            Self::None => "None",
            Self::BasePeak => "Base peak",
            Self::Unit => "Unit",
            Self::Relative => "Relative",
        }
    }
}

/// Intensities divided by the norm of `mode`. A spectrum whose norm is not a positive finite number
/// (empty, all zero) stays unchanged: the old code divided by zero and plotted NaN.
pub fn normalize(intensity: &[f64], mode: Normalization) -> Vec<f64> {
    let values = intensity.iter().copied();
    let norm = match mode {
        Normalization::None => 1.0,
        Normalization::BasePeak => values.clone().fold(0.0, f64::max),
        Normalization::Unit => values.clone().map(|v| v * v).sum::<f64>().sqrt(),
        Normalization::Relative => values.clone().sum(),
    };
    let norm = if norm.is_finite() && norm > 0.0 {
        norm
    } else {
        1.0
    };
    values.map(|v| v / norm).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SPECTRUM: [f64; 3] = [3.0, 4.0, 0.0];

    #[test]
    fn each_mode_divides_by_its_norm() {
        let n = |m| normalize(&SPECTRUM, m);
        assert_eq!(n(Normalization::None), [3.0, 4.0, 0.0]);
        assert_eq!(n(Normalization::BasePeak), [0.75, 1.0, 0.0]);
        assert_eq!(n(Normalization::Unit), [0.6, 0.8, 0.0]);
        let rel = n(Normalization::Relative);
        assert!((rel[0] - 3.0 / 7.0).abs() < 1e-12 && (rel[1] - 4.0 / 7.0).abs() < 1e-12);
        assert!((rel.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn degenerate_spectra_stay_unchanged() {
        for mode in Normalization::ALL {
            assert!(normalize(&[], mode).is_empty());
            assert_eq!(normalize(&[0.0, 0.0], mode), [0.0, 0.0]);
        }
    }
}
