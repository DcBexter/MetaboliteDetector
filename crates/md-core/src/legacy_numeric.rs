//! The deliberate roundings of MetaboliteDetector 3.3, in one place. Each operation reproduces what
//! the old x86-64 Linux build computes bit for bit, so the port lands on the old peaks and spectra;
//! the evidence is in docs/reverse-engineering/algorithms.md and gsl-clean-room.md. Callers use
//! these instead of inline `trunc`/`as` casts, so a platform difference shows up here only.

/// C `(int)x` as x86-64 `cvttsd2si` computes it: toward zero, and `INT_MIN` for NaN and for values
/// outside the `int` range (Rust's `as` saturates instead). The old program stores intensities,
/// filter values and heights as `int`.
pub(crate) fn c_int(x: f64) -> i32 {
    if x > -2_147_483_649.0 && x < 2_147_483_648.0 {
        x as i32
    } else {
        i32::MIN
    }
}

/// One point of the old `SavitzkyGolayFilter<int>` (window 5, quadratic): the smoothed value of
/// the window's centre, truncated to `int` (algorithms.md §3.1).
pub(crate) fn sg_trunc(w: [f32; 5]) -> f32 {
    let [a, b, c, d, e] = w.map(f64::from);
    c_int((-3.0 * (a + e) + 12.0 * (b + d) + 17.0 * c) / 35.0) as f32
}

/// The old `PeakEssential::roundSignificantFigures(x, 8)` (0x6fa301): rounds half up to one digit
/// more, then floors to the last one, with the digit count from `(int)log10`. So it floors rather
/// than rounds, keeps 8 significant digits below 1 and 9 from 1 up, and works on decimal places for
/// 0 and negative `x` (`log10` gives no exponent there). From 1e9 up the second digit count is
/// negative, which gives NaN like the old code.
pub(crate) fn round_sig8(x: f64) -> f64 {
    // `(int)log10(v)` via the 64-bit `cvttsd2si`: its low 32 bits are 0 for -inf and NaN.
    let exponent = |v: f64| {
        let l = v.log10();
        if l.is_finite() { l as i32 } else { 0 }
    };
    // `std::pow(10.0, unsigned)`: a negative count wraps to a huge one, giving inf. shortcut: powi
    // is exact up to 10^22 like glibc pow; beyond (|x| < 1e-14) the last bit may differ.
    let pow10 = |p: i32| {
        if p < 0 { f64::INFINITY } else { 10f64.powi(p) }
    };
    let p = 8 - exponent(x) + 1;
    let y = (0.5 + pow10(p) * x).floor() / pow10(p);
    let p = 8 - exponent(y);
    (pow10(p) * y).floor() / pow10(p)
}

/// glibc `pow(x, 2)`, which the old apex fit uses for the `x²` column. shortcut: `x * x` equals it
/// while `x² < 2^53` (bins × scan < 9.4e7); beyond, pow may round half-way cases the other way.
/// Port pow if such runs appear (gsl-clean-room.md, "Known limit").
pub(crate) fn glibc_pow2(x: f64) -> f64 {
    x * x
}

/// glibc `hypot`, which the GSL fit of the old program calls: `sqrt(x² + y²)` with the correction
/// step of Borges (2019, arXiv:1904.09481, form without FMA). The platform `hypot` (e.g. MSVC)
/// differs in the last bit, so it is written out. Inputs here are far from overflow.
pub(crate) fn glibc_hypot(x: f64, y: f64) -> f64 {
    let (a, b) = (x.abs().max(y.abs()), x.abs().min(y.abs()));
    if a == 0.0 {
        return 0.0;
    }
    let h = (a * a + b * b).sqrt();
    let (t1, t2) = if h <= 2.0 * b {
        let delta = h - b;
        (a * (2.0 * delta - a), (delta - 2.0 * (a - b)) * delta)
    } else {
        let delta = h - a;
        (
            2.0 * delta * (a - 2.0 * b),
            (4.0 * delta - b) * b + delta * delta,
        )
    };
    h - (t1 + t2) / (2.0 * h)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn c_int_truncates_like_x86() {
        assert_eq!(c_int(2.9), 2);
        assert_eq!(c_int(-2.9), -2);
        // Half-bin tie of the precise apex: `(int)(t + 0.5)` rounds 4.5 up and -4.5 toward zero.
        assert_eq!(c_int(4.5 + 0.5), 5);
        assert_eq!(c_int(-4.5 + 0.5), -4);
        assert_eq!(c_int(2_147_483_647.9), i32::MAX);
        assert_eq!(c_int(-2_147_483_648.9), i32::MIN);
        // Out of range (above 2^53 too) and NaN give the x86 "integer indefinite", not saturation.
        for x in [
            2_147_483_648.0,
            -2_147_483_649.0,
            2f64.powi(53) + 2.0,
            f64::NAN,
        ] {
            assert_eq!(c_int(x), i32::MIN, "{x}");
        }
    }

    #[test]
    fn sg_trunc_truncates_toward_zero() {
        // A flat top keeps its value; 17/35 truncates to 0, -17/35 too.
        assert_eq!(sg_trunc([7.0; 5]), 7.0);
        assert_eq!(sg_trunc([0.0, 0.0, 1.0, 0.0, 0.0]), 0.0);
        assert_eq!(sg_trunc([0.0, 0.0, -1.0, 0.0, 0.0]), 0.0);
        assert_eq!(sg_trunc([0.0, 1.0, 3.0, 1.0, 0.0]), 2.0); // 75/35
    }

    #[test]
    fn round_sig8_floors_like_the_old_code() {
        assert_eq!(round_sig8(0.148429759), 0.14842975);
        assert_eq!(round_sig8(0.0123456789), 0.012345678);
        // From 1 up the old code keeps 9 digits.
        assert_eq!(round_sig8(1.234567891), 1.23456789);
        assert_eq!(round_sig8(12345.678912), 12345.6789);
        // The first step rounds half up one digit further, which can carry into the kept digit.
        assert_eq!(round_sig8(0.123456789999), 0.12345679);
        assert_eq!(round_sig8(0.0), 0.0);
        // Negative: 9 decimals rounded, then 8 floored (toward -inf).
        assert_eq!(round_sig8(-0.123456781), -0.12345679);
        // From 1e9 the second digit count is -1, which `unsigned` wraps: inf * y, floor, inf / inf.
        assert_eq!(round_sig8(987654321.4), 987654321.0);
        assert!(round_sig8(1.5e9).is_nan());
        assert!(round_sig8(f64::NAN).is_nan());
    }

    #[test]
    fn glibc_pow2_is_exact_below_2_53() {
        // The largest bins × scan below the limit: its square is still an exact integer.
        let x = 94_906_265u64;
        assert!(glibc_pow2(x as f64) < 2f64.powi(53));
        assert_eq!(glibc_pow2(x as f64) as u64, x * x);
        assert_eq!(glibc_pow2(-3.5), 12.25);
    }

    #[test]
    fn glibc_hypot_is_close_and_symmetric() {
        assert_eq!(glibc_hypot(0.0, 0.0), 0.0);
        assert_eq!(glibc_hypot(3.0, 4.0), 5.0);
        assert_eq!(glibc_hypot(-3.0, 0.0), 3.0);
        for (x, y) in [
            (1.0, 1e-9),
            (-5.0, 12.0),
            (3.3, 7.1),
            (1e100, 2e99),
            (0.1, 0.2),
        ] {
            assert_eq!(glibc_hypot(x, y), glibc_hypot(y, x));
            let h = glibc_hypot(x, y);
            assert!((h - (x * x + y * y).sqrt()).abs() <= 1e-15 * h);
        }
    }
}
