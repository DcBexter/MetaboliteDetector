//! Least-squares parabola through three points, bit for bit like the old
//! `tools::Statistics::polynomialRegression(x, y, 2)` (§3.6). The old program solves the
//! Vandermonde system `[1, x, x²]` with GSL's SVD least squares; the steps below were rebuilt from
//! black-box outputs only (no GSL source), see
//! `docs/reverse-engineering/gsl-clean-room.md`. Every operation and its order matters: the
//! rounding of the result decides ties of the precise apex.

type M3 = [[f64; 3]; 3];

const I3: M3 = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

/// Coefficients `c` of `y = c0 + c1 x + c2 x²` through the three points.
pub(crate) fn quadratic(x: [f64; 3], y: [f64; 3]) -> [f64; 3] {
    // shortcut: the old code takes glibc pow(x, 2), equal to x * x while x² < 2^53 (bins × scan
    // < 9.4e7); beyond, pow may round half-way cases the other way. Port pow if such runs appear.
    let mut a = x.map(|x| [1.0, x, x * x]);
    // Column balance by powers of two.
    let mut d = [0.0; 3];
    for j in 0..3 {
        let f = power_of_two(a.iter().fold(0.0, |s, r| s + r[j].abs()));
        for r in &mut a {
            r[j] /= f;
        }
        d[j] = f;
    }
    let (u, s, v) = svd_mod(a);
    let mut t = [0.0; 3];
    for (j, t) in t.iter_mut().enumerate() {
        *t = (0..3).fold(0.0, |acc, i| acc + y[i] * u[i][j]);
    }
    let mut c = [0.0; 3];
    for i in 0..3 {
        let acc = (0..3).fold(0.0, |acc, j| {
            let q = if s[j] > f64::EPSILON * s[0] {
                v[i][j] * (1.0 / s[j])
            } else {
                0.0
            };
            acc + q * t[j]
        });
        c[i] = acc / d[i];
    }
    c
}

/// Power of two `f` with `s / f` in `[0.5, 1]`; 1 for 0.
fn power_of_two(mut s: f64) -> f64 {
    let mut f = 1.0;
    if s == 0.0 {
        return f;
    }
    if s > 1.0 {
        while s > 1.0 {
            f *= 2.0;
            s /= 2.0;
        }
    } else {
        while s < 0.5 {
            f /= 2.0;
            s *= 2.0;
        }
    }
    f
}

/// `sqrt(x² + y²)` with the correction step of Borges (2019, arXiv:1904.09481, form without FMA);
/// equals the glibc `hypot` the old program calls. The platform `hypot` (e.g. MSVC) differs in the
/// last bit, so it is written out. Inputs here are far from overflow.
fn hypot(x: f64, y: f64) -> f64 {
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

/// Euclidean norm with running scale (reference BLAS `dnrm2`).
fn nrm2(v: &[f64]) -> f64 {
    let (mut scale, mut ssq) = (0.0f64, 1.0f64);
    for x in v.iter().filter(|x| **x != 0.0) {
        let ax = x.abs();
        if scale < ax {
            ssq = 1.0 + ssq * (scale / ax) * (scale / ax);
            scale = ax;
        } else {
            ssq += (ax / scale) * (ax / scale);
        }
    }
    scale * ssq.sqrt()
}

/// Householder reflector `I - tau w wᵀ` with `w = [1, v1..]` that maps `v` to `[beta, 0..]`:
/// returns `tau`, leaves `beta` in `v[0]` and `w1..` in `v[1..]`.
fn householder(v: &mut [f64]) -> f64 {
    if v.len() == 1 {
        return 0.0;
    }
    let alpha = v[0];
    let xnorm = nrm2(&v[1..]);
    if xnorm == 0.0 {
        return 0.0;
    }
    let beta = if alpha >= 0.0 {
        -hypot(alpha, xnorm)
    } else {
        hypot(alpha, xnorm)
    };
    let tau = (beta - alpha) / beta;
    let r = 1.0 / (alpha - beta);
    for x in &mut v[1..] {
        *x *= r;
    }
    v[0] = beta;
    tau
}

/// `A[r0.., c0..] = (I - tau w wᵀ) A[r0.., c0..]`, `w = [1, v[1..]]`.
#[allow(
    clippy::needless_range_loop,
    reason = "column-wise like the old order of operations"
)]
fn reflect_rows(tau: f64, v: &[f64], a: &mut M3, r0: usize, c0: usize) {
    if tau == 0.0 {
        return;
    }
    for j in c0..3 {
        let w = (r0 + 1..3).fold(a[r0][j], |w, i| w + a[i][j] * v[i - r0]);
        a[r0][j] -= tau * w;
        for i in r0 + 1..3 {
            a[i][j] -= tau * v[i - r0] * w;
        }
    }
}

/// `A[r0.., c0..] = A[r0.., c0..] (I - tau w wᵀ)`, `w = [1, v[1..]]`.
fn reflect_cols(tau: f64, v: &[f64], a: &mut M3, r0: usize, c0: usize) {
    if tau == 0.0 {
        return;
    }
    for row in &mut a[r0..] {
        let w = (c0 + 1..3).fold(row[c0], |w, j| w + row[j] * v[j - c0]);
        row[c0] -= tau * w;
        for j in c0 + 1..3 {
            row[j] -= tau * w * v[j - c0];
        }
    }
}

/// Column `j` of `a` from row `r0` down.
fn column(a: &M3, j: usize, r0: usize) -> ([f64; 3], usize) {
    let mut v = [0.0; 3];
    for i in r0..3 {
        v[i - r0] = a[i][j];
    }
    (v, 3 - r0)
}

/// SVD `A = U diag(S) Vᵀ` through a QR decomposition first (`U = Q U_R`).
fn svd_mod(mut a: M3) -> (M3, [f64; 3], M3) {
    let mut taus = [0.0; 3];
    let mut vs = [[0.0; 3]; 3];
    for j in 0..3 {
        let (mut v, n) = column(&a, j, j);
        taus[j] = householder(&mut v[..n]);
        for i in j..3 {
            a[i][j] = v[i - j];
        }
        reflect_rows(taus[j], &v[..n], &mut a, j, j + 1);
        vs[j] = v;
    }
    let mut r = a;
    for (i, row) in r.iter_mut().enumerate() {
        row[..i].fill(0.0);
    }
    let (ur, s, v) = svd(r);
    let mut q = I3;
    for j in (0..3).rev() {
        reflect_rows(taus[j], &vs[j][..3 - j], &mut q, j, 0);
    }
    let mut u = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            u[i][j] = (0..3).fold(0.0, |acc, k| acc + q[i][k] * ur[k][j]);
        }
    }
    (u, s, v)
}

/// SVD `A = U diag(S) Vᵀ`, `S` descending: bidiagonalisation, then Golub-Kahan iterations.
fn svd(mut a: M3) -> (M3, [f64; 3], M3) {
    let (mut tu, mut tv) = ([0.0; 3], [0.0; 2]);
    for j in 0..3 {
        let (mut v, n) = column(&a, j, j);
        tu[j] = householder(&mut v[..n]);
        for i in j..3 {
            a[i][j] = v[i - j];
        }
        reflect_rows(tu[j], &v[..n], &mut a, j, j + 1);
        if j + 1 < 3 {
            let mut row = [0.0; 3];
            let n = 2 - j;
            row[..n].copy_from_slice(&a[j][j + 1..]);
            tv[j] = householder(&mut row[..n]);
            a[j][j + 1..].copy_from_slice(&row[..n]);
            reflect_cols(tv[j], &row[..n], &mut a, j + 1, j + 1);
        }
    }
    let mut d = [a[0][0], a[1][1], a[2][2]];
    let mut e = [a[0][1], a[1][2]];
    let mut u = I3;
    for j in (0..3).rev() {
        let (mut w, n) = column(&a, j, j);
        w[0] = 1.0;
        reflect_rows(tu[j], &w[..n], &mut u, j, j);
    }
    let mut v = I3;
    for j in (0..2).rev() {
        let mut w = [1.0, 0.0, 0.0];
        let n = 2 - j;
        w[1..n].copy_from_slice(&a[j][j + 2..]);
        reflect_rows(tv[j], &w[..n], &mut v, j + 1, j + 1);
    }
    bidiagonal_svd(&mut d, &mut e, &mut u, &mut v);
    (u, d, v)
}

/// Rotation `(c, s)` with `c a - s b = r`, `s a + c b = 0` (Golub & Van Loan, Algorithm 5.1.3).
fn givens(a: f64, b: f64) -> (f64, f64) {
    if b == 0.0 {
        (1.0, 0.0)
    } else if b.abs() > a.abs() {
        let t = -a / b;
        let s = 1.0 / (1.0 + t * t).sqrt();
        (s * t, s)
    } else {
        let t = -b / a;
        let c = 1.0 / (1.0 + t * t).sqrt();
        (c, c * t)
    }
}

fn rotate_cols(m: &mut M3, i: usize, j: usize, c: f64, s: f64) {
    for row in m {
        let (x, y) = (row[i], row[j]);
        row[i] = c * x - s * y;
        row[j] = s * x + c * y;
    }
}

/// SVD of the upper bidiagonal `diag(d) + superdiag(e)`, accumulated into `u` and `v`; leaves the
/// singular values in `d`, non-negative and descending.
fn bidiagonal_svd(d: &mut [f64; 3], e: &mut [f64; 2], u: &mut M3, v: &mut M3) {
    for _ in 0..300 {
        for i in 0..2 {
            if e[i].abs() <= f64::EPSILON * (d[i].abs() + d[i + 1].abs()) {
                e[i] = 0.0;
            }
        }
        let mut b = 2;
        while b > 0 && e[b - 1] == 0.0 {
            b -= 1;
        }
        if b == 0 {
            break;
        }
        let mut a = b - 1;
        while a > 0 && e[a - 1] != 0.0 {
            a -= 1;
        }
        if b - a == 1 && d[a] != 0.0 && d[b] != 0.0 {
            let (u2, s2, v2) = svd2(d[a], e[a], d[b]);
            (d[a], d[b], e[a]) = (s2[0], s2[1], 0.0);
            rotate_2x2(u, a, b, u2);
            rotate_2x2(v, a, b, v2);
        } else {
            qr_step(d, e, a, b, u, v);
        }
    }
    sort_descending(d, u, v);
}

/// Columns `a`, `b` of `m` times the 2×2 `m2`.
fn rotate_2x2(m: &mut M3, a: usize, b: usize, m2: [[f64; 2]; 2]) {
    for row in m {
        let (x, y) = (row[a], row[b]);
        row[a] = x * m2[0][0] + y * m2[1][0];
        row[b] = x * m2[0][1] + y * m2[1][1];
    }
}

/// Makes the singular values non-negative (negating the column of `v`) and sorts them descending
/// with the columns of `u` and `v`; ties keep their order.
fn sort_descending(d: &mut [f64; 3], u: &mut M3, v: &mut M3) {
    for j in 0..3 {
        if d[j] < 0.0 {
            d[j] = -d[j];
            for row in v.iter_mut() {
                row[j] = -row[j];
            }
        }
    }
    for i in 0..3 {
        let mut m = i;
        for k in i + 1..3 {
            if d[k] > d[m] {
                m = k;
            }
        }
        if m != i {
            d.swap(i, m);
            for row in u.iter_mut().chain(v.iter_mut()) {
                row.swap(i, m);
            }
        }
    }
}

/// SVD of the 2×2 block `[[f, g], [0, h]] = U diag(S) Vᵀ`, `S` signed.
fn svd2(f: f64, g: f64, h: f64) -> ([[f64; 2]; 2], [f64; 2], [[f64; 2]; 2]) {
    let tau = ((h - f) * (h + f) + g * g) / (2.0 * f * g);
    let r = hypot(1.0, tau);
    let t = if tau >= 0.0 {
        1.0 / (tau + r)
    } else {
        -1.0 / (-tau + r)
    };
    let c = 1.0 / hypot(1.0, t);
    let s = t * c;
    let w0 = (c * f - s * g, -s * h);
    let w1 = (s * f + c * g, c * h);
    let (cu, su, sv, v) = if hypot(w0.0, w0.1) > hypot(w1.0, w1.1) {
        let (cu, su) = givens(w0.0, w0.1);
        let sv = [cu * w0.0 - su * w0.1, su * w1.0 + cu * w1.1];
        (cu, su, sv, [[c, s], [-s, c]])
    } else {
        let (cu, su) = givens(w1.0, w1.1);
        let sv = [cu * w1.0 - su * w1.1, su * w0.0 + cu * w0.1];
        (cu, su, sv, [[s, c], [c, -s]])
    };
    ([[cu, su], [-su, cu]], sv, v)
}

/// One implicit-shift QR sweep over the unreduced block `a..=b` (Wilkinson shift of the trailing
/// 2×2 of `BᵀB`).
fn qr_step(d: &mut [f64; 3], e: &mut [f64; 2], a: usize, b: usize, u: &mut M3, v: &mut M3) {
    let ta = d[b - 1] * d[b - 1] + if b >= a + 2 { e[b - 2] * e[b - 2] } else { 0.0 };
    let tb = d[b] * d[b] + e[b - 1] * e[b - 1];
    let tab = d[b - 1] * e[b - 1];
    let dt = (ta - tb) / 2.0;
    let sign = if dt >= 0.0 { 1.0 } else { -1.0 };
    let mu = tb - tab * tab / (dt + sign * hypot(dt, tab));
    let (mut y, mut z) = (d[a] * d[a] - mu, d[a] * e[a]);
    let mut bulge = 0.0;
    for k in a..b {
        let (c, s) = givens(y, z);
        if k > a {
            e[k - 1] = c * e[k - 1] - s * bulge;
        }
        let (dk, ek) = (d[k], e[k]);
        d[k] = c * dk - s * ek;
        e[k] = s * dk + c * ek;
        bulge = -s * d[k + 1];
        d[k + 1] *= c;
        rotate_cols(v, k, k + 1, c, s);
        let (c, s) = givens(d[k], bulge);
        d[k] = c * d[k] - s * bulge;
        let (ek, dk1) = (e[k], d[k + 1]);
        e[k] = c * ek - s * dk1;
        d[k + 1] = s * ek + c * dk1;
        if k + 1 < b {
            bulge = -s * e[k + 1];
            e[k + 1] *= c;
        }
        rotate_cols(u, k, k + 1, c, s);
        (y, z) = (e[k], bulge);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_an_exact_parabola() {
        let x = [100.0, 110.0, 120.0];
        let c = quadratic(x, x.map(|x| 3.0 + 2.0 * x + 0.5 * x * x));
        for (got, want) in c.iter().zip([3.0, 2.0, 0.5]) {
            assert!((got - want).abs() < 1e-6 * want.abs().max(1.0), "{c:?}");
        }
    }

    /// Replays black-box cases of the old program (`docs/reverse-engineering/gsl-clean-room.md`,
    /// `tools/gsl-oracle/make-cases.py`):
    /// `MD_QUADFIT_CASES=cases.txt cargo test -p md-core quadfit -- --ignored --nocapture`, one case
    /// per line as nine IEEE bit patterns in hex: x0 x1 x2 y0 y1 y2 c0 c1 c2. The cases stay local.
    /// Every case with the default 10 bins per scan must match; other bin counts differ in about 1
    /// of 10 000 design matrices by 1–2 ulp (known, counted only; issue 290).
    #[test]
    #[ignore]
    fn matches_the_old_fit_bit_for_bit() {
        let path = std::env::var("MD_QUADFIT_CASES").expect("MD_QUADFIT_CASES");
        let text = std::fs::read_to_string(path).unwrap();
        let (mut n, mut bad, mut other_bins) = (0, 0, 0);
        for line in text.lines() {
            let v: Vec<f64> = line
                .split_whitespace()
                .map(|h| f64::from_bits(u64::from_str_radix(h, 16).unwrap()))
                .collect();
            let c = quadratic([v[0], v[1], v[2]], [v[3], v[4], v[5]]);
            n += 1;
            if c.map(f64::to_bits) != [v[6], v[7], v[8]].map(f64::to_bits) {
                if v[1] - v[0] != 10.0 || v[2] >= 9.4e7 {
                    other_bins += 1;
                    continue;
                }
                bad += 1;
                if bad <= 5 {
                    eprintln!("mismatch: {line} -> {c:?}");
                }
            }
        }
        eprintln!(
            "{n} cases, {bad} differ, {other_bins} known differences (other bins, x² ≥ 2^53)"
        );
        assert!(n > 0 && bad == 0, "{bad} of {n} cases differ");
        // About 1-2 in 10 000 with make-cases.py; more means a regression outside the 10-bin path.
        assert!(
            other_bins * 1000 <= n,
            "{other_bins} known differences in {n} cases"
        );
    }

    #[test]
    fn hypot_is_close_and_symmetric() {
        for (x, y) in [(3.0, 4.0), (1.0, 1e-9), (0.0, 0.0), (-5.0, 12.0)] {
            assert_eq!(hypot(x, y), hypot(y, x));
            assert!((hypot(x, y) - (x * x + y * y).sqrt()).abs() <= 1e-15 * hypot(x, y));
        }
    }
}
