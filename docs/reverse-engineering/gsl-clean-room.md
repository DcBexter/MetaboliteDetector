# Clean-room rebuild of the old quadratic fit (issue 250)

The old program finds the precise apex of a peak (algorithms.md §3.6) with
`tools::Statistics::polynomialRegression(X, Y, 2)`, which calls `gsl_multifit_linear` of the GNU
Scientific Library, statically linked into the binary. At a half-bin tie or a flat top, the last bits
of that fit decide the result, so the port must reproduce them. GSL is GPL and this project is MIT,
so `crates/md-core/src/quadfit.rs` rebuilds the fit without GSL code. This note records how, so the
boundary can be checked.

## The boundary

Decided in issue 250:

- **Allowed:** reverse engineering the old author's own code from the binary:
  `determinePreciseMaximizationTime` (all template variants), `polynomialRegression`, units and
  scaling. These were decompiled with Ghidra (REA, WSL).
- **Not allowed:** reading or transcribing GSL internals. No GSL source was opened. No GSL function
  in the binary was decompiled or disassembled.
- **GSL as a black box:** the oracle below calls GSL functions of the old binary through their
  public API (prototypes from the GSL reference manual). It feeds inputs and records outputs bit
  for bit. The rebuild is checked against those outputs only.

## What the old author's code does (decompiled, allowed)

- `determinePreciseMaximizationTime`: `X_i = bins * (maxScan - 1 + i)` in absolute scans,
  `Y_i` = the peak values at those scans. Then `c = polynomialRegression(X, Y, 2)` and
  `t = (double)(int)(c1 / (c2 * -2.0) + 0.5) / bins`. It is set only if `|t - maxScan| < 1`.
  Otherwise, and at the peak edge, `interpolatedMaxScan` keeps its initial 0.
- `polynomialRegression`: the matrix `[1, x, std::pow(x, j)]`, then `gsl_multifit_linear` and
  the coefficients `c`. `pow` is glibc's, dynamically linked.

## The oracle

`docs/reverse-engineering/tools/gsl-oracle/`:

- `gsl_oracle.cpp` is an `LD_PRELOAD` driver like the golden driver. It replaces `main` and calls
  functions of the old binary by their `nm` addresses: `polynomialRegression`,
  `gsl_multifit_linear`, `gsl_linalg_balance_columns`, `gsl_linalg_SV_decomp(_mod)`,
  `gsl_linalg_bidiag_decomp`, `gsl_linalg_householder_transform`/`_hm`/`_mh`,
  `gsl_linalg_SV_solve`, `cblas_dnrm2`/`ddot`, and the binary's internal hypot helpers. It prints
  each result double as 16 hex digits.
- `run.sh` builds the driver (it checks the binary's SHA-256 first) and runs it on stdin.
- `make-cases.py` writes random fit cases for the ignored Rust test
  `quadfit::tests::matches_the_old_fit_bit_for_bit`.

The outputs come from the old program. They stay local and are never committed.

## How the steps were found

Each step started from a textbook method: Golub & Van Loan, *Matrix Computations*, and the
reference BLAS. Each was accepted only when it reproduced the oracle bit for bit on many inputs.

1. **The fit is a plain least-squares solve.** `fit` equals `mfl` on the same matrix. Only `X`
   enters the decomposition; `y` enters only the solve.
2. **Column balancing** (`bal`): each column is divided by a power of two `f`, with the column's
   absolute sum over `f` in [0.5, 1].
3. **Householder reflections** (`hh`, `hm`, `mh`): the tail norm is reference-BLAS `dnrm2` (running
   scale). `beta = -sign(alpha) * hypot(alpha, xnorm)`, `tau = (beta - alpha) / beta`, tail scaled
   by `1 / (alpha - beta)`. They are applied column by column (`hm`) and row by row (`mh`).
4. **hypot:** the binary imports `hypot` from glibc. A correctly rounded hypot differs on 318 of
   185 000 calls. The correction step of Borges (2019, arXiv:1904.09481, form without FMA) matched
   all 570 197 calls on the fit path and 200 000 random pairs (1e-150 to 1e150). The port writes
   it out because MSVC's `hypot` differs in the last bit.
5. **Bidiagonalisation** (`bid`): alternating column and row reflections.
6. **`SV_decomp_mod`** (`svdm` against `svd`): a QR decomposition first, then the SVD of `R`, and
   `U = Q U_R`.
7. **The SVD of the bidiagonal:** Golub-Kahan iteration (G&VL §8.6). Deflation where
   `|e_i| <= eps (|d_i| + |d_{i+1}|)`. A Wilkinson shift from the trailing 2×2 of `BᵀB`. A
   bulge chase with the Givens rotation of G&VL Algorithm 5.1.3. A 2×2 block is solved in closed
   form: the Jacobi angle on `BᵀB`, then `U` from the larger column (2 926 of 2 926 random 2×2
   cases). At the end, signs are made non-negative and a selection sort orders the values
   descending.
8. **The solve:** `t = Uᵀy`, then `c_i = (Σ_j V_ij (1 / S_j) t_j) / D_i` for `S_j > eps S_0`.

## Verification

- SVD of the balanced design matrix, bit for bit:
  - 10 bins per scan, every scan 1 to 1 000 000: no difference.
  - 1, 2, 3, 4, 5, 7, 20, 25, 50 and 100 bins, scans 1 to 50 000: no difference.
- Whole fit:
  - 300 269 cases, among them the 269 real precise-apex ties of the 15 example runs (217 half-step
    ties rounded down, 50 rounded up, 2 flat tops).
  - With 10 bins and `x² < 2^53`, every case matches.
  - The 269 ties give the old apex.
- Golden test:
  - 246 of 246 components (default detector).
  - 2 807 of 2 807 components (matched filter).
  - Checked against a matched-filter dump of the old peaks (golden driver with `MD_MATCHED`) for
    eg07, eg08, eg09, eg15 and eg25.

## Known limits (issue 290)

- **Other bin counts:** about 1 in 10 000 design matrices differs by 1–2 ulp (11 of 99 900 for 1
  to 100 bins and scans 1 to 999). It differs inside the QR sweeps; the path and the deflation
  decisions are the same. Random 3×3 bidiagonals match only about 70 % of the time, so a rounding
  detail of the sweep is still unknown. It matters only together with a half-bin tie.
- **Very large x:** `x * x` stands in for glibc `pow(x, 2)`. The two are equal while `x² < 2^53`,
  i.e. while bins × scan stays below 9.4e7.

## Reproduce

```sh
# WSL with the old .deb installed:
docs/reverse-engineering/tools/gsl-oracle/make-cases.py 100000 > cases.txt
MD_QUADFIT_CASES=cases.txt cargo test -p md-core quadfit -- --ignored --nocapture
```
