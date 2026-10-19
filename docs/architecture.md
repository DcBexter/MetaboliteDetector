# Architecture

MetaboliteDetector is being rebuilt in Rust as a native desktop app for Windows, Linux and macOS.
The old Qt5/C++ version (3.3, 2020) is the behavioural reference; its features are listed in
[reverse-engineering/feature-inventory.md](reverse-engineering/feature-inventory.md).

## Decisions

| Topic | Decision | Why |
|-------|----------|-----|
| Language | Rust (edition 2024) | Memory-safe, fast numerics, one static binary per OS, no runtime to install. |
| GUI | [egui](https://github.com/emilk/egui)/eframe + egui_plot | Pure Rust, single executable, interactive plots that handle large series. Chosen by the maintainer 2026-10-19. |
| File dialogs | rfd | Native dialogs on every OS. |
| NetCDF | `netcdf3` (pure Rust) | ANDI-MS files are NetCDF-3; no libnetcdf/HDF5 C dependency to ship on Windows. |
| Parallelism | rayon | Data-parallel loops over scans, ions and files. |
| R console | Not ported directly | RInside tied the old app to a Linux R install. Replace with CSV/R-script export first; optional embedded scripting later. |

## Layout

```
Cargo.toml              workspace, release profile (LTO, 1 codegen unit)
src/md-core/            library: data model, import/export, analysis algorithms (no GUI deps)
  src/lib.rs            Chromatogram (flat ANDI-style arrays), TIC/SIC, plot decimation
  src/andi.rs           ANDI-MS NetCDF import with input validation
  examples/perf.rs      performance smoke test on synthetic large runs
src/md-app/             egui desktop app (binary `metabolite-detector`), drawing only
docs/reverse-engineering/  evidence from the old binary (features, algorithms, file formats)
```

Planned additions: `src/md-cli` (headless batch runs for big studies) and modules in `md-core` for
peak detection, deconvolution, library search, RI calibration, quantification, MID and statistics.
All analysis code lives in `md-core` so the GUI, CLI and tests share it.

## Performance rules

Large studies (hundreds of runs, tens of millions of points each) are the normal case.

1. **Compact storage.** Spectra are flat `f32` arrays with scan offsets (≈ 8 bytes per data point),
   matching the old `GCMSDiskScan<int, float>`. Accumulate sums in `f64`.
2. **Never block the UI.** File loading and every analysis step run on worker threads; the UI polls
   results and shows progress. Long jobs must be cancellable.
3. **Draw proportional to pixels, not data.** Line plots go through `md_core::decimate`
   (min/max per pixel bucket of the visible range), so peaks stay visible at any zoom.
4. **Parallel by default.** Per-scan, per-ion and per-file work uses rayon.
5. **Index once, query many.** m/z is sorted within each scan (binary search for SICs). Deconvolution
   will build an ion-major index once per run instead of rescanning spectra per ion.
6. **Measure.** `examples/perf.rs` reports timings; each new algorithm adds a case, and PRs that
   touch hot paths report before/after numbers against the table below.
7. Memory-mapped or streamed access (like the old disk scan) is added when real runs exceed RAM;
   the cut-off is recorded as a `shortcut:` comment in `andi.rs`.

Measured 2026-10-19, release build, 36 000 scans × 500 ions (18 M points, 137 MiB):

| Operation | Time |
|-----------|------|
| TIC | 2.0 ms |
| SIC (one m/z) | 1.0 ms |
| Decimate TIC for drawing | 0.05 ms |

Real runs (TargetSearchData, 4 400 scans, ~745 000 points, 6 MB `.cdf`): read + TIC 13–25 ms per
file; all 15 runs in 0.27 s (`cargo test --release -p md-core --test example_data -- --nocapture`).
The whole file is read into memory first because `netcdf3` seeks per record (1.8 s per file on a
plain `File`).

## Build and test

```sh
cargo build --release            # binary: target/release/metabolite-detector(.exe)
cargo test --workspace
cargo clippy --workspace --all-targets
cargo run --release -p md-core --example perf
cargo run --release -p metabolite-detector -- path/to/run.cdf
example-data/fetch.sh             # Windows: .\example-data\fetch.ps1
```

## Compatibility goals

- Read the old binary files students already have: `.lbr` libraries first, then `.cmp`, `.calib`,
  `.bin`/`.idx` (see [reverse-engineering/file-formats.md](reverse-engineering/file-formats.md)).
- Reproduce the old numerical results (peaks, compounds, scores, areas, MIDs) within a documented
  tolerance, verified against output of the old app on the same input.
