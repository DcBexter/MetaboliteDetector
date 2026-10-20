# MetaboliteDetector 3.3 – Feature inventory

Source: the shipped Linux build `MetaboliteDetector-3.3.20200313-Linux_0.deb` (local `Old/`, not committed).
The ELF `usr/bin/metabolitedetector_app` (SHA-256 `7883d70c…4794`) is a Qt5/C++ GCC build with DWARF
debug info and full symbols. The source code is lost.

Evidence types: **S** = symbol/class names (`nm -C`), **T** = UI/log strings in the binary, **D** = Ghidra
decompilation, **P** = package contents. Algorithm details: [algorithms.md](algorithms.md). File layouts:
[file-formats.md](file-formats.md).

## Package and dependencies (P)

- Depends on Qt5 (core, gui, printsupport, concurrent), zlib, `r-cran-rcpp`; ships `RInside_0.2.13.tar.gz`
  (embedded R console), Qwt and QCustomPlot (plots), Boost (regex), ALGLIB (linear algebra).
- Ships one sample library `usr/share/metabolitedetector/lib_aa.lbr` (gzip, custom binary; amino acids).
- Maintainer `k.hiller@tu-bs.de`, TU Braunschweig.

## Workflow and features

| F | Feature | Old entry point (T/S) | Notes | New app |
|---|---------|----------------------|-------|---------|
| 1 | NetCDF (ANDI-MS) import | `&NetCDF import...`, `CDFBatchImporter`, `gcms::CdfImporter` | Batch converts `.cdf` to the internal `.bin`/`.idx` format; optional spectrum compression; optional TOF calibration. | done (read directly, no conversion; NetCDF export #30) |
| 2 | FastFlight import | `&FastFlight Import...`, `gcms::FastFlightImporter`, `MsData.FFT` | Raw TOF data; needs `.mzcalib` TOF calibration; can write `_centroid/_profile/_nominal.cdf`. | done: raw import (#29), centroid/nominal CDF export and TOF calibration window (#148); automatic calibration and compress-on-import not ported |
| 3 | TOF m/z calibration | `TOFCalibrationWindow`, `TOFCalibrationCalculator` | Polynomial fit of calibration pairs from a reference spectrum; saved as `.mzcalib`. | done: `.mzcalib` read/write and the fit in `md_core::tof` (#29), calibration window from reference ions (#148); automatic calibration not ported |
| 4 | Chromatogram viewer | `ChromatogramPlot`, `SpectrumPlot`, `MainWindow` | TIC, single-ion chromatogram, scan/time/ion spinners, slider, 1st-derivative plot, zoom, CSV/clipboard/image export, label compounds by score (≤0.7, >0.7, >0.9). | done: TIC/SIC, zoom, slider, unload, derivative, time unit, run info, score filter, auto-Y (#79), CSV and PNG (#9); PDF/SVG/clipboard #204 |
| 5 | Spectrum view | `SpectrumPlot`, `SplitSpectrumPlot` | Profile/centroid, normalisation (none, base peak, unit, relative), split/difference mode for detected vs library spectrum. | done: normalisation (#9), library comparison and component panel (#77, #126); no profile mode |
| 6 | Chromatogram overlay | `ChromatogramOverlayDialog`, `MultiChromatogramPlot` | Overlay TIC or SIC of several runs. | done (#8) |
| 7 | Component detection (deconvolution) | `&Detect Components...`, `GCMSAnalyzer`, `PeakAnalyzer`, `DefaultPeakDetector`, `MatchedPeakDetector`, `ModelPeak`, `BaselineStatistics` | Peak detection per ion (derivative thresholds or matched filtering), baseline adjustment, noise units, binning (bins/scan), deconvolution width, min peaks, min % of base peak, model-peak integration. Results stored as `.cmp`. | done: default and matched detector (#82) + deconvolution; `.cmp` read/write (#13), save/load and `.cmp` as batch input (#107) |
| 8 | Library identification | `gcms::LibrarySearch`, `SpectrumFunctions`, settings `LS_*` | Score = spectrum similarity, RI similarity or combined; pure/impure factors, cutoff, max ΔRI, mass filter ranges, scaled library, fragment pre-filter (`SQL_*`). | done (scaled library always on); `SQL_*` pre-filter in the settings (#131); `LS_IDENTIFICATION_ALGORITHM` never set in the old app, not ported |
| 9 | RI calibration wizard | `RICalibrationWizard`, `RICalculator` | Alkane (or other) calibration mix from a library, RT→RI table per chromatogram (`.calib`), optional "advanced" RT shift correction using a calibration compound; RT-only alignment as fallback. | basics + reference shift done (#20, #45); mix from library, editable table, per-run calibration, RT-only (#81) |
| 10 | Batch compound detection wizard | `Compound Detection Wizard...` | Runs #7 for many files after RI calibration. | done: single run (#78, #107), many runs to `.cmp` in the batch (#205) |
| 11 | Batch quantification | `BatchQuantificationWizard`, `BatchQuantificator` | Targeted (reference library only) or non-targeted (temporary library collecting all spectra); replicate groups via regex; compound filter (reproducibility, discrepancy index, min ions); extended SIC rescan for missing quant ions; result table export. | done: targeted + groups (#22), non-targeted, rescan, S/N and discrepancy filters (#58, [batch-quantification.md](batch-quantification.md)), result views (#83, #160); ANOVA column #203 |
| 12 | Quantification ions | `Compound::assignQuantificationIons`, `QuantIonChecker`, `QUANT_*` | N ions per compound, min distance, quality cut-off, excluded ions; manual/automatic ions in the library. | done: automatic ions in batch, store in library (#80), editor (#12); automatic ions in the editor (#208) |
| 13 | Normalisation | `CompoundSeriesNormalizer`, `PoolNormalizer` | By internal standard compound, summed sample signal, pool measurement, z-score. | done (#23, #206) |
| 14 | Statistics | `BiomarkerDialog`, `tools::Statistics`, `BoxPlot`, `HistogramPlot` | Two-group biomarker table (means, ratio, p-value, hits), independent/paired samples, box plots. | done (#24, #189); ANOVA column #203 |
| 15 | PCA | `tools::PCA`, `PCAPlot` | Score plot with selectable PCs, loadings ("impact of metabolites"), cumulative variance. | done (#24) |
| 16 | MID wizard | `MIDWizard`, `tools::MIDSolver`, `NaturalMIDCalculator`, `MIDCompoundSeries` | Mass isotopomer distributions of labelled vs unlabelled runs; correction by sum formula, by unlabelled reference runs, by library spectrum, or uncorrected; M1 correction; tracer natural abundance; fragment syntax `(174/8)(261/5)`. | done (#28, #141, #156): four corrections, M1, tracer, result tree, bar plots, group statistics with StdErr, text exports with clipboard copy, FRAG export; column/diagram order #207; live test with labelled runs #147 |
| 17 | Natural isotope calculator | `NaturalIsotopeWidget`, `:/isotopes.csv` | Isotope pattern for a sum formula. | done (#27) |
| 18 | Library editor | `LibraryEditor`, `DialogAddLib`, `DialogImportLibrary`, `DialogExportLibrary` | Edit compounds, RI, quant ions, features (FRAG, comments); import/export AMDIS MSL (nominal/profile); add detected compound to library; copy RTs to a target library. | done: `.lbr` v2 to v5 read (#41), MSL import (#11), MSL export nominal / profile / FWHM (#84, #114), editor (#12; details #208), add compound / copy RTs (#80); "Modify entry" spectrum merge not ported |
| 19 | NIST search | `&NIST search...`, `gcms::NISTParser` | Search a user-supplied NIST library file. | NIST text `.sdf` with the old fragment pre-filter and NIST MSP `.msp` read and searched (PR #123); live test with a real NIST file open (#21) |
| 20 | R console | `RConsole`, `MDRConsole`, RInside | Results as R `data.frame data`, scripts, SVG/PDF plots, `.RData` export. | done: CSV + R templates in `scripts/r/` (#26); embedded R not ported on purpose |
| 21 | Exports | `NetCDFExporter`, `CSVChromatogramExporter`, `CompoundExporter` | Chromatogram to NetCDF/CSV, compound tables to CSV, plots to PDF/PNG. | done: NetCDF (#30), CSV/PNG (#9), compound CSV (#78); PDF/SVG/clipboard #204 |

## Parity with the new app

Second audit of 2026-10-20 against `main` (e7c390e). "done" means the new app covers it; `#N` is the
open issue that owns the rest. Sub-features come from the binary's menus, dialogs and settings.

| F | Sub-feature | New app |
|---|-------------|---------|
| 1 | Batch conversion, spectrum compression, recursive search | not needed: runs are read directly |
| 2, 3 | FastFlight `.FFT`/`.FFC` import, centroiding (`SF_*`), CDF export, TOF calibration window, `.mzcalib`, Calibration Setter | done (#29, #148); automatic TOF calibration not ported |
| 4 | First derivative, time unit (ms/s/min), run info, label filter by score, auto-scale Y | done (#79, #111) |
| 4, 21 | Chromatogram/spectrum CSV, PNG | done (#9) |
| 4, 14–16, 21 | PDF/SVG export, copy image/data to clipboard, image size, PNG of MID plots | #204 |
| 5 | Split/difference/normal comparison with the library, component panel (S/N, sharpness, peaks) | done (#77, #126) |
| 5 | Profile (centroid toggle): the rewrite reads centroid spectra only | not planned |
| 7 | Matched-filter peak detection for components | done (#82) |
| 7, 10 | Results of a run: table, find compound, CSV, save/load `.cmp`, `.cmp` as batch input | done (#78, #107, #169) |
| 10 | Compound Detection Wizard: detect many runs and save a `.cmp` for each | done (#205) |
| 8, 19 | NIST text/MSP import, fragment pre-filter, `SQL_*` settings, second NIST library, NIST search for one component | done (PR #123, #131); NIST path not kept between sessions |
| 8, 19 | Live test with a real (licensed) NIST library | #21 (needs-human-input) |
| 9 | Calibration mix from a library, editable RT/RI table, RT-only fallback, calibration per run | done (#81) |
| 11 | Non-targeted batch, rescan for missing quant ions, excluded ions, discrepancy index, S/N and required-peaks filters | done (#58, #110, #124) |
| 11 | Result tree, alternatives, RT table, rel. std. error, blanks / P(blank), delete rows/columns, total intensities, all-ions CSV | done (#83, #160, #190); cell details, full settings file, R script #193; ANOVA p-value #203 |
| 12, 18 | Store detected compound and quant ions in the library, copy RTs to a target library | done (#80) |
| 18 | Editor: compounds, RI, quant ions, features dialog, `FRAG` | done (#12) |
| 18 | Editor: automatic quant ions, ion frequency, resolution, v5 MID fragments (`formulas`, `mid_sizes`) | done (#208) |
| 13 | Internal standard, summed signal, pool, z-score | done (#23, #188) |
| 14, 15 | Group assignment, paste sample list, p-value filter, box plots, PCA axes/loadings/variance, PNG | done (#24, #187, #189) |
| 16 | MID wizard: methods, tracer, M1, result tree, StdErr, CSV/clipboard, MID plots, group statistics | done (#28, #141, #156, #183) |
| 16 | Change column order, change order of diagrams | #207 |
| 16 | Check with real labelled runs | #147 (needs-human-input) |
| 20 | CSV results + R template instead of the embedded R console | done (#26); manual check in Excel and R #154 (needs-human-input) |
| files | `.cmp` v4, `.calib` text and binary, `.lbr` v2–v5, `.mzcalib`, `.bin`/`.idx` | done (#13, #29, #30, #41) |
| files | `.cmp` versions 1 to 3 (older releases, no sample) | not planned |
| files | High-resolution runs: exact-mass traces, spectra, `reorderPeaks` | done (#34, #173); exact masses in `.cmp` and search, real FastFlight noise factor #198 (needs-human-input) |

Open bugs and acceptance outside parity: #196, #200 (UI tests and layout), #129 (ruleset), #130 (manual
acceptance with sample data).

## Internal file types (T)

`.bin` + `.idx` converted chromatogram · `.cmp` detected compounds · `.lbr` compound library (gzip) ·
`.calib` RI calibration · `.mzcalib` TOF calibration · `.msl` AMDIS library (import/export) ·
NIST library (import for search) · settings stored by `GCMSSettings::toDisk`.

## Observations about the old code quality (T/S)

- Long-running work on Qt threads with ad-hoc locks (`CRITICAL: Was not able to obtain lock in MainWindow::setLibrary()`).
- Mixed German/English UI text (`Keine Lib`), hard-coded home paths (`/home/kh`).
- R integration via RInside makes Linux packaging depend on a matching R installation – the main
  reason the app does not run on Windows.
