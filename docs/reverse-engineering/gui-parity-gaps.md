# GUI parity gaps: MetaboliteDetector 3.3 vs. the rewrite

Date: 2026-10-20. Sources: an REA inventory of the old ELF (`metabolitedetector_app`, SHA-256 `7883d70c…4794`), a source inventory of `crates/md-app`, and a side-by-side check of both apps on `7235eg06`.

Rule: parity means every old **function** can be reached somewhere. It does not mean the old layout or the old Qt widgets. See `PRODUCT.md` ("Parity first, then better").

## 1. Data accuracy (status)

| Check | Result |
|---|---|
| `golden_old_md` / `example_data` / `old_results` (old `.cmp`, `.calib`, `raw.tsv` written by the old binary) | all pass |
| GUI spot check `7235eg06` | old 18 components = new 18 components; TIC heights agree (about 0.8 M at 5.92 min) |
| Default detector | 246/246 components identical (clean-room fit, #250) |
| Matched filter | 2807/2807 identical (#250) |
| Library search scores | open: #242 (scores 0.05–0.15, few hits); #248 (TargetSearch library does not open) comes first |

## 2. Main window: what users miss first

These are the views the old main window shows at all times.

| Old element | New state | Gap |
|---|---|---|
| **Analysis Results** dock: compound-name list of all components, "No match: …" rows without a library | Only under Analysis → "Results of the active run…", and **empty without a library** (`results.rs:56-62`), although `md_core::results::result_rows` already builds "No match: N" rows | **Bug / gap**: the list must always be visible and must work without a library |
| **Detected Peaks** plot: overlaid ion profiles of the selected component over scans | Missing. The right panel shows only a peak table; a click shows one SIC in the main plot | **Gap** |
| **Spectra** plot: compound + library + centroided library in one plot | Present below the TIC (deconvoluted spectrum after a click on a component dot; Split/Difference with a library hit). In a small window it falls below the visible area | Layout: keep the spectrum visible at small window sizes |
| **Compound Information**: detected compound (RI, RT, ions, centroid ions, S/N, peak sharpness) next to library compound (RI, RT, ions, spec. sim., RI sim.) | Right panel: S/N, sharpness, discrepancy, ions, scores, peaks, hits. RT/RI of the component are not shown side by side with the library entry | Small gap: side-by-side detected vs. library |
| Spectrum dock with scan spinbox and RT field | Scan only through click or ←/→ | Small gap: direct entry of scan/RT |
| View menu: show/hide docks | Fixed panels | Replaced by a modern layout; a "View" menu with panel toggles is still useful |

Status: #251 (spec S1) closed every gap of this table: the component list sits below the runs and works
without a library, Detected Peaks is beside the spectrum, the component info shows detected and library
values side by side, scan and RT can be typed in, and **View** hides each resizable area. Saving the layout
is S4 (#253).

## 3. Menus, toolbar, shell

| Old function | New | Gap |
|---|---|---|
| File: Open (Ctrl+O), Save results (Ctrl+S), GCMS info, Unload | Open chromatograms…, Results → Save components…, Run info…, "×" in the run list | Shortcuts missing |
| Import: NetCDF import dialog (batch conversion, compression), FastFlight import (Ctrl+I, TOF calib., recursive), Import MSL library | NetCDF/FastFlight open directly; MSL opens as a library | Batch conversion of many files to the native format: check whether still needed |
| Export: results CSV (Ctrl+E), NetCDF export, chromatogram CSV, library export (MSL nominal/profile/FWHM) | Results Export CSV, Export as NetCDF, Export chromatogram/spectrum, library editor Export MSL | Library export: check that all three MSL types are there; shortcut missing |
| Tools: RI-Calibration Wizard, Compound Detection Wizard, Batch quantification wizard, MID wizard | Single windows: RI calibration, Batch (option "Save components as .cmp"), MID | Guided step flows since #256 (Analysis → Step by step) |
| Tools: Detect Components (Ctrl+D) | Runs automatically on load / after Settings Apply | Fine; shortcut not needed |
| Tools: Library Editor, Split spectrum window, Chromatogram overlay, Find Compound, Natural Isotope Calculator, Copy RT to target lib | Library editor, Split/Difference, run-list checkboxes (overlay), Results "Find substance", Natural isotopes, Copy retention times | Covered |
| Toolbar: NIST search | Button in the right panel, only after "Open NIST library" | Covered; the NIST path is also chosen in Settings (#254) |
| Toolbar: Add selected compound to library (Ctrl+L) | Library editor → add component | Shortcut missing |
| Toolbar: About | None | **Gap**: Help → About (version, licence, authors, source) |
| Quit | None (window close only) | **Gap**: File → Quit (Ctrl+Q) |
| Settings dialog: Appearance (time scale, axes labels, label highest peaks) | Time unit in the main view; Settings → Appearance: axes labels, m/z of the 5 highest peaks in view | Covered (#254) |
| Settings: Identification: Compound Lib path, NIST path, Cutoff score (0.85), Similarity-score radio (combined / spectrum / RI), Use scaled lib | Settings → Libraries (both paths), Score (combined / spectrum / RI), Cutoff score; scaled lib "stored, not used" | Covered (#254); scaled lib is dead in 3.3 too (algorithms.md §4b.4) |
| Settings: Max. RI difference shown as 10.0 in the old dialog | New RI tolerance 20 | Check: same key `LS_RI_DIFF`? The RI wizard sets it to 20 during calibration |
| Settings: Centroid (SF_*), baseline adjustment | SF_CENTROID_* shown as used (centroid export); SF_PEAK_*/SF_MAX_BASE_DEV dead in 3.3 too; SF_FWHM and baseline adjustment not wired | Gap for SF_FWHM and baseline adjustment (follow-ups of #254, algorithms.md §4b.4) |
| QSettings: window size/position, dock state, recent files, library path persist | Settings, panels, panel sizes, window, recent runs and library paths persist (not the Qt dock blobs); 3.3 values are imported once (#253) | Fine |
| Help button (F1) in Settings | Hover texts | Fine |
| Context menus, drag and drop | None | Modern additions, not parity |
| Language: English only | English only | New request: DE/EN switch |

## 4. Wizard dialogs not covered as a step flow

- **RI Calibration Wizard**: introduction → reference chromatogram (alkane mix) / RT-only option → calibration mix from library → advanced settings (calibration compound, required score) → table RT↔RI → chromatogram selection → compound detection (only if needed) → start. New: one window with all the parts; the order of steps is not shown.
- **Compound Detection Wizard**: chromatogram selection → detection settings → start. New: Batch option "Save components as .cmp".
- **Batch Quantification Wizard**: replicate groups (manual, auto-regex `_(.*)_`), blanks → targeted / non-targeted → settings → start. New: Batch window; replicate groups through a regex or a sample table.
- **MID Wizard**: method (sum formula / reference / library / uncorrected) → labeled / unlabeled chromatograms → settings. New: MID window.

## 5. Not yet covered by any issue

1. ~~Results list always visible and usable without a library (bug, section 2).~~ Done in #251.
2. ~~Detected Peaks plot (ion profiles of the selected component).~~ Done in #251.
3. ~~Detected vs. library compound side by side, scan/RT entry.~~ Done in #251.
4. Help → About, File → Quit, keyboard shortcuts.
5. ~~Persisted settings, window layout, recent files, library paths.~~ Done in #253.
6. ~~Settings scope: cutoff score, score mode, scaled lib, centroid values, appearance options.~~ Done in #254; SF_FWHM and baseline adjustment are follow-ups.
7. ~~Simple mode with guided step flows (RI calibration, detection, batch, MID) next to the expert windows.~~ Done in #256.
8. DE/EN language switch.
9. ~~GSL tie deviations (section 1): decide whether to close them with a clean-room SVD.~~ Done in #250 and #290.
10. ~~Component markers: keep the dots on the peaks, and add the old triangles on the x axis as larger click targets. Put a view toolbar above the plots, as in the old app (zoom/pan/reset, crosshair, info, marker filters by score).~~ Done in #243.
11. ~~TIC zoom stays stale after switching runs; a click does not snap to a nearby component when it is just over 3 scans away.~~ Done in #251.

Already covered elsewhere: #241, #242, #243 (zoom buttons, "?" help, MSL export, bigger Settings window), #244, #245, #246, #247, #248.
