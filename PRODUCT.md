# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

<!-- Native desktop app built with egui/eframe (Windows, Linux, macOS). "web" is the closest
schema value: egui draws its own widgets on every OS, there is no iOS/Android design language. -->

## Users

Students and researchers in GC-MS metabolomics labs. Two groups weigh equally:

- **Experienced users of MetaboliteDetector 3.3** (the old Qt app, TU Braunschweig). They know the
  workflow, the old setting keys (`AN_*`, `LS_*`, `QUANT_*`) and the old file formats, and expect the
  same results and fast, dense screens.
- **New students** who meet the tool for the first time, often for one thesis study. They need every
  control to say what it does without reading a manual.

Typical session: open a study (one to several hundred runs), look at chromatograms and spectra,
tune peak detection and deconvolution, identify substances against a library with retention index,
quantify all runs in a batch, export the table for statistics.

## Product Purpose

Rewrite of the lost-source MetaboliteDetector 3.3: deconvolution, identification and quantification
of GC/MS metabolomics data, as one native binary per OS without R or Qt installs. Success means full
feature parity with 3.3 (see `docs/reverse-engineering/feature-inventory.md`), identical numbers where
the old algorithm is known, and responsive work on huge studies.

## Positioning

The only maintained continuation of MetaboliteDetector's algorithms (reverse-engineered from the
3.3 binary), reading its old files (`.lbr`, `.bin/.idx`, `.calib`, settings) directly, and fast enough
that hundreds of runs stay interactive.

## Operating Context

- Desktop or laptop at the lab or office, long sessions, often two monitors; mouse and keyboard.
- Input files: ANDI-MS NetCDF (`.cdf`), old `.bin/.idx` caches, libraries (`.mdlib`, `.lbr`, AMDIS
  `.msl`), RI calibration (`.calib`), the old settings file.
- Output: CSV tables loaded into R (`scripts/r/load_results.R`) or Excel, PNG plots for theses and
  papers, NetCDF and library files.
- Headless big studies go through `md-cli`; the app is for looking, tuning and checking.

## Capabilities and Constraints

- Chromatogram viewer (TIC/SIC, overlay of many runs), spectrum view with normalisation, component
  detection, library search with RI, RI calibration, batch quantification, library editor, natural
  isotope calculator, exports. Remaining parity work is owned by the issues in the feature inventory.
- **Performance:** loading and analysis never block the UI; plots draw proportional to pixels, not
  data (`md_core::decimate`).
- **UI language:** English and German (S8, #257): the system language decides, switchable in Help
  and Settings. Labels stay short, complete strings (no sentence assembled from fragments), each
  written as an English/German pair (`lang::t`, `lang::tf!`). CSV headers and file formats stay English.
- **Setting names:** each setting shows a plain-language name; the old 3.3 key (e.g.
  `AN_PEAK_THRESHOLD_BEGIN`) stays reachable on hover so 3.3 users and settings files still match.
- **Explanations everywhere:** every button, field and setting carries a short hover explanation of
  what it does (and for settings, what the old key was). Users asked for this explicitly.
- Analysis code lives in `md-core`; `md-app` only draws.

## Brand Commitments

- Name: **MetaboliteDetector** (one word, as in 3.3).
- Score colours for library hits follow 3.3: above 0.9, above 0.7, the rest (green, orange, red).

## Evidence on Hand

- 15 public example runs (`example-data/fetch.ps1`, not in git).
- Screenshot of a 15-run overlay: `docs/images/overlay-15-runs.png`.
- Reverse-engineering notes of 3.3: `docs/reverse-engineering/`.
- No user studies, testimonials or usage numbers exist; do not invent them.

## Product Principles

1. **Parity first, then better.** Every 3.3 workflow exists and produces the same numbers before it
   gets reshaped.
2. **Self-explaining for newcomers, dense for experts.** Defaults and hover explanations serve the
   student; compact layouts and old key names serve the 3.3 user. Neither is sacrificed.
3. **Never block, never lose.** Long work runs in the background with visible progress; files are
   written atomically; unsaved edits are guarded.
4. **Data first.** Plots and tables get the space; chrome stays out of the way.

## Accessibility & Inclusion

No product-specific standard is set. Hover explanations must not be the only way to learn critical
state (errors and warnings are shown inline). Score colours always come with the number.
