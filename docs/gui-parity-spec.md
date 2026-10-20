# GUI parity and guided mode: spec

Status: agreed with the owner on 2026-10-20 (grilling session). Gap list: `docs/reverse-engineering/gui-parity-gaps.md`.

Principle: every function of MetaboliteDetector 3.3 can be reached. The look may be modern. Correct data comes first.

## Order

1. #248 (TargetSearch library), then #242 (library-search scores).
2. Clean-room SVD for the GSL ties (S6, #250), same priority as #242.
3. Main window (S1 #251, S2 #243).
4. App basics: Help, Quit, shortcuts (S3 #252) and persistence (S4 #253).
5. Settings scope (S5 #254).
6. Guided mode (S7 #256).
7. Language (S8 #257); multi-run NetCDF export (#255) any time.

## S1 Main window

- Layout:
  - **Left:** runs, and the component list below them.
  - **Centre:** the TIC with the view toolbar (S2) at the top. Below it, the spectrum and Detected Peaks side by side.
  - **Right:** component info, with the detected component and its library entry side by side (RT, RI, ions, centroid ions, S/N, sharpness against RI, RT, ions, spectrum score, RI score). Below that, the peaks and the hits.
- Every area can be resized and hidden from a **View** menu. The layout is saved (S4).
- **Component list:** always visible, and usable without a library.
  - Without a library, each row is "No match: N" with RT, RI (when calibrated), area, S/N and number of ions.
  - With a library, name, score, spectrum score and RI score are added.
  - Sortable. A click selects the component; ↑/↓ steps through the list.
- **Detected Peaks:** for the selected component, the intensity of each of its ions over the scans of the component, overlaid. A click on a curve, or on an ion in the peak table, highlights it. The modelled peak profile can be switched on.
- **Spectrum:** stays visible at small window sizes.
- **Direct entry:** the scan and the RT can be typed in.
- **Bugs:**
  - The TIC zoom range must reset when the run changes.
  - A click near a component must snap to it within a sensible pixel radius, not only within 3 scans.

## S2 Component markers and view toolbar (extends #243)

- Keep the dots on the peaks.
- Add triangles on the x axis as larger click targets.
- View toolbar above the plots:
  - zoom in, zoom out, show all
  - pan
  - crosshair
  - info
  - score filter as coloured triangles (≤ 0.7, > 0.7, > 0.9), replacing today's checkboxes

## S3 App basics

- **Help menu:**
  - About (version, MIT licence, authors, citation, source link)
  - Manual (README with screenshots)
  - Open log folder
  - Check for updates (existing)
  - Guided mode on/off
- **File:** Quit (Ctrl+Q).
- **Shortcuts,** shown in the menus:

  | Shortcut | Action |
  |---|---|
  | Ctrl+O | open |
  | Ctrl+S | save components |
  | Ctrl+E | export results CSV |
  | Ctrl+L | add component to library |
  | Ctrl+Q | quit |
  | F1 | help |
  | ←/→ | scan |
  | ↑/↓ | component |

## S4 Persistence

- Saved to a file in the user profile:
  - settings
  - window layout
  - recent files
  - compound library and NIST library paths
  - mode (S7)
  - language (S8)
- First start: the app imports the old `MetaboliteDetector.conf` (Linux, `~/.config/TU-Braunschweig/`) or the old registry keys (Windows) once, if present, and says so.
  - The old app's built-in `LS_RI_DIFF` is 20. A user's saved value (for example 10) wins.

## S5 Settings scope

- **Add:**
  - Cutoff score
  - Score mode: combined, spectrum or RI (replaces the two use-RI/use-spectrum checkboxes)
  - Use scaled library
  - Compound library and NIST library paths
  - Axes labels and labels for the highest peaks
  - Guided mode switch (also in Help)
- **Baseline adjustment and the centroid values (SF_\*):** first confirm by RE that the old detector uses them; add them only if it does.
- Check against #242 before changing the scoring.
- **Result (#254):** RE shows `SF_PEAK_*`, `SF_MAX_BASE_DEV` and `LS_USE_SCALED_LIB` are dead in 3.3 (stored only);
  `SF_CENTROID_*` steer the centroid export, `SF_FWHM` every centroiding, and baseline adjustment only type-1
  caches (algorithms.md §4b.4). The guided-mode switch is left to S7 (#256).

## S6 GSL ties

- A clean-room SVD under MIT, with no GSL code: textbook methods checked against the old binary's GSL as a
  black box (method and boundary: `docs/reverse-engineering/gsl-clean-room.md`).
- It reproduces the old rounding on precise-apex ties.
- **Goal:** default detector 246/246 and matched filter 2807/2807 identical in the golden tests. Reached in #250.

## S7 Guided mode

- A global switch "Guided / Expert".
  - The first start asks, with Guided preselected.
  - The switch is in Help and in Settings, and is saved.
- **Four step flows:** RI calibration, component detection (batch with .cmp), batch quantification and MID.
  - Each shows "Step x of n" and can be skipped.
  - Each ends in the expert window with the values set.
- The step flows reuse the expert windows' building blocks; there is no second implementation.
- In Guided mode, expert options sit under a collapsed "Advanced". In Expert mode, it is open. Every function stays reachable in both modes.
- **Result (#256):** Analysis → Step by step starts the flows; component detection is the batch with
  "Save components as .cmp" switched on. Design and files: [architecture.md](architecture.md#guided-mode).

## S8 Language

- The default follows the system language, falling back to English. Available: English and German.
- Switchable in Help and Settings. It applies immediately and is saved.
- Only UI texts are translated. CSV headers and file formats stay English so existing scripts keep working.

## Dropped

- The old NetCDF/FastFlight batch conversion to `.bin`: the app opens both formats directly.
- "Export as NetCDF" gains export of several runs at once.
