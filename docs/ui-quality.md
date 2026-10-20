# UI quality

How the egui app is reviewed. The design rules are in `DESIGN.md`, and users and principles are in `PRODUCT.md`.

## Seeing the app without a window

```powershell
cargo test --release -p metabolite-detector screenshots -- --ignored
```

This test needs the example data (`example-data/fetch.*`) and a GPU. It writes `target/ui-snapshots/`:
- one PNG per screen: start and main window in dark and light, 15-run overlay, TIC + SIC, both menus,
  settings, calibration, batch, isotopes and library editor;
- one `*.unnamed.txt` per PNG, listing the inputs a screen reader cannot name.

Pixels differ between GPUs and operating systems, so the images are never committed. Compare them by eye
before and after a change.

`every_ui_symbol_has_a_glyph` runs in CI. It fails when a UI string uses a symbol that no bundled font
can draw.

## Checklist for every UI change

- [ ] Check it in dark and light (both screenshot sets).
- [ ] Nothing overlaps at 800 × 600 (overlap test, #89).
- [ ] Every new button, toggle and setting has a hover text that says what it does.
- [ ] Disabled controls say why on hover (`on_disabled_hover_text`).
- [ ] New inputs have a name (`labelled_by` or a label inside the widget); the `unnamed.txt` list does
      not grow.
- [ ] Colours come from `ui.visuals()` or the run and score palettes; no new colours, fonts or radii.
- [ ] Labels are short whole English phrases (German comes later).
- [ ] Errors appear next to their cause and can be cleared. Background work shows a spinner.

## Re-check after the open UI pull requests merge

Second review (#168, 2026-10-20) of results, chromatogram, calibration, isotopes, statistics, library
editor, MID and window placement. Fixed in #168: Run info wraps a long path in a readable column and
names the peak settings like the Settings dialog; statistics explains why it needs two replicate groups;
the status bar shows load times in whole ms. Calibration, isotopes, chromatogram and component detail
had no findings, and every screen except batch has no unnamed input.

| Area | Issue | What to check |
|---|---|---|
| Results list | #169 | Columns are not cut off; headers in sentence case |
| Library editor | #170 | No "RT 0.00" or "RI known 0" for compounds without RT/RI (needs a decision on how "unknown" is stored) |
| Screenshot test | #171 | The editor shot selects a compound; Run info, Results and Statistics have shots |
| MID | #141 | "Rescan missing ions" is named once; column headers without `.cdf`; R² and negative values explained on hover |
| Batch | #175 | Every input has a name and a hover text (8 unnamed) |
