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

Findings from the first review (#102) that live in files other sessions are changing:

| Area | Issue | What to check |
|---|---|---|
| Windows (layout) | #89, done in #127 | ✓ Windows open over the plots beside the run list, a step apart; the library editor uses `×` and `KNOWN` is gone |
| TIC plot | #111 | Hit labels do not overlap; y axis is readable |
| Spectrum, hit table | #77 / #78 | The spectrum says whether it shows a scan or a deconvoluted component; table headers are `strong` |
| Library editor | #80 / #84 | No "RT 0.00" or "RI known 0" for compounds without RT/RI; inputs have names |
| Settings, Analysis menu | #113 | Plain names, old key and explanation on hover; inputs have names |
| Batch | #88 | A disabled Start button looks like a button and says why; inputs have names |
