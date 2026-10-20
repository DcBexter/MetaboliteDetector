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

Findings from the first review (#102) that are still open. #111, #112, #77/#78 and #113 were re-checked
on 2026-10-20 and are done; #127 (windows beside the run list, `×`) is done too.

| Area | Issue | What to check |
|---|---|---|
| Library editor | #80 / #84 | No "RT 0.00" or "RI known 0" for compounds without RT/RI. `11-library-editor` now selects the first compound of the editor list |
| Batch | #88 | A disabled Start button looks like a button and says why; inputs have names (8 unnamed) |
