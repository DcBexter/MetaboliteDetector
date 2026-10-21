# UI quality

How the egui app is reviewed. The design rules are in `DESIGN.md`, and users and principles are in `PRODUCT.md`.

## Seeing the app without a window

```powershell
cargo test --release -p metabolite-detector screenshots -- --ignored
```

This test needs the example data (`example-data/fetch.*`) and a GPU. It writes `target/ui-snapshots/`:
- one PNG per screen: start and main window (with and without a library) in dark and light, the main
  window at 800 × 600, 15-run overlay, TIC + SIC, both menus,
  settings, calibration, batch, isotopes, library editor with a selected compound, run info,
  MID settings and results, and a batch table and statistics (group 0 against group 2) from made-up
  quantities (`Batch::finished`); the first-start mode question in dark and light, and in Guided mode
  (dark and light) the Settings window with "Advanced" closed and one step of each of the four step flows;
- one `*.unnamed.txt` per PNG, listing the inputs a screen reader cannot name.

Pixels differ between GPUs and operating systems, so the images are never committed, except the copies the
README shows in `docs/images/` (dark set; recopy them when a shown screen changes). Compare them by eye
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
- [ ] Labels are short whole phrases with a German pair (`lang::t`, `lang::tf!`); the `*.de.png`
      screenshots show them fitting.
- [ ] Errors appear next to their cause and can be cleared. Background work shows a spinner.
