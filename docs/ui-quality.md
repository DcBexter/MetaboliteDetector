# UI quality

How the egui app is reviewed. The design rules are in `DESIGN.md`, and users and principles are in `PRODUCT.md`.

## Seeing the app without a window

```powershell
cargo test --release -p metabolite-detector screenshots -- --ignored
```

This test needs the example data (`example-data/fetch.*`) and a GPU. It writes `target/ui-snapshots/`:
- one PNG per screen: start and main window in dark and light, 15-run overlay, TIC + SIC, both menus,
  settings, calibration, batch, isotopes, library editor with a selected compound, results, run info,
  MID settings and results, and a batch table and statistics (group 0 against group 2) from made-up
  quantities (`Batch::finished`);
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
- [ ] Labels are short whole English phrases (German comes later).
- [ ] Errors appear next to their cause and can be cleared. Background work shows a spinner.

## Re-check after the open UI pull requests merge

Second review (#168, 2026-10-20) of results, chromatogram, calibration, isotopes, statistics, library
editor, MID and window placement. Fixed in #168: Run info wraps a long path in a readable column and
names the peak settings like the Settings dialog; statistics explains why it needs two replicate groups;
the status bar shows load times as "80 ms" or "1.2 s". #183 fixed the MID findings (one rescan label, hovers for
R² and negative values). Calibration, isotopes, chromatogram and component detail had no findings, and
every screen except batch had no unnamed input.

#169, #170, #171 and #187 are done. #175 named every batch input, gave each a hover text and renamed
`QUANT_EXCL_IONS` to a plain name; it also named the library editor's feature fields. No screenshot
lists an unnamed input now. Nothing from the second review is open.

Third review (#219, 2026-10-20) of every screenshot on `main` after #175, #187, #189, #190 and #196.
Fixed in #219: the results table gives the compound name the room the fixed columns leave, and its score
headers read "Score" and "Spectrum" like the library hits; biomarker statistics shows means and standard
deviations with 6 significant digits like the batch table ("7700", not "7700.0000"); RI calibration
opens wide enough that its first checkbox fits on one line. Start screen, menus, settings, isotopes, run
info and the batch, library editor and MID windows had no findings, and no screenshot lists an unnamed
input.
