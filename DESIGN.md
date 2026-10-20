---
name: MetaboliteDetector
description: Native GC-MS metabolomics workbench — chromatograms, spectra, identification and quantification.
colors:
  panel-dark: "#1b1b1b"
  plot-well-dark: "#0a0a0a"
  stroke-dark: "#3c3c3c"
  panel-light: "#f8f8f8"
  plot-well-light: "#ffffff"
  selection-dark: "#005c80"
  selection-light: "#90d1ff"
  warn-dark: "#ff8f00"
  warn-light: "#ff6400"
  error: "#ff0000"
  score-high: "#009600"
  score-mid: "#dc8c00"
  score-low: "#d22828"
typography:
  body:
    fontFamily: "Ubuntu Light, Hack (fallback for symbols)"
    fontSize: "13px"
    fontWeight: 300
  button:
    fontFamily: "Ubuntu Light"
    fontSize: "13px"
    fontWeight: 300
  heading:
    fontFamily: "Ubuntu Light"
    fontSize: "18px"
    fontWeight: 300
  small:
    fontFamily: "Ubuntu Light"
    fontSize: "9px"
    fontWeight: 300
  mono:
    fontFamily: "Hack"
    fontSize: "13px"
    fontWeight: 400
rounded:
  static: "4px"
  control: "6px"
  window: "6px"
spacing:
  item-x: "8px"
  item-y: "3px"
  window-margin: "10px"
  button-x: "8px"
  button-y: "2px"
  group-gap: "8px"
components:
  button:
    textColor: "#d7d7d7"
    rounded: "{rounded.control}"
    padding: "2px 8px"
  toggle-selected:
    backgroundColor: "{colors.selection-dark}"
    rounded: "{rounded.control}"
  window:
    backgroundColor: "{colors.panel-dark}"
    rounded: "{rounded.window}"
    padding: "6px"
---

# Design System: MetaboliteDetector

## Overview

**Creative North Star: "The Instrument Panel"**

The app sits next to a GC-MS instrument and should feel like the best software that ships with one:
dense, calm and exact. The data — chromatograms, spectra, score tables — owns the space; the chrome
around it is quiet grey, small and consistent, and never competes with a peak. Colour is reserved for
data (one hue per run, score colours from 3.3) and for state (selection, warning, error).

The visual base is stock egui (0.37) on purpose: native-feeling, fast, and themed by the OS (dark and
light are both first-class, every change is checked in both). The app adds only what egui lacks:
a symbol font fallback (`crates/md-app/src/style.rs`) and shared layout rules. New screens extend this
system; they do not introduce their own colours, fonts or corner radii.

**Key Characteristics:**
- Plots fill the centre; side panels (runs left, hits right) and a one-line status bar frame them.
- Tools open as floating, resizable windows over the plots (settings, calibration, batch, editor).
  The first time, each opens right of the run list and a step below and right of the one before,
  with its starting size capped to the room there (`layout::ToolWindow`). This holds down to 800 px
  width; in a narrower window, content wider than the room still covers the list. egui keeps where
  the user moves it.
- Every control carries a short hover explanation; settings show plain names with the 3.3 key.
- Light and dark themes follow the OS.

## Colors

Neutral greys from egui's visuals for all chrome; saturated colour only where it carries meaning.

### Primary
- **Selection Blue** (dark `{colors.selection-dark}`, light `{colors.selection-light}`): the selected
  toggle (TIC/SIC, normalisation), the active run, the selected list row. egui's own selection colour.

### Data
- **Run hues**: run *i* gets `Hsva(i × 0.618034, 0.85, 0.5)` — the golden-ratio hue walk egui_plot uses —
  identical in the run list swatch, the TIC, the SIC and any export (`color()` in `main.rs`).
- **Score High / Mid / Low** (`{colors.score-high}`, `{colors.score-mid}`, `{colors.score-low}`):
  library hit scores above 0.9, above 0.7, the rest, as in 3.3. Always printed next to the number.

### Neutral
- **Panel** (dark `{colors.panel-dark}`, light `{colors.panel-light}`): panels and windows.
- **Plot Well** (dark `{colors.plot-well-dark}`, light `{colors.plot-well-light}`): plot backgrounds and
  text fields (egui `extreme_bg_color`), one step darker/lighter than the panel.
- **Stroke** (`{colors.stroke-dark}`): separators, window borders.

### Named Rules
**The Data Owns Colour Rule.** Chrome is grey. A hue on screen is a run, a score or a state — never
decoration.

**The Two Signals Rule.** Warnings use `warn_fg_color`, errors `error_fg_color` from the active visuals,
inline next to the thing that caused them. No other alert colours.

## Typography

**Body Font:** Ubuntu Light (egui default), with Hack as a fallback for symbols (`→ ← ■`).
**Mono Font:** Hack, for numbers in code-like places only.

**Character:** One light humanist sans at one size. Hierarchy comes from `strong` (brighter text) and
from grouping, not from size jumps.

### Hierarchy
- **Heading** (18px): window titles only (egui draws them).
- **Body / Button** (13px): everything else, including table cells and plot axes.
- **Strong**: section titles inside windows ("Peak detection and deconvolution"), table headers,
  the assigned hit.
- **Small** (9px): not used; too small for lab monitors.

### Named Rules
**The No Tofu Rule.** Every non-ASCII symbol in a UI string must have a glyph in the bundled fonts;
`every_ui_symbol_has_a_glyph` in `ui_tests.rs` enforces it.

**The Translatable Label Rule.** Labels are whole, short English phrases; a German translation is
planned, so no sentence is assembled from fragments.

## Layout

- Top: menu bar (`File`, `Analysis`). Bottom: one status line. Left: run list (when runs are open).
  Right: hit list of the selected component (when a library is loaded). Centre: trace toggles, TIC
  and/or SIC plot(s), normalisation toggles, spectrum plot.
- Chromatogram plots share the x axis (linked) and a fixed 90 px y-axis width so axes line up.
- Spacing is egui's: 8 px between items in a row, 3 px between rows, 8 px (`add_space(8.0)`) between
  groups inside a window. Window margin 10 px, the corner resize grip, so no control sits under it;
  scroll bars are solid and take their own room (`style::apply`).
- Forms are two-column `egui::Grid`s: label left, control right. Tables are striped grids.
- Default window size 1200 × 800; nothing may overlap or be cut off at the side down to
  800 × 600 and at 520 px width (overlap test, #89). Rows of controls use `horizontal_wrapped`;
  the run list and the component panel grow to their content up to 30 % of the width, then scroll
  (`layout::side_panel`); wide tables scroll sideways.

## Elevation & Depth

Flat. Panels sit edge to edge; floating windows get egui's soft shadow and 1 px stroke. No other
shadows or layering.

## Shapes

Small, quiet corners: 4 px on static widgets, 6 px on controls, windows and menus. Plot markers are filled circles
(component apexes, radius 3). The run swatch is a filled square.

## Components

### Buttons
- **Shape:** 6 px corners, 8 × 2 px padding, egui's grey fill.
- **States:** hover brightens; disabled is dimmed and should say why on hover
  (`on_disabled_hover_text`).
- **Labels:** verbs or nouns of the action; `…` when a dialog follows.

### Toggles (selectable values)
- TIC / SIC / TIC + SIC and None / Base peak / Unit / Relative: text-only, selected one filled with
  Selection Blue.

### Inputs
- DragValues for numbers (with sensible ranges), single-line text edits with hint text showing an
  example (`e.g. 174 248`).

### Lists and tables
- Run list rows: visibility checkbox, colour swatch, selectable name, close button; full path on hover.
- Tables: striped `Grid`, `strong` header row, numbers with fixed decimals.

### Plots (signature component)
- egui_plot with axis titles in sentence case ("Retention time (min)", "m/z").
- Lines decimated to ~2 points per pixel; one colour per run; a vertical line marks the scan cursor.
- Hit labels on the TIC in score colours.

## Do's and Don'ts

### Do:
- **Do** give every button, toggle and setting a hover explanation of what it does.
- **Do** show settings by plain name, with the old 3.3 key on hover.
- **Do** take colours from `ui.visuals()` (selection, warn, error) or the run/score palettes above.
- **Do** check every change in dark and light (`cargo test --release -p metabolite-detector screenshots -- --ignored`).
- **Do** show background work in the status bar with a spinner, and errors inline where they happen.

### Don't:
- **Don't** add custom colours, fonts, corner radii or shadows outside this file.
- **Don't** use a symbol without a glyph in the bundled fonts (boxes in the UI).
- **Don't** shrink text below 13 px or use colour as the only carrier of a score or state.
- **Don't** let a control cover another one at any window size.
