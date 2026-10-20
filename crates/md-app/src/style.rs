//! App-wide look that egui's defaults do not give us.

use eframe::egui;

/// egui's Ubuntu font has no arrows or squares (`→`, `←`, `■`); the bundled Hack font has them, so it
/// backs up the proportional text instead of drawing boxes.
pub fn apply(ctx: &egui::Context) {
    let mut fonts = egui::FontDefinitions::default();
    if let Some(family) = fonts.families.get_mut(&egui::FontFamily::Proportional) {
        family.push("Hack".into());
    }
    ctx.set_fonts(fonts);
}
