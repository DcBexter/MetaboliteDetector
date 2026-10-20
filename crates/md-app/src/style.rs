//! App-wide look that egui's defaults do not give us.

use eframe::egui;

/// egui's Ubuntu font has no arrows or squares (`→`, `←`, `■`); the bundled Hack font has them, so it
/// backs up the proportional text instead of drawing boxes. Scroll bars take their own room instead
/// of floating over the controls, and window content starts past the corner resize grips.
pub fn apply(ctx: &egui::Context) {
    ctx.all_styles_mut(|s| {
        s.spacing.scroll = egui::style::ScrollStyle::solid();
        s.spacing.window_margin =
            egui::Margin::same(s.interaction.resize_grab_radius_corner.ceil() as i8);
    });
    let mut fonts = egui::FontDefinitions::default();
    if let Some(family) = fonts.families.get_mut(&egui::FontFamily::Proportional) {
        family.push("Hack".into());
    }
    ctx.set_fonts(fonts);
}
