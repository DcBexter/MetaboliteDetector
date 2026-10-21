//! Save a plot as PNG, SVG or PDF, or copy it as an image. One menu for every plot: the caller
//! hands over the screen area of the plot and its layer; the file is written after the next frame
//! (PNG, image: egui's screenshot) or at the end of the current one (SVG, PDF: the frame's shapes).

use std::path::PathBuf;

use eframe::egui;

use crate::export::{crop, not_visible, write_file, write_png};
use crate::lang::t;
use crate::vector;

/// Where a plot was drawn: the screen rectangle (points) and the layer that holds its shapes.
#[derive(Clone, Copy)]
pub struct Area {
    pub rect: egui::Rect,
    pub layer: egui::LayerId,
}

impl Default for Area {
    fn default() -> Self {
        Self {
            rect: egui::Rect::NOTHING,
            layer: egui::LayerId::background(),
        }
    }
}

impl Area {
    /// Everything the ui drew since the cursor stood at `top`, over its full width: a plot with its axes.
    pub fn since(ui: &egui::Ui, top: egui::Pos2) -> Self {
        Self::of(
            ui,
            egui::Rect::from_min_max(top, egui::pos2(ui.max_rect().right(), ui.cursor().min.y)),
        )
    }

    pub fn of(ui: &egui::Ui, rect: egui::Rect) -> Self {
        Self {
            rect,
            layer: ui.layer_id(),
        }
    }
}

#[derive(Clone, Copy, PartialEq, Debug)]
enum Kind {
    Png,
    Svg,
    Pdf,
    Image,
}

#[derive(Clone)]
struct Pending {
    area: Area,
    kind: Kind,
    path: Option<PathBuf>,
    /// Requested width and height of SVG and PDF; 0 keeps the shown size.
    size: [f32; 2],
}

fn pending_id() -> egui::Id {
    egui::Id::NULL.with("plot export pending")
}

fn error_id() -> egui::Id {
    egui::Id::NULL.with("plot export error")
}

fn size_id() -> egui::Id {
    egui::Id::NULL.with("plot export size")
}

/// A "Save plot" menu button for a plot's own window.
pub fn button(ui: &mut egui::Ui, area: Area, stem: &str) {
    ui.menu_button(t("Save plot…", "Grafik speichern…"), |ui| {
        entries(ui, area, stem);
    });
}

/// The menu entries: pictures of the plot and their size.
pub fn entries(ui: &mut egui::Ui, area: Area, stem: &str) {
    if !area.rect.is_positive() {
        ui.weak(t(
            "Show the plot to save it as a picture (View menu).",
            "Die Grafik einblenden, um sie als Bild zu speichern (Menü Ansicht).",
        ));
        return;
    }
    let ctx = ui.ctx().clone();
    let mut size: [f32; 2] = ctx.data(|d| d.get_temp(size_id())).unwrap_or_default();
    for (label, hint, kind, ext, name) in [
        (
            t("As PNG…", "Als PNG…"),
            t(
                "Save the plot as it is shown on screen",
                "Die Grafik so speichern, wie sie am Bildschirm aussieht",
            ),
            Kind::Png,
            "png",
            t("PNG image", "PNG-Bild"),
        ),
        (
            t("As SVG…", "Als SVG…"),
            t(
                "Vector graphic; stays sharp at any size",
                "Vektorgrafik; bleibt in jeder Größe scharf",
            ),
            Kind::Svg,
            "svg",
            t("SVG image", "SVG-Bild"),
        ),
        (
            t("As PDF…", "Als PDF…"),
            t(
                "One-page vector graphic for reports",
                "Einseitige Vektorgrafik für Berichte",
            ),
            Kind::Pdf,
            "pdf",
            t("PDF document", "PDF-Dokument"),
        ),
    ] {
        if ui.button(label).on_hover_text(hint).clicked()
            && let Some(path) = rfd::FileDialog::new()
                .add_filter(name, &[ext])
                .set_file_name(format!("{stem}.{ext}"))
                .save_file()
        {
            ctx.data_mut(|d| {
                d.insert_temp(
                    pending_id(),
                    Some(Pending {
                        area,
                        kind,
                        path: Some(path),
                        size,
                    }),
                );
            });
            ctx.request_repaint();
            ui.close();
        }
    }
    if ui
        .button(t("Copy image", "Bild kopieren"))
        .on_hover_text(t(
            "Copy the plot as a picture",
            "Die Grafik als Bild kopieren",
        ))
        .clicked()
    {
        ctx.data_mut(|d| {
            d.insert_temp(
                pending_id(),
                Some(Pending {
                    area,
                    kind: Kind::Image,
                    path: None,
                    size,
                }),
            );
        });
        ctx.request_repaint();
        ui.close();
    }
    ui.separator();
    ui.label(t("Size of SVG and PDF", "Größe von SVG und PDF"));
    ui.horizontal(|ui| {
        let mut changed = false;
        let names = [t("Width", "Breite"), t("Height", "Höhe")];
        for (v, what) in size.iter_mut().zip(names) {
            let label = ui.label(what);
            changed |= ui
                .add(
                    egui::DragValue::new(v)
                        .range(0.0..=10_000.0)
                        .speed(10.0)
                        .custom_formatter(|n, _| {
                            if n == 0.0 {
                                t("as shown", "wie gezeigt").to_owned()
                            } else {
                                format!("{n:.0} pt")
                            }
                        }),
                )
                .labelled_by(label.id)
                .on_hover_text(t(
                    "0 keeps the shown size; with both set the plot keeps its proportions",
                    "0 behält die gezeigte Größe; sind beide gesetzt, bleiben die Proportionen erhalten",
                ))
                .changed();
        }
        if changed {
            ctx.data_mut(|d| d.insert_temp(size_id(), size));
        }
    });
}

/// Runs the export asked for in the last frame; returns the error of a finished one.
pub fn finish(ctx: &egui::Context) -> Option<String> {
    if let Some(p) = ctx
        .data_mut(|d| d.remove_temp::<Option<Pending>>(pending_id()))
        .flatten()
    {
        match run(ctx, p) {
            Ok(()) => {}
            Err(e) => return Some(e),
        }
    }
    ctx.data_mut(|d| d.remove_temp::<String>(error_id()))
}

fn run(ctx: &egui::Context, p: Pending) -> Result<(), String> {
    let rect = p.area.rect;
    match p.kind {
        Kind::Svg | Kind::Pdf => {
            let scene = vector::capture(ctx, p.area.layer, rect);
            if scene.items.is_empty() {
                return Err(not_visible().into());
            }
            let path = p.path.ok_or(t("no file chosen", "keine Datei gewählt"))?;
            let [w, h] = p.size;
            if p.kind == Kind::Svg {
                let svg = vector::svg(&scene, w, h);
                write_file(&path, |out| std::io::Write::write_all(out, svg.as_bytes()))
            } else {
                let pdf = vector::pdf(&scene, w, h);
                write_file(&path, |out| std::io::Write::write_all(out, &pdf))
            }
        }
        Kind::Png | Kind::Image => {
            // The screenshot is read after the next frame, so the file dialog is gone from it.
            let ppp = ctx.pixels_per_point();
            let done = ctx.clone();
            ctx.request_screenshot(move |shot| {
                let result = match &p.path {
                    Some(path) => write_png(path, &shot, rect, ppp),
                    None => crop(&shot, rect, ppp).map(|img| done.copy_image(img)),
                };
                if let Err(e) = result {
                    done.data_mut(|d| d.insert_temp(error_id(), e));
                    done.request_repaint();
                }
            });
            Ok(())
        }
    }
}
