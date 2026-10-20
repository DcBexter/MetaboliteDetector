//! Vector copy of what egui drew in one screen area, written as SVG or PDF. It reads the frame's
//! own shapes (lines, rectangles, circles, text), so axes and curves match the display exactly.
//! Meshes (images, gradients) and ellipses are not exported; the plots draw none.

use std::fmt::Write as _;

use eframe::egui::{self, Color32, Pos2, Rect, Shape, Vec2, epaint};

#[derive(Clone, Debug, PartialEq)]
pub enum Item {
    Poly {
        points: Vec<Pos2>,
        closed: bool,
        fill: Color32,
        stroke: (f32, Color32),
    },
    Circle {
        center: Pos2,
        radius: f32,
        fill: Color32,
        stroke: (f32, Color32),
    },
    Text {
        pos: Pos2,
        /// Clockwise radians around `pos`.
        angle: f32,
        /// Offset of the baseline start from `pos`.
        base: Vec2,
        text: String,
        size: f32,
        color: Color32,
    },
}

/// Items in area coordinates (origin at the area's top-left), each with the clip rectangle it was
/// drawn under.
pub struct Scene {
    pub size: Vec2,
    pub background: Color32,
    pub items: Vec<(Rect, Item)>,
}

/// Shapes of `layer` that fall into `area`, back to front.
pub fn capture(ctx: &egui::Context, layer: egui::LayerId, area: Rect) -> Scene {
    let mut raw = Vec::new();
    ctx.graphics(|g| {
        if let Some(list) = g.get(layer) {
            for s in list.all_entries() {
                walk(&s.shape, s.clip_rect.intersect(area), &mut raw);
            }
        }
    });
    let shift = -area.min.to_vec2();
    let items = raw
        .into_iter()
        .filter(|(clip, item)| {
            clip.width() > 0.0 && clip.height() > 0.0 && bounds(item).intersects(*clip)
        })
        .map(|(clip, item)| (clip.translate(shift), moved(item, shift)))
        .collect();
    let v = &ctx.global_style().visuals;
    let background = if layer.order == egui::Order::Background {
        v.panel_fill
    } else {
        v.window_fill
    };
    Scene {
        size: area.size(),
        background,
        items,
    }
}

fn walk(shape: &Shape, clip: Rect, out: &mut Vec<(Rect, Item)>) {
    let stroke = |s: &epaint::Stroke| (s.width, s.color);
    let item = match shape {
        Shape::Vec(v) => return v.iter().for_each(|s| walk(s, clip, out)),
        Shape::LineSegment { points, stroke: s } => Item::Poly {
            points: points.to_vec(),
            closed: false,
            fill: Color32::TRANSPARENT,
            stroke: stroke(s),
        },
        Shape::Path(p) => Item::Poly {
            points: p.points.clone(),
            closed: p.closed,
            fill: p.fill,
            stroke: match &p.stroke.color {
                epaint::ColorMode::Solid(c) => (p.stroke.width, *c),
                _ => (0.0, Color32::TRANSPARENT),
            },
        },
        Shape::Rect(r) => Item::Poly {
            points: vec![
                r.rect.left_top(),
                r.rect.right_top(),
                r.rect.right_bottom(),
                r.rect.left_bottom(),
            ],
            closed: true,
            fill: r.fill,
            stroke: stroke(&r.stroke),
        },
        Shape::Circle(c) => Item::Circle {
            center: c.center,
            radius: c.radius,
            fill: c.fill,
            stroke: stroke(&c.stroke),
        },
        Shape::QuadraticBezier(b) => Item::Poly {
            points: b.flatten(None),
            closed: b.closed,
            fill: b.fill,
            stroke: match &b.stroke.color {
                epaint::ColorMode::Solid(c) => (b.stroke.width, *c),
                _ => (0.0, Color32::TRANSPARENT),
            },
        },
        Shape::CubicBezier(b) => Item::Poly {
            points: b.flatten(None),
            closed: b.closed,
            fill: b.fill,
            stroke: match &b.stroke.color {
                epaint::ColorMode::Solid(c) => (b.stroke.width, *c),
                _ => (0.0, Color32::TRANSPARENT),
            },
        },
        Shape::Text(t) => return text(t, clip, out),
        _ => return,
    };
    let mut item = item;
    if let Item::Poly { points, .. } = &mut item {
        points.retain(|p| p.x.is_finite() && p.y.is_finite());
        if points.len() < 2 {
            return;
        }
    }
    let visible = match &item {
        Item::Poly { fill, stroke, .. } | Item::Circle { fill, stroke, .. } => {
            *fill != Color32::TRANSPARENT || (stroke.0 > 0.0 && stroke.1 != Color32::TRANSPARENT)
        }
        Item::Text { .. } => true,
    };
    if visible {
        out.push((clip, item));
    }
}

fn text(t: &epaint::TextShape, clip: Rect, out: &mut Vec<(Rect, Item)>) {
    let Some(format) = t.galley.job.sections.first().map(|s| &s.format) else {
        return;
    };
    let color = match (t.override_text_color, format.color) {
        (Some(c), _) => c,
        (None, Color32::PLACEHOLDER) => t.fallback_color,
        (None, c) => c,
    };
    for row in &t.galley.rows {
        let text = row.row.text();
        let Some(g) = row.row.glyphs.first() else {
            continue;
        };
        if text.trim().is_empty() {
            continue;
        }
        out.push((
            clip,
            Item::Text {
                pos: t.pos,
                angle: t.angle,
                base: egui::vec2(row.pos.x + g.pos.x, row.pos.y + g.pos.y),
                text,
                size: format.font_id.size,
                color: color.gamma_multiply(t.opacity_factor),
            },
        ));
    }
}

fn rot(v: Vec2, angle: f32) -> Vec2 {
    let (s, c) = angle.sin_cos();
    egui::vec2(v.x * c - v.y * s, v.x * s + v.y * c)
}

fn bounds(item: &Item) -> Rect {
    match item {
        Item::Poly { points, stroke, .. } => Rect::from_points(points).expand(stroke.0),
        Item::Circle {
            center,
            radius,
            stroke,
            ..
        } => Rect::from_center_size(*center, Vec2::splat(2.0 * (radius + stroke.0))),
        Item::Text {
            pos,
            angle,
            base,
            text,
            size,
            ..
        } => {
            let w = text.chars().count() as f32 * size * 0.6;
            let a = *pos + rot(*base - egui::vec2(0.0, *size), *angle);
            let b = *pos + rot(*base + egui::vec2(w, *size * 0.3), *angle);
            Rect::from_two_pos(a, b).expand(*size)
        }
    }
}

fn moved(item: Item, by: Vec2) -> Item {
    match item {
        Item::Poly {
            points,
            closed,
            fill,
            stroke,
        } => Item::Poly {
            points: points.into_iter().map(|p| p + by).collect(),
            closed,
            fill,
            stroke,
        },
        Item::Circle {
            center,
            radius,
            fill,
            stroke,
        } => Item::Circle {
            center: center + by,
            radius,
            fill,
            stroke,
        },
        Item::Text {
            pos,
            angle,
            base,
            text,
            size,
            color,
        } => Item::Text {
            pos: pos + by,
            angle,
            base,
            text,
            size,
            color,
        },
    }
}

/// Output size in points: the area itself, or scaled to `width` and/or `height` (0 = free). With
/// both set the picture keeps its proportions and is centred.
pub fn fit(area: Vec2, width: f32, height: f32) -> (Vec2, f32, Vec2) {
    let scale = match (width > 0.0, height > 0.0) {
        (false, false) => 1.0,
        (true, false) => width / area.x,
        (false, true) => height / area.y,
        (true, true) => (width / area.x).min(height / area.y),
    };
    let canvas = egui::vec2(
        if width > 0.0 { width } else { area.x * scale },
        if height > 0.0 { height } else { area.y * scale },
    );
    (canvas, scale, (canvas - area * scale) / 2.0)
}

/// egui colours carry their alpha premultiplied; SVG and PDF want the plain colour and the alpha apart.
fn rgb(c: Color32) -> String {
    let [r, g, b, _] = c.to_srgba_unmultiplied();
    format!("#{r:02x}{g:02x}{b:02x}")
}

fn opacity(c: Color32) -> f32 {
    f32::from(c.a()) / 255.0
}

fn xml(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

fn fill_attr(c: Color32) -> String {
    if c == Color32::TRANSPARENT {
        "fill=\"none\"".into()
    } else {
        format!("fill=\"{}\" fill-opacity=\"{:.3}\"", rgb(c), opacity(c))
    }
}

fn stroke_attr((w, c): (f32, Color32)) -> String {
    if w <= 0.0 || c == Color32::TRANSPARENT {
        String::new()
    } else {
        format!(
            " stroke=\"{}\" stroke-opacity=\"{:.3}\" stroke-width=\"{w:.2}\" stroke-linejoin=\"round\" stroke-linecap=\"round\"",
            rgb(c),
            opacity(c)
        )
    }
}

pub fn svg(scene: &Scene, width: f32, height: f32) -> String {
    let (canvas, scale, offset) = fit(scene.size, width, height);
    let mut s = String::new();
    let _ = write!(
        s,
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{w:.1}\" height=\"{h:.1}\" viewBox=\"0 0 {w:.1} {h:.1}\">\n\
         <rect width=\"100%\" height=\"100%\" {bg}/>\n<g transform=\"translate({ox:.2},{oy:.2}) scale({scale:.5})\">\n",
        w = canvas.x,
        h = canvas.y,
        bg = fill_attr(scene.background),
        ox = offset.x,
        oy = offset.y,
    );
    let mut clip = None;
    let mut ids = 0;
    for (c, item) in &scene.items {
        if clip != Some(*c) {
            if clip.is_some() {
                s.push_str("</g>\n");
            }
            let _ = write!(
                s,
                "<clipPath id=\"c{ids}\"><rect x=\"{:.2}\" y=\"{:.2}\" width=\"{:.2}\" height=\"{:.2}\"/></clipPath>\n<g clip-path=\"url(#c{ids})\">\n",
                c.min.x,
                c.min.y,
                c.width(),
                c.height()
            );
            ids += 1;
            clip = Some(*c);
        }
        match item {
            Item::Poly {
                points,
                closed,
                fill,
                stroke,
            } => {
                let pts: Vec<String> = points
                    .iter()
                    .map(|p| format!("{:.2},{:.2}", p.x, p.y))
                    .collect();
                let tag = if *closed || *fill != Color32::TRANSPARENT {
                    "polygon"
                } else {
                    "polyline"
                };
                let _ = writeln!(
                    s,
                    "<{tag} points=\"{}\" {}{}/>",
                    pts.join(" "),
                    fill_attr(if tag == "polygon" {
                        *fill
                    } else {
                        Color32::TRANSPARENT
                    }),
                    stroke_attr(*stroke)
                );
            }
            Item::Circle {
                center,
                radius,
                fill,
                stroke,
            } => {
                let _ = writeln!(
                    s,
                    "<circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{radius:.2}\" {}{}/>",
                    center.x,
                    center.y,
                    fill_attr(*fill),
                    stroke_attr(*stroke)
                );
            }
            Item::Text {
                pos,
                angle,
                base,
                text,
                size,
                color,
            } => {
                let _ = writeln!(
                    s,
                    "<text transform=\"translate({:.2},{:.2}) rotate({:.3})\" x=\"{:.2}\" y=\"{:.2}\" font-family=\"Helvetica, Arial, sans-serif\" font-size=\"{size:.2}\" fill=\"{}\" fill-opacity=\"{:.3}\" xml:space=\"preserve\">{}</text>",
                    pos.x,
                    pos.y,
                    angle.to_degrees(),
                    base.x,
                    base.y,
                    rgb(*color),
                    opacity(*color),
                    xml(text)
                );
            }
        }
    }
    if clip.is_some() {
        s.push_str("</g>\n");
    }
    s.push_str("</g>\n</svg>\n");
    s
}

/// `c` over `bg`, as opaque RGB in 0..=1 (PDF fills here carry no alpha).
fn flat(c: Color32, bg: Color32) -> [f32; 3] {
    let a = opacity(c);
    let [r, g, b, _] = c.to_srgba_unmultiplied();
    let mix = |f: u8, back: u8| (f32::from(f) * a + f32::from(back) * (1.0 - a)) / 255.0;
    [mix(r, bg.r()), mix(g, bg.g()), mix(b, bg.b())]
}

/// WinAnsi bytes of `s` as a PDF string literal; characters outside Latin-1 become `?`.
fn pdf_text(s: &str) -> Vec<u8> {
    let mut out = vec![b'('];
    for ch in s.chars() {
        let b = u8::try_from(u32::from(ch)).unwrap_or(b'?');
        if matches!(b, b'(' | b')' | b'\\') {
            out.push(b'\\');
        }
        out.push(b);
    }
    out.push(b')');
    out
}

/// One-page PDF of the scene with the standard Helvetica font (nothing is embedded).
pub fn pdf(scene: &Scene, width: f32, height: f32) -> Vec<u8> {
    let (canvas, scale, offset) = fit(scene.size, width, height);
    let bg = scene.background;
    // Area point (x, y) to page (x', H - y'): PDF's y axis points up.
    let at = |p: Pos2| (offset.x + p.x * scale, canvas.y - (offset.y + p.y * scale));
    let mut c: Vec<u8> = Vec::new();
    let put = |c: &mut Vec<u8>, s: String| c.extend_from_slice(s.as_bytes());
    let col = |c: Color32| {
        let f = flat(c, bg);
        format!("{:.3} {:.3} {:.3}", f[0], f[1], f[2])
    };
    put(
        &mut c,
        format!("{} rg 0 0 {:.2} {:.2} re f\n", col(bg), canvas.x, canvas.y),
    );
    // Clip every item to the canvas so scaled content never spills over the page.
    for (clip, item) in &scene.items {
        let (x0, y0) = at(clip.left_bottom());
        let (x1, y1) = at(clip.right_top());
        put(
            &mut c,
            format!("q {x0:.2} {y0:.2} {:.2} {:.2} re W n\n", x1 - x0, y1 - y0),
        );
        match item {
            Item::Poly {
                points,
                closed,
                fill,
                stroke,
            } => {
                let Some((first, rest)) = points.split_first() else {
                    put(&mut c, "Q\n".into());
                    continue;
                };
                let (x, y) = at(*first);
                let mut path = format!("{x:.2} {y:.2} m ");
                for p in rest {
                    let (x, y) = at(*p);
                    let _ = write!(path, "{x:.2} {y:.2} l ");
                }
                let filled = *fill != Color32::TRANSPARENT && *closed;
                let stroked = stroke.0 > 0.0 && stroke.1 != Color32::TRANSPARENT;
                if *closed {
                    path.push_str("h ");
                }
                paint(
                    &mut c,
                    &path,
                    filled.then(|| col(*fill)),
                    stroked.then(|| (stroke.0 * scale, col(stroke.1))),
                );
            }
            Item::Circle {
                center,
                radius,
                fill,
                stroke,
            } => {
                // Four Bézier arcs.
                let k = 0.552_285 * radius * scale;
                let r = radius * scale;
                let (x, y) = at(*center);
                let path = format!(
                    "{:.2} {y:.2} m {:.2} {:.2} {:.2} {:.2} {x:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {y:.2} c {:.2} {:.2} {:.2} {:.2} {x:.2} {:.2} c {:.2} {:.2} {:.2} {:.2} {:.2} {y:.2} c h ",
                    x + r,
                    x + r,
                    y + k,
                    x + k,
                    y + r,
                    y + r,
                    x - k,
                    y + r,
                    x - r,
                    y + k,
                    x - r,
                    x - r,
                    y - k,
                    x - k,
                    y - r,
                    y - r,
                    x + k,
                    y - r,
                    x + r,
                    y - k,
                    x + r,
                );
                let stroked = stroke.0 > 0.0 && stroke.1 != Color32::TRANSPARENT;
                paint(
                    &mut c,
                    &path,
                    (*fill != Color32::TRANSPARENT).then(|| col(*fill)),
                    stroked.then(|| (stroke.0 * scale, col(stroke.1))),
                );
            }
            Item::Text {
                pos,
                angle,
                base,
                text,
                size,
                color,
            } => {
                let origin = *pos + rot(*base, *angle);
                let (x, y) = at(origin);
                // Page angle is counter-clockwise, the screen's is clockwise.
                let (s, co) = (-angle).sin_cos();
                put(
                    &mut c,
                    format!(
                        "BT {} rg /F1 {:.2} Tf {co:.5} {s:.5} {:.5} {co:.5} {x:.2} {y:.2} Tm ",
                        col(*color),
                        size * scale,
                        -s
                    ),
                );
                c.extend(pdf_text(text));
                put(&mut c, " Tj ET\n".into());
            }
        }
        put(&mut c, "Q\n".into());
    }
    let objects: [Vec<u8>; 5] = [
        b"<< /Type /Catalog /Pages 2 0 R >>".to_vec(),
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>".to_vec(),
        format!(
            "<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {:.2} {:.2}] /Contents 4 0 R /Resources << /Font << /F1 5 0 R >> >> >>",
            canvas.x, canvas.y
        )
        .into_bytes(),
        [
            format!("<< /Length {} >>\nstream\n", c.len()).into_bytes(),
            c,
            b"\nendstream".to_vec(),
        ]
        .concat(),
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica /Encoding /WinAnsiEncoding >>".to_vec(),
    ];
    let mut out = b"%PDF-1.4\n".to_vec();
    let mut offsets = Vec::new();
    for (i, o) in objects.iter().enumerate() {
        offsets.push(out.len());
        out.extend(format!("{} 0 obj\n", i + 1).into_bytes());
        out.extend(o);
        out.extend(b"\nendobj\n");
    }
    let xref = out.len();
    out.extend(format!("xref\n0 {}\n0000000000 65535 f \n", objects.len() + 1).into_bytes());
    for o in offsets {
        out.extend(format!("{o:010} 00000 n \n").into_bytes());
    }
    out.extend(
        format!(
            "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n",
            objects.len() + 1
        )
        .into_bytes(),
    );
    out
}

fn paint(c: &mut Vec<u8>, path: &str, fill: Option<String>, stroke: Option<(f32, String)>) {
    let mut s = String::new();
    if let Some(f) = &fill {
        let _ = write!(s, "{f} rg ");
    }
    if let Some((w, col)) = &stroke {
        let _ = write!(s, "{col} RG {w:.2} w 1 j 1 J ");
    }
    let op = match (fill.is_some(), stroke.is_some()) {
        (true, true) => "B",
        (true, false) => "f",
        (false, true) => "S",
        (false, false) => "n",
    };
    let _ = writeln!(s, "{path}{op}");
    c.extend(s.as_bytes());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scene() -> Scene {
        let clip = Rect::from_min_size(Pos2::ZERO, egui::vec2(200.0, 100.0));
        Scene {
            size: clip.size(),
            background: Color32::WHITE,
            items: vec![
                (
                    clip,
                    Item::Poly {
                        points: vec![
                            Pos2::new(0.0, 100.0),
                            Pos2::new(50.0, 20.5),
                            Pos2::new(200.0, 80.0),
                        ],
                        closed: false,
                        fill: Color32::TRANSPARENT,
                        stroke: (1.5, Color32::RED),
                    },
                ),
                (
                    clip,
                    Item::Text {
                        pos: Pos2::new(10.0, 10.0),
                        angle: 0.0,
                        base: egui::vec2(0.0, 10.0),
                        text: "m/z (a) <b> µ".into(),
                        size: 12.0,
                        color: Color32::BLACK,
                    },
                ),
            ],
        }
    }

    #[test]
    fn svg_keeps_points_text_and_escapes() {
        let s = svg(&scene(), 0.0, 0.0);
        assert!(s.contains("points=\"0.00,100.00 50.00,20.50 200.00,80.00\""));
        assert!(s.contains("stroke=\"#ff0000\""));
        assert!(s.contains("m/z (a) &lt;b&gt; µ"));
        assert!(s.contains("viewBox=\"0 0 200.0 100.0\""));
    }

    #[test]
    fn size_scales_uniformly_and_centres() {
        let (canvas, scale, offset) = fit(egui::vec2(200.0, 100.0), 400.0, 400.0);
        assert_eq!((canvas.x, canvas.y, scale), (400.0, 400.0, 2.0));
        assert_eq!((offset.x, offset.y), (0.0, 100.0));
        let (canvas, scale, _) = fit(egui::vec2(200.0, 100.0), 100.0, 0.0);
        assert_eq!((canvas.y, scale), (50.0, 0.5));
    }

    #[test]
    fn pdf_has_valid_structure() {
        let p = pdf(&scene(), 0.0, 0.0);
        // The micro sign is the single WinAnsi byte 0xB5, so the file is not UTF-8: read it as ASCII.
        let text: String = p
            .iter()
            .map(|b| if b.is_ascii() { char::from(*b) } else { '?' })
            .collect();
        assert!(text.starts_with("%PDF-1.4"));
        assert!(text.contains("/MediaBox [0 0 200.00 100.00]"));
        assert!(text.contains("(m/z \\(a\\) <b> ?) Tj"));
        // The cross-reference table points at every object.
        let start: usize = text
            .rsplit("startxref\n")
            .next()
            .and_then(|t| t.lines().next())
            .and_then(|n| n.parse().ok())
            .unwrap();
        assert!(text[start..].starts_with("xref"));
        let table: Vec<usize> = text[start..]
            .lines()
            .skip(3)
            .take(5)
            .map(|l| l[..10].parse().unwrap())
            .collect();
        for (i, o) in table.iter().enumerate() {
            assert!(text[*o..].starts_with(&format!("{} 0 obj", i + 1)));
        }
    }

    /// A real egui_plot frame: the exported line has the plotted points and the axes carry text.
    #[test]
    fn captures_a_plot_drawn_by_egui() {
        use egui_plot::{Line, Plot};
        let ctx = egui::Context::default();
        let input = egui::RawInput {
            screen_rect: Some(Rect::from_min_size(Pos2::ZERO, egui::vec2(400.0, 300.0))),
            ..Default::default()
        };
        let mut scene = None;
        // The second frame lays text out with loaded fonts.
        for _ in 0..2 {
            let mut out = ctx.run_ui(input.clone(), |ui| {
                let top = ui.cursor().min;
                Plot::new("p").x_axis_label("Time (min)").show(ui, |p| {
                    p.line(Line::new("t", vec![[0.0, 0.0], [1.0, 2.0], [2.0, 1.0]]));
                });
                let area =
                    Rect::from_min_max(top, egui::pos2(ui.max_rect().right(), ui.cursor().min.y));
                scene = Some(capture(ui.ctx(), ui.layer_id(), area));
            });
            out.textures_delta.clear();
        }
        let scene = scene.unwrap();
        let has_curve = scene.items.iter().any(|(_, i)| {
            matches!(i, Item::Poly { points, stroke, .. } if points.len() == 3 && stroke.0 > 0.0)
        });
        assert!(has_curve, "the three plotted points are one polyline");
        // Baselines lie inside the picture: the axis title is not pushed out of the bottom.
        let title = scene.items.iter().find_map(|(_, i)| match i {
            Item::Text {
                pos, base, text, ..
            } if text == "Time (min)" => Some(pos.y + base.y),
            _ => None,
        });
        assert!(title.is_some_and(|y| y <= scene.size.y), "{title:?}");
        let svg = svg(&scene, 0.0, 0.0);
        assert!(svg.contains("Time (min)"));
        assert!(pdf(&scene, 300.0, 0.0).starts_with(b"%PDF-1.4"));
        if let Some(dir) = std::env::var_os("MD_PLOT_DIR") {
            let dir = std::path::PathBuf::from(dir);
            std::fs::write(dir.join("plot.svg"), &svg).unwrap();
            std::fs::write(dir.join("plot.pdf"), pdf(&scene, 0.0, 0.0)).unwrap();
        }
    }

    #[test]
    fn transparent_colours_are_unpremultiplied_and_nan_points_dropped() {
        // 50 % white is stored as (128, 128, 128, 128); the file wants white at 50 %.
        let c = Color32::from_white_alpha(128);
        assert_eq!(rgb(c), "#ffffff");
        let mut out = Vec::new();
        let line = Shape::line(
            vec![
                Pos2::new(0.0, 1.0),
                Pos2::new(f32::NAN, 2.0),
                Pos2::new(5.0, 5.0),
            ],
            epaint::Stroke::new(1.0, Color32::RED),
        );
        walk(&line, Rect::EVERYTHING, &mut out);
        assert!(matches!(&out[0].1, Item::Poly { points, .. } if points.len() == 2));
    }
}
