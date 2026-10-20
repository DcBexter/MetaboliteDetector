//! Spectrum normalisation and export of the shown chromatogram and spectrum as CSV or PNG.
//! CSV rows come from the same numbers the plots draw; PNG is a crop of egui's own screenshot.

use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use eframe::egui;
use md_core::csv::write_columns;
use md_core::spectrum::{Normalization, normalize};

use crate::{App, Run, Trace};

pub struct Export {
    pub norm: Normalization,
    /// Screen areas of the plots in the last frame, in points.
    pub chromatograms: egui::Rect,
    pub spectrum: egui::Rect,
    /// Errors of screenshot callbacks, which run after the frame that asked for them.
    errors: (mpsc::Sender<String>, mpsc::Receiver<String>),
}

impl Default for Export {
    fn default() -> Self {
        Self {
            norm: Normalization::default(),
            chromatograms: egui::Rect::NOTHING,
            spectrum: egui::Rect::NOTHING,
            errors: mpsc::channel(),
        }
    }
}

/// The spectrum on screen: title, m/z and (normalised) intensities. The deconvoluted spectrum of a
/// selected component replaces the raw scan.
pub fn shown_spectrum(run: &Run, norm: Normalization) -> (String, Vec<f64>, Vec<f64>) {
    let (title, mz, intensity): (_, Vec<f64>, Vec<f64>) = match run.component {
        Some(k) => {
            let s = &run.components[k].spectrum;
            (
                format!("Component {} (deconvoluted)", k + 1),
                s.iter().map(|(m, _)| f64::from(*m)).collect(),
                s.iter().map(|(_, h)| f64::from(*h)).collect(),
            )
        }
        None => {
            let (mz, i) = run.data.spectrum(run.scan);
            (
                format!("Scan {}", run.scan + 1),
                mz.iter().map(|m| f64::from(*m)).collect(),
                i.iter().map(|v| f64::from(*v)).collect(),
            )
        }
    };
    (title, mz, normalize(&intensity, norm))
}

/// Writes through `<path>.part` and renames, so a failed export never leaves a half-written file.
pub(crate) fn write_file(
    path: &Path,
    write: impl FnOnce(&mut BufWriter<std::fs::File>) -> std::io::Result<()>,
) -> Result<(), String> {
    let mut part = path.as_os_str().to_owned();
    part.push(".part");
    let part = PathBuf::from(part);
    let run = || -> std::io::Result<()> {
        let mut w = BufWriter::new(std::fs::File::create(&part)?);
        write(&mut w)?;
        w.flush()?;
        std::fs::rename(&part, path)
    };
    run().map_err(|e| {
        std::fs::remove_file(&part).ok();
        format!("{}: {e}", path.display())
    })
}

/// PNG of `rect` (points) out of a full-window screenshot taken at `ppp` pixels per point.
pub(crate) fn write_png(
    path: &Path,
    shot: &egui::ColorImage,
    rect: egui::Rect,
    ppp: f32,
) -> Result<(), String> {
    let size = egui::vec2(shot.width() as f32, shot.height() as f32) / ppp;
    let rect = rect.intersect(egui::Rect::from_min_size(egui::Pos2::ZERO, size));
    if rect.width() < 1.0 || rect.height() < 1.0 {
        return Err("nothing to export: the plot is not visible".into());
    }
    let img = shot.region(&rect, Some(ppp));
    write_file(path, |w| {
        let mut enc = png::Encoder::new(&mut *w, img.width() as u32, img.height() as u32);
        enc.set_color(png::ColorType::Rgba);
        enc.set_depth(png::BitDepth::Eight);
        let mut enc = enc.write_header()?;
        let bytes: Vec<u8> = img.pixels.iter().flat_map(|p| p.to_array()).collect();
        enc.write_image_data(&bytes)?;
        enc.finish().map_err(Into::into)
    })
}

impl App {
    /// Shown chromatogram of the active run as CSV: time and the traces the plots draw.
    fn chromatogram_csv(&self, path: &Path) -> Result<(), String> {
        let run = self.runs.get(self.active).ok_or("no run open")?;
        let mut header = vec![format!("Retention time ({})", self.view.unit.label())];
        let times: Vec<f64> = run.minutes.iter().map(|m| self.view.time(*m)).collect();
        // The plots show the derivative when it is on, so the file does too.
        let shape = |v: &[f64]| {
            if self.view.derivative {
                md_core::first_derivative(v)
            } else {
                v.to_vec()
            }
        };
        let mut owned = Vec::new();
        if self.trace != Trace::Sic {
            header.push(self.derivative_label("TIC".into()));
            owned.push(shape(&run.tic));
        }
        if self.trace != Trace::Tic {
            if run.sic_for != Some((self.mz, self.tolerance)) || run.sic.len() != run.minutes.len()
            {
                return Err("the SIC is still being computed, try again in a moment".into());
            }
            header.push(self.derivative_label(format!("SIC m/z {} ± {}", self.mz, self.tolerance)));
            owned.push(shape(&run.sic));
        }
        let mut columns: Vec<&[f64]> = vec![&times];
        columns.extend(owned.iter().map(Vec::as_slice));
        let header: Vec<&str> = header.iter().map(String::as_str).collect();
        write_file(path, |w| write_columns(w, &header, &columns))
    }

    fn spectrum_csv(&self, path: &Path) -> Result<(), String> {
        let run = self.runs.get(self.active).ok_or("no run open")?;
        let (_, mz, y) = shown_spectrum(run, self.export.norm);
        let y_name = match self.export.norm {
            Normalization::None => "Intensity".to_owned(),
            n => format!("Intensity ({})", n.label().to_lowercase()),
        };
        write_file(path, |w| write_columns(w, &["m/z", &y_name], &[&mz, &y]))
    }

    /// File → Export items for the active run; save dialogs run on the UI thread like the NetCDF export.
    pub fn export_menu(&mut self, ui: &mut egui::Ui) {
        let Some(name) = self.runs.get(self.active).map(|r| r.name.clone()) else {
            return;
        };
        let save = |what: &str, kind: &str, ext: &str| {
            rfd::FileDialog::new()
                .add_filter(kind, &[ext])
                .set_file_name(format!("{name}-{what}.{ext}"))
                .save_file()
        };
        let mut result = None;
        if ui.button("Export chromatogram as CSV…").clicked()
            && let Some(p) = save("chromatogram", "CSV", "csv")
        {
            result = Some(self.chromatogram_csv(&p));
        }
        if ui.button("Export spectrum as CSV…").clicked()
            && let Some(p) = save("spectrum", "CSV", "csv")
        {
            result = Some(self.spectrum_csv(&p));
        }
        for (label, what, rect) in [
            (
                "Export chromatogram as PNG…",
                "chromatogram",
                self.export.chromatograms,
            ),
            ("Export spectrum as PNG…", "spectrum", self.export.spectrum),
        ] {
            if ui.button(label).clicked()
                && let Some(path) = save(what, "PNG image", "png")
            {
                // The screenshot is read after the next frame, so the dialog is gone from it.
                let (ctx, errors) = (ui.ctx().clone(), self.export.errors.0.clone());
                let ppp = ctx.pixels_per_point();
                ctx.request_screenshot(move |shot| {
                    if let Err(e) = write_png(&path, &shot, rect, ppp) {
                        errors.send(e).ok();
                    }
                });
            }
        }
        if let Some(Err(e)) = result {
            self.error = Some(e);
        }
    }

    /// Shows errors of finished screenshot exports.
    pub fn poll_export(&mut self) {
        if let Ok(e) = self.export.errors.1.try_recv() {
            self.error = Some(e);
        }
    }
}

/// Normalisation picker above the spectrum plot.
pub fn norm_picker(ui: &mut egui::Ui, norm: &mut Normalization) {
    ui.horizontal_wrapped(|ui| {
        ui.label("Normalise spectrum:");
        for n in Normalization::ALL {
            ui.selectable_value(norm, n, n.label());
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn png_is_a_crop_of_the_screenshot() {
        // 2x2 points at 2 pixels per point: left half red, right half blue.
        let (r, b) = (egui::Color32::RED, egui::Color32::BLUE);
        let shot = egui::ColorImage::new([4, 4], [r, r, b, b].repeat(4));
        let dir = std::env::temp_dir().join(format!("md-export-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("crop.png");
        let right = egui::Rect::from_min_max(egui::pos2(1.0, 0.0), egui::pos2(2.0, 2.0));
        write_png(&path, &shot, right, 2.0).unwrap();
        let file = std::io::BufReader::new(std::fs::File::open(&path).unwrap());
        let mut reader = png::Decoder::new(file).read_info().unwrap();
        let mut buf = vec![0; reader.output_buffer_size().unwrap()];
        let info = reader.next_frame(&mut buf).unwrap();
        assert_eq!((info.width, info.height), (2, 4));
        assert!(buf.chunks(4).all(|p| p == [0, 0, 255, 255]));
        assert!(write_png(&path, &shot, egui::Rect::NOTHING, 2.0).is_err());
        std::fs::remove_dir_all(dir).ok();
    }
}
