//! Component list of the active run, below the runs (old "Analysis Results" dock with
//! `DataTableModel`, `CompoundExporter` and Find Compound): one row per component, also without a
//! library, sortable; a click or ↑/↓ selects a component, a substance can be looked up in the library
//! and the list exported as CSV.

use eframe::egui;
use egui_extras::{Column, TableBuilder};
use md_core::Chromatogram;
use md_core::results::{ResultRow, find_substance, result_rows, write_csv};
use md_core::ri::RiCalibration;
use md_core::search::Identification;
use std::sync::{Arc, Weak};

use crate::Run;
use crate::identify::Library;
use crate::lang::{Msg, msg, t, tf};

/// A column of the list in the app's spelling; the CSV keeps the old program's
/// (`md_core::results::write_csv`).
#[derive(Clone, Copy, Debug, PartialEq)]
enum Col {
    Name,
    Rt,
    Ri,
    Area,
    Score,
    Spectrum,
    RiScore,
    Sn,
    Ions,
    QuantIons,
}

/// Without a library there are no scores and no quant ions.
const PLAIN: [Col; 6] = [Col::Name, Col::Rt, Col::Ri, Col::Area, Col::Sn, Col::Ions];
const WITH_LIBRARY: [Col; 10] = [
    Col::Name,
    Col::Rt,
    Col::Ri,
    Col::Area,
    Col::Score,
    Col::Spectrum,
    Col::RiScore,
    Col::Sn,
    Col::Ions,
    Col::QuantIons,
];

impl Col {
    fn title(self) -> &'static str {
        match self {
            Col::Name => t("Compound name", "Substanzname"),
            Col::Rt => "RT (min)",
            Col::Ri => "RI",
            Col::Area => t("Peak area", "Peakfläche"),
            Col::Score => "Score",
            Col::Spectrum => t("Spectrum", "Spektrum"),
            Col::RiScore => t("RI score", "RI-Score"),
            Col::Sn => "S/N",
            Col::Ions => t("Ions", "Ionen"),
            Col::QuantIons => t("Quant ions", "Quant-Ionen"),
        }
    }

    fn column(self) -> Column {
        match self {
            // The name takes the room the fixed columns leave; a longer one is cut (the row's
            // tooltip has it whole) so the scores stay in view. Not resizable: egui keeps a resizable
            // column's first width, so it would stop following the panel.
            Col::Name => Column::remainder()
                .at_least(106.0)
                .clip(true)
                .resizable(false),
            // Each with room for the sort arrow (a column never shrinks again); right-aligned cells
            // would make automatic columns grow to the panel.
            Col::Rt => Column::initial(66.0).at_least(40.0),
            Col::Ri | Col::Sn | Col::Ions => Column::initial(44.0).at_least(40.0),
            Col::Area | Col::Spectrum => Column::initial(80.0).at_least(40.0),
            Col::Score => Column::initial(60.0).at_least(40.0),
            Col::RiScore => Column::initial(70.0).at_least(40.0),
            Col::QuantIons => Column::initial(90.0).at_least(60.0).clip(true),
        }
    }

    /// Numbers are right-aligned so digits line up; name and quant ions stay left.
    fn align(self) -> egui::Layout {
        match matches!(self, Col::Name | Col::QuantIons) {
            false => egui::Layout::right_to_left(egui::Align::Center),
            true => egui::Layout::left_to_right(egui::Align::Center),
        }
    }
}

/// A result row and what the list shows beyond the CSV.
struct Row {
    res: ResultRow,
    /// Also without a match: then the three largest peaks (`md_core::results::area`).
    area: f64,
    /// Mean S/N of the ion peaks; `None` for a component without peaks (saved by the old program).
    sn: Option<f64>,
    ions: usize,
}

/// What the cached rows were built from: the run's data, its `rev`, the calibration used without
/// a library and the cutoff. The data is held weakly, so a closed run never matches a new one at the
/// same address.
type Stamp = (Weak<Chromatogram>, u64, Option<RiCalibration>, f64);

fn same(a: &Stamp, b: &Stamp) -> bool {
    a.0.ptr_eq(&b.0) && a.1 == b.1 && a.2 == b.2 && a.3 == b.3
}

#[derive(Default)]
pub struct Results {
    /// Only rows whose name contains this (case-insensitive), and the library names offered for a lookup.
    filter: String,
    /// Column and descending.
    sort: Option<(Col, bool)>,
    /// Library compound whose components are listed.
    looked_up: Option<usize>,
    note: Msg,
    cache: Option<(Stamp, Vec<Row>)>,
}

impl Results {
    /// The list of the active run; `cal` (the run's calibration) gives the RI while no library
    /// identified the run. Hits below `cutoff` (`LS_CUTOFF`) are listed as no match.
    pub fn panel(
        &mut self,
        ui: &mut egui::Ui,
        run: Option<&mut Run>,
        cal: Option<&RiCalibration>,
        cutoff: f64,
    ) {
        ui.strong(t("Component list", "Komponentenliste"));
        let Some(run) = run else {
            return;
        };
        let lib = (run.ident.as_ref()).map(|((lib, ..), ..)| lib.clone());
        let cal = cal.filter(|_| lib.is_none());
        let stamp = (Arc::downgrade(&run.data), run.rev, cal.cloned(), cutoff);
        let mut cache = self.cache.take().filter(|c| same(&c.0, &stamp));
        if cache.is_none() {
            cache = Some((stamp, build_rows(run, cal, cutoff)));
            self.looked_up = None;
        }
        let (stamp, rows) = cache.expect("filled above");

        let mut loaded = false;
        ui.horizontal_wrapped(|ui| {
            if ui
                .button(t("Export CSV…", "CSV exportieren…"))
                .on_hover_text(t(
                    "Save the list as CSV, in the columns of the old program",
                    "Speichert die Liste als CSV, in den Spalten des alten Programms",
                ))
                .clicked()
            {
                let res: Vec<ResultRow> = rows.iter().map(|r| r.res.clone()).collect();
                self.note = export(&run.name, &res);
            }
            loaded = self.saved_buttons(ui, run, cal);
        });
        // The list below was built from the components that were just replaced.
        if loaded {
            ui.ctx().request_repaint();
            return;
        }
        ui.horizontal_wrapped(|ui| {
            let label = ui.label(t("Find substance", "Substanz suchen"));
            ui.add(
                egui::TextEdit::singleline(&mut self.filter).hint_text(match lib {
                    Some(_) => t(
                        "Name in the list and the library",
                        "Name in der Liste und der Bibliothek",
                    ),
                    None => t("Name in the list", "Name in der Liste"),
                }),
            )
            .labelled_by(label.id);
            if !self.filter.is_empty()
                && (ui.button(t("Clear", "Leeren")))
                    .on_hover_text(t(
                        "Clear the search and show every component",
                        "Leert die Suche und zeigt alle Komponenten",
                    ))
                    .clicked()
            {
                self.filter.clear();
                self.looked_up = None;
            }
        });
        let want = self.filter.to_lowercase();
        let mut selected = None;
        if let Some(lib) = &lib {
            if !want.is_empty() {
                self.lookup(ui, lib, &want);
            }
            selected = self.found_list(ui, run, lib, cutoff);
        }
        if !self.note.is_empty() {
            ui.label(&*self.note);
        }

        let mut order: Vec<usize> = (0..rows.len())
            .filter(|i| want.is_empty() || shown_name(&rows[*i].res).to_lowercase().contains(&want))
            .collect();
        if let Some((col, desc)) = self.sort {
            order.sort_by(|a, b| {
                let o = compare(&rows[*a], &rows[*b], col);
                if desc { o.reverse() } else { o }
            });
        }
        ui.horizontal_wrapped(|ui| {
            let (s, n) = (order.len(), rows.len());
            ui.label(tf!("{s} of {n} components", "{s} von {n} Komponenten"));
            if lib.is_none() {
                ui.weak(t(
                    "Open a library to name them.",
                    "Eine Bibliothek öffnen, um sie zu benennen.",
                ));
            }
        });
        if order.is_empty() && !rows.is_empty() {
            ui.weak(t(
                "No component matches. Clear the search to see all.",
                "Keine Komponente passt. Die Suche leeren, um alle zu sehen.",
            ));
        }
        // ↑/↓ step through the list as it is shown, unless a text field takes the keys.
        let step = ui.input(|i| {
            i32::from(i.key_pressed(egui::Key::ArrowDown))
                - i32::from(i.key_pressed(egui::Key::ArrowUp))
        });
        let at =
            (run.component).and_then(|k| order.iter().position(|i| rows[*i].res.component == k));
        let to = (step != 0 && !ui.ctx().text_edit_focused())
            .then(|| stepped(order.len(), at, step))
            .flatten();
        let clicked = table(ui, &mut self.sort, &rows, &order, run, lib.as_deref(), to)
            .or(selected)
            .or(to.map(|p| rows[order[p]].res.component));
        if let Some(k) = clicked {
            run.component = Some(k);
            run.scan = run.components[k].apex_scan as usize;
        }
        self.cache = Some((stamp, rows));
    }

    /// Two buttons: save the components of the run as a `.cmp` file the old program reads, or load such
    /// a file in place of them. Returns whether components were loaded.
    fn saved_buttons(
        &mut self,
        ui: &mut egui::Ui,
        run: &mut Run,
        cal: Option<&RiCalibration>,
    ) -> bool {
        if ui
            .button(t("Save components…", "Komponenten speichern…"))
            .on_hover_text(t(
                "Writes the components and their peaks as a .cmp file (old format)",
                "Schreibt die Komponenten und ihre Peaks als .cmp-Datei (altes Format)",
            ))
            .clicked()
        {
            self.save(run, cal);
        }
        let busy = run.analysis.busy();
        if ui
            .add_enabled(
                !busy,
                egui::Button::new(t("Load components…", "Komponenten laden…")),
            )
            .on_hover_text(t(
                "Replaces the components with those of a .cmp file of this run",
                "Ersetzt die Komponenten durch die einer .cmp-Datei dieses Laufs",
            ))
            .on_disabled_hover_text(t(
                "Wait until the analysis has finished",
                "Warten, bis die Analyse fertig ist",
            ))
            .clicked()
            && let Some(p) = dialog(run).pick_file()
        {
            let result = load_components(run, &p);
            self.note = match &result {
                Ok(n) => {
                    let p = p.display();
                    msg!(
                        "Loaded {n} components from {p}",
                        "{n} Komponenten aus {p} geladen"
                    )
                }
                Err(e) => msg!("Load failed: {e}", "Laden fehlgeschlagen: {e}"),
            };
            return result.is_ok();
        }
        false
    }

    pub fn note(&self) -> &str {
        &self.note
    }

    pub fn set_note(&mut self, note: Msg) {
        self.note = note;
    }

    /// Asks for a `.cmp` file and writes the components of `run` there (button and Ctrl+S).
    pub fn save(&mut self, run: &Run, cal: Option<&RiCalibration>) {
        if let Some(p) = dialog(run).save_file() {
            self.note = match save_components(run, &p, cal) {
                Ok(n) => {
                    let p = p.display();
                    msg!(
                        "Saved {n} components to {p}",
                        "{n} Komponenten in {p} gespeichert"
                    )
                }
                Err(e) => msg!("Save failed: {e}", "Speichern fehlgeschlagen: {e}"),
            };
        } else {
            self.note = Msg::default();
        }
    }

    /// Ctrl+E: the whole list of `run` as CSV, unfiltered, whether or not the list is on screen.
    pub fn export_csv(&mut self, run: &Run, cal: Option<&RiCalibration>, cutoff: f64) {
        let rows: Vec<ResultRow> = build_rows(run, cal.filter(|_| run.ident.is_none()), cutoff)
            .into_iter()
            .map(|r| r.res)
            .collect();
        self.note = export(&run.name, &rows);
    }

    /// Library names that contain the search text; a click lists the components that have it as a hit.
    fn lookup(&mut self, ui: &mut egui::Ui, lib: &Library, want: &str) {
        let names = lib
            .compounds
            .iter()
            .enumerate()
            .filter(|(_, c)| c.name.to_lowercase().contains(want))
            .take(8);
        ui.horizontal_wrapped(|ui| {
            ui.label(t("Library:", "Bibliothek:"));
            for (i, c) in names {
                if ui
                    .selectable_label(self.looked_up == Some(i), &c.name)
                    .clicked()
                {
                    self.looked_up = Some(i);
                }
            }
        });
    }

    /// Components with a hit for the looked-up substance, best first; returns the clicked one.
    fn found_list(
        &self,
        ui: &mut egui::Ui,
        run: &Run,
        lib: &Library,
        cutoff: f64,
    ) -> Option<usize> {
        let (i, ((..), id, _)) = (self.looked_up?, run.ident.as_ref()?);
        let found = find_substance(id, &lib.compounds, i, cutoff);
        let mut clicked = None;
        ui.horizontal_wrapped(|ui| {
            ui.strong(&lib.compounds[i].name);
            if found.is_empty() {
                ui.label(t("not found in this run", "in diesem Lauf nicht gefunden"));
            }
            for (k, score) in found.iter().take(12) {
                let at = run.data.scan_times_s[run.components[*k].apex_scan as usize] / 60.0;
                if ui.link(format!("{at:.2} min: {score:.3}")).clicked() {
                    clicked = Some(*k);
                }
            }
        });
        clicked
    }
}

/// Position in a list of `len` after one step of ↑ (-1) or ↓ (+1) from position `at`; without a
/// selection ↓ starts at the top and ↑ at the bottom.
pub(crate) fn stepped(len: usize, at: Option<usize>, step: i32) -> Option<usize> {
    let last = len.checked_sub(1)?;
    Some(match at {
        None if step > 0 => 0,
        None => last,
        Some(p) => p.saturating_add_signed(step as isize).min(last),
    })
}

/// The calibration whose RI the list shows: that of the identification, else `cal`.
fn calibration<'a>(run: &'a Run, cal: Option<&'a RiCalibration>) -> Option<&'a RiCalibration> {
    match &run.ident {
        Some(((_, used, _), _, shifted)) => shifted.as_ref().or(used.as_deref()),
        None => cal,
    }
}

/// RI of component `k` as the list shows it; see `calibration`.
pub fn component_ri(run: &Run, k: usize, cal: Option<&RiCalibration>) -> Option<f64> {
    let rt = run.data.scan_times_s[run.components[k].apex_scan as usize];
    calibration(run, cal).map(|c| c.ri(rt))
}

/// Rows of the active run. Without a library every row is "No match: N" (`result_rows`).
fn build_rows(run: &Run, cal: Option<&RiCalibration>, cutoff: f64) -> Vec<Row> {
    let unnamed;
    let (id, compounds) = match &run.ident {
        Some(((lib, ..), id, _)) => (id, &lib.compounds[..]),
        None => {
            let n = run.components.len();
            unnamed = Identification {
                hits: vec![Vec::new(); n],
                assigned: vec![None; n],
            };
            (&unnamed, &[][..])
        }
    };
    let times = &run.data.scan_times_s;
    let ris = calibration(run, cal).map(|c| c.component_ris(times, &run.components));
    let rows = result_rows(
        &run.components,
        times,
        ris.as_deref(),
        &run.peak_masses,
        id,
        compounds,
        cutoff,
    );
    (rows.into_iter())
        .map(|res| {
            let c = &run.components[res.component];
            let d = md_core::detail::component_detail(c, &run.peaks, &run.peak_masses);
            Row {
                area: (res.area)
                    .unwrap_or_else(|| md_core::results::area(c, &[], &run.peak_masses)),
                sn: (!d.peaks.is_empty()).then_some(d.mean_signal_to_noise),
                ions: c.spectrum.len(),
                res,
            }
        })
        .collect()
}

/// Click order of a column header: ascending, descending, off.
fn next_sort(sort: Option<(Col, bool)>, col: Col) -> Option<(Col, bool)> {
    match sort {
        Some((s, false)) if s == col => Some((col, true)),
        Some((s, true)) if s == col => None,
        _ => Some((col, false)),
    }
}

/// The table of `order` (indices into `rows`), scrolled to position `scroll` of `order`; returns the
/// component clicked in it.
fn table(
    ui: &mut egui::Ui,
    sort: &mut Option<(Col, bool)>,
    rows: &[Row],
    order: &[usize],
    run: &Run,
    lib: Option<&Library>,
    scroll: Option<usize>,
) -> Option<usize> {
    let cols = match lib {
        Some(_) => &WITH_LIBRARY[..],
        None => &PLAIN[..],
    };
    let mut clicked = None;
    // Scrolls sideways when the panel is narrower than the columns, as the batch table does.
    egui::ScrollArea::horizontal().show(ui, |ui| {
        ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
        let mut t = TableBuilder::new(ui)
            .id_salt("results table")
            .striped(true)
            .sense(egui::Sense::click())
            .resizable(true);
        if let Some(row) = scroll {
            t = t.scroll_to_row(row, None);
        }
        for c in cols {
            t = t.column(c.column());
        }
        t.header(20.0, |mut h| {
            for c in cols {
                h.col(|ui| {
                    ui.with_layout(c.align(), |ui| header_cell(ui, sort, *c));
                });
            }
        })
        .body(|body| {
            body.rows(18.0, order.len(), |mut row| {
                let r = &rows[order[row.index()]];
                row.set_selected(run.component == Some(r.res.component));
                for c in cols {
                    row.col(|ui| {
                        ui.with_layout(c.align(), |ui| {
                            let text = cell(r, *c);
                            // An unmatched component reads as a placeholder.
                            match r.res.compound.is_none() && *c == Col::Name {
                                true => ui.weak(text),
                                false => ui.label(text),
                            };
                        });
                    });
                }
                let resp = row.response();
                if resp.clicked() {
                    clicked = Some(r.res.component);
                }
                resp.on_hover_ui(|ui| hover(ui, &r.res, run, lib));
            });
        });
    });
    clicked
}

fn header_cell(ui: &mut egui::Ui, sort: &mut Option<(Col, bool)>, col: Col) {
    let arrow = match *sort {
        Some((s, d)) if s == col => [" ▲", " ▼"][d as usize],
        _ => "",
    };
    let text = egui::RichText::new(format!("{}{arrow}", col.title())).strong();
    let resp = ui
        .add(egui::Label::new(text).sense(egui::Sense::click()))
        .on_hover_cursor(egui::CursorIcon::PointingHand)
        .on_hover_text(t(
            "Click to sort: ascending, descending, off",
            "Klicken zum Sortieren: aufsteigend, absteigend, aus",
        ));
    if resp.clicked() {
        *sort = next_sort(*sort, col);
    }
}

/// Cell text as the old table showed it (a missing value as a dash).
fn cell(r: &Row, col: Col) -> String {
    use crate::detail::num;
    let res = &r.res;
    match col {
        Col::Name => shown_name(res),
        Col::Rt => format!("{:.2}", res.rt_s / 60.0),
        Col::Ri => num(res.ri, 1),
        Col::Area => format!("{:.0}", r.area),
        Col::Score => num(res.overall, 3),
        Col::Spectrum => num(res.spec_score, 3),
        Col::RiScore => num(res.ri_score, 3),
        Col::Sn => num(r.sn, 1),
        Col::Ions => r.ions.to_string(),
        Col::QuantIons => {
            let ions: Vec<String> = res.quant_ions.iter().map(f32::to_string).collect();
            ions.join(" ")
        }
    }
}

/// Missing values sort first.
fn compare(a: &Row, b: &Row, col: Col) -> std::cmp::Ordering {
    let num = |v: Option<f64>| v.unwrap_or(f64::NEG_INFINITY);
    let (x, y) = (&a.res, &b.res);
    match col {
        Col::Name => shown_name(x)
            .to_lowercase()
            .cmp(&shown_name(y).to_lowercase()),
        Col::Rt => x.rt_s.total_cmp(&y.rt_s),
        Col::Ri => num(x.ri).total_cmp(&num(y.ri)),
        Col::Area => a.area.total_cmp(&b.area),
        Col::Score => num(x.overall).total_cmp(&num(y.overall)),
        Col::Spectrum => num(x.spec_score).total_cmp(&num(y.spec_score)),
        Col::RiScore => num(x.ri_score).total_cmp(&num(y.ri_score)),
        Col::Sn => num(a.sn).total_cmp(&num(b.sn)),
        Col::Ions => a.ions.cmp(&b.ions),
        Col::QuantIons => x.quant_ions.len().cmp(&y.quant_ions.len()),
    }
}

/// md-core names a component without match "No match: N" for the CSV too; on screen it is translated.
fn shown_name(r: &ResultRow) -> String {
    match r.compound {
        None => tf!("No match: {}", "Kein Treffer: {}", r.component + 1),
        Some(_) => r.name.clone(),
    }
}

/// Old tooltip: scores, RI, RT and the alternative hits of the component.
fn hover(ui: &mut egui::Ui, r: &ResultRow, run: &Run, lib: Option<&Library>) {
    ui.strong(shown_name(r));
    let rt = r.rt_s / 60.0;
    ui.label(tf!("RT {rt:.2} min", "RT {rt:.2} min"));
    if let Some(ri) = r.ri {
        ui.label(format!("RI {ri:.1}"));
    }
    let (Some(((..), id, _)), Some(lib)) = (&run.ident, lib) else {
        return;
    };
    for h in id.hits[r.component].iter().take(5) {
        let mark = if Some(h.compound) == r.compound {
            "→ "
        } else {
            "   "
        };
        ui.label(format!(
            "{mark}{:.3}  {}",
            h.overall, lib.compounds[h.compound].name
        ));
    }
}

/// The file dialog for the saved components of `run`: `.cmp` next to the run, where it is found again.
fn dialog(run: &Run) -> rfd::FileDialog {
    let d = rfd::FileDialog::new()
        .add_filter(t("Saved components", "Gespeicherte Komponenten"), &["cmp"])
        .set_file_name(format!("{}.cmp", run.name));
    match run.path.parent() {
        Some(dir) => d.set_directory(dir),
        None => d,
    }
}

/// Writes the components of `run` to `path`; returns how many.
fn save_components(
    run: &Run,
    path: &std::path::Path,
    cal: Option<&RiCalibration>,
) -> Result<usize, String> {
    // The RI the list shows, so the old program lists it too.
    let mut ris = vec![-1.0; run.components.len()];
    for (k, ri) in ris.iter_mut().enumerate() {
        *ri = component_ri(run, k, cal).unwrap_or(-1.0);
    }
    let compounds = md_core::cmp::compounds_of_components(
        &run.components,
        &run.peak_masses,
        &run.peaks,
        &run.data.scan_times_s,
        &ris,
    );
    md_core::cmp::write_cmp(path, &compounds).map_err(|e| e.to_string())?;
    Ok(compounds.len())
}

/// Replaces the components of `run` by those in `path`, which must be of this run; returns how many.
/// The selection and hits belong to the old components and go.
fn load_components(run: &mut Run, path: &std::path::Path) -> Result<usize, String> {
    let (components, masses, peaks) =
        md_core::cmp::components_from_file(path, &run.data.scan_times_s)
            .map_err(|e| e.to_string())?;
    let n = components.len();
    (run.components, run.peak_masses, run.peaks) = (components, masses, peaks);
    (run.component, run.ident, run.ident_job) = (None, None, None);
    run.rev += 1;
    Ok(n)
}

/// Asks for a file and writes the table; returns the note to show.
fn export(run: &str, rows: &[ResultRow]) -> Msg {
    let Some(p) = rfd::FileDialog::new()
        .add_filter("CSV", &["csv"])
        .set_file_name(format!("{run}-results.csv"))
        .save_file()
    else {
        return Msg::default();
    };
    match std::fs::File::create(&p)
        .map(std::io::BufWriter::new)
        .and_then(|f| write_csv(f, rows))
    {
        Ok(()) => {
            let p = p.display();
            msg!("Saved {p}", "{p} gespeichert")
        }
        Err(e) => msg!("Export failed: {e}", "Export fehlgeschlagen: {e}"),
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use std::sync::Arc;

    use md_core::Chromatogram;
    use md_core::deconv::{Component, FittedPeak};
    use md_core::library::LibraryCompound;
    use md_core::search::{Identification, LibraryHit, PreparedLibrary, SearchSettings};

    use super::*;
    use crate::settings;

    pub(crate) fn run() -> Run {
        let s = SearchSettings::default();
        let compounds = vec![LibraryCompound::new("Alanine")];
        let lib = Arc::new(Library {
            path: "lib.mdlib".into(),
            name: "lib.mdlib".into(),
            prepared: PreparedLibrary::new(&compounds, &s),
            compounds,
            settings: s,
        });
        let hit = LibraryHit {
            compound: 0,
            spec_sim: 0.9,
            ri_sim: None,
            overall: 0.9,
        };
        let comp = |apex| Component {
            first_scan: 0,
            apex_scan: apex,
            last_scan: 2,
            model: Vec::new(),
            spectrum: Vec::new(),
            peaks: vec![FittedPeak {
                peak: 0,
                model_share: 1.0,
                integral: 5.0,
            }],
            exact_spectrum: Vec::new(),
        };
        Run {
            path: "a.cdf".into(),
            name: "a".into(),
            data: Arc::new(Chromatogram {
                scan_times_s: vec![0.0, 60.0, 120.0],
                offsets: vec![0; 4],
                mz: Vec::new(),
                intensity: Vec::new(),
                exact_masses: false,
                baseline: Default::default(),
            }),
            minutes: vec![0.0, 1.0, 2.0],
            tic: vec![0.0; 3],
            sic: Vec::new(),
            sic_for: None,
            sic_job: None,
            components: vec![comp(1), comp(2)],
            peak_masses: vec![73],
            // The peak both components point to.
            peaks: vec![md_core::peaks::Peak {
                row: 0,
                first_scan: 0,
                apex_scan: 1,
                last_scan: 2,
                precise_apex: None,
                signal_to_noise: 12.0,
                quality: 0.5,
                sharpness: 3.0,
            }],
            saved_calibration: None,
            analysis: settings::Analysis::default(),
            scan: 0,
            component: None,
            visible: true,
            ident: Some((
                (lib, None, None),
                Identification {
                    hits: vec![vec![hit], Vec::new()],
                    assigned: vec![Some(hit), None],
                },
                None,
            )),
            ident_job: None,
            rev: 0,
        }
    }

    fn draw(r: &mut Results, run: &mut Run, cal: Option<&RiCalibration>, events: Vec<egui::Event>) {
        let ctx = egui::Context::default();
        let input = egui::RawInput {
            events,
            ..Default::default()
        };
        ctx.run_ui(input, |ui| r.panel(ui, Some(run), cal, 0.0))
            .textures_delta
            .clear();
    }

    #[test]
    fn panel_draws_every_row_and_caches_them() {
        let mut run = run();
        let mut r = Results {
            filter: "ala".into(),
            looked_up: Some(0),
            ..Default::default()
        };
        for _ in 0..3 {
            draw(&mut r, &mut run, None, Vec::new());
        }
        let rows = &r.cache.as_ref().unwrap().1;
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].res.name, "Alanine");
        assert_eq!(rows[0].area, 5.0);
    }

    #[test]
    fn without_a_library_every_component_is_listed_with_its_ri() {
        let mut run = run();
        run.ident = None;
        let cal = RiCalibration::new([(0.0, 1000.0), (120.0, 1200.0)]).unwrap();
        let mut r = Results::default();
        draw(&mut r, &mut run, Some(&cal), Vec::new());
        let rows = &r.cache.as_ref().unwrap().1;
        let names: Vec<_> = rows.iter().map(|r| r.res.name.as_str()).collect();
        assert_eq!(names, ["No match: 1", "No match: 2"]);
        // The only peak (5.0) counts as the area also without a match.
        assert_eq!((rows[0].res.ri, rows[0].area), (Some(1100.0), 5.0));
    }

    #[test]
    fn arrow_keys_step_through_the_list() {
        assert_eq!(stepped(3, None, 1), Some(0));
        assert_eq!(stepped(3, None, -1), Some(2));
        assert_eq!(stepped(3, Some(2), 1), Some(2));
        assert_eq!(stepped(3, Some(0), -1), Some(0));
        assert_eq!(stepped(0, None, 1), None);
        let mut run = run();
        let mut r = Results::default();
        let down = egui::Event::Key {
            key: egui::Key::ArrowDown,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        };
        draw(&mut r, &mut run, None, vec![down.clone()]);
        draw(&mut r, &mut run, None, vec![down]);
        assert_eq!((run.component, run.scan), (Some(1), 2));
    }

    #[test]
    fn saved_components_load_again_with_their_peaks() {
        use md_core::peaks::Peak;

        let mut a = run();
        for k in &mut a.components {
            k.model = vec![1.0, 4.0, 1.0];
        }
        a.peaks = vec![Peak {
            row: 3,
            first_scan: 0,
            apex_scan: 1,
            last_scan: 2,
            precise_apex: Some(1.25),
            signal_to_noise: 12.0,
            quality: 0.5,
            sharpness: 3.0,
        }];
        let path = std::env::temp_dir().join(format!("md-app-save-{}.cmp", std::process::id()));
        assert_eq!(save_components(&a, &path, None), Ok(2));

        let mut b = run();
        b.component = Some(1);
        let rev = b.rev;
        assert_eq!(load_components(&mut b, &path), Ok(2));
        std::fs::remove_file(&path).ok();
        // Both components shared one peak; the file gives each its own.
        let mut want = a.components.clone();
        want[1].peaks[0].peak = 1;
        assert_eq!(b.components, want);
        assert_eq!((b.peak_masses.clone(), b.peaks.len()), (vec![73, 73], 2));
        // "Add component" takes its ions from these peaks (#312), whatever the current width.
        let ions = md_core::batch::ion_peaks_by(&b.components[1], &b.peaks, |i| b.peak_masses[i]);
        assert_eq!(
            ions.iter()
                .map(|p| (p.mass, p.integral))
                .collect::<Vec<_>>(),
            [(73, 5.0)]
        );
        assert_eq!(
            (b.peaks[0].signal_to_noise, b.peaks[0].precise_apex),
            (12.0, Some(1.25))
        );
        assert!(b.component.is_none() && b.ident.is_none() && b.rev == rev + 1);

        // A file of another run is refused and leaves the run as it was.
        let mut other = run();
        other.data = std::sync::Arc::new(Chromatogram {
            scan_times_s: vec![0.0, 61.0, 120.0],
            ..(*other.data).clone()
        });
        let path = std::env::temp_dir().join(format!("md-app-save2-{}.cmp", std::process::id()));
        save_components(&a, &path, None).unwrap();
        assert!(load_components(&mut other, &path).is_err());
        std::fs::remove_file(&path).ok();
        assert_eq!(other.components.len(), 2);
    }

    #[test]
    fn a_new_result_rebuilds_the_rows() {
        let mut run = run();
        let mut r = Results::default();
        draw(&mut r, &mut run, None, Vec::new());
        assert_eq!(r.cache.as_ref().unwrap().1.len(), 2);
        // The hits shrink with the components; a stale cache would index past them.
        run.components.truncate(1);
        let id = &mut run.ident.as_mut().unwrap().1;
        id.hits.truncate(1);
        id.assigned.truncate(1);
        run.rev += 1;
        draw(&mut r, &mut run, None, Vec::new());
        assert_eq!(r.cache.as_ref().unwrap().1.len(), 1);
        // Another run with the same `rev` (the active one was closed) gets its own rows.
        let mut other = super::tests::run();
        other.rev = run.rev;
        draw(&mut r, &mut other, None, Vec::new());
        assert_eq!(r.cache.as_ref().unwrap().1.len(), 2);
    }
}
