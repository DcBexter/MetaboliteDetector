# Glossar

Domänenbegriffe des Projekts, wie sie in Code, Issues und docs/architecture.md benutzt werden.
Architekturbegriffe (module, interface, seam, adapter, deep/shallow, leverage, locality)
folgen dem Skill `codebase-design` und werden nicht übersetzt.

| Begriff | Bedeutung | Code |
|---------|-----------|------|
| **Run** | Ein GC-MS-Lauf: eine Datei (ANDI-MS `.cdf`, FastFlight, alter Cache) mit Scans, Zeiten und Spektren. | `md_core::Chromatogram` |
| **Scan** | Ein Spektrum zu einem Zeitpunkt des Runs; m/z innerhalb des Scans sortiert. | `Chromatogram` (flache Arrays mit Scan-Offsets) |
| **Ion trace** | Intensität einer (nominalen oder exakten) Masse über alle Scans. Einmal pro Run als ion-major Index gebaut. | `md_core::peaks::IonTraces`, `trace(row)` |
| **Peak** | Ein Maximum in einem Ion trace mit Grenzen, Apex, interpoliertem Apex und S/N. | `md_core::peaks::Peak` |
| **Peak pass** | Welche Peak-Erkennung ein Aufrufer will: *Reordered* (Hauptfenster und Batch, `detect_all_peaks` mit `reorder_peaks` auf Exact-Mass-Runs) oder *Fresh* (Rescan, `detect_peaks` wie das alte `ResultSeries::scanIons`). | `md_core::analysis::PeakPass` |
| **Component** | Ein deconvolvierter Stoff: Peaks mehrerer Ion traces mit gemeinsamer Retentionszeit, Modellpeak und Spektrum. | `md_core::deconv::Component` |
| **Analyze run** | Die Kette Ion traces, Peaks, Deconvolution für einen Run hinter einer Interface. | `md_core::analysis::analyze_run`; Peak-Stufe allein: `find_run_peaks` |
| **Library** | Sammlung von Referenzspektren mit RI, Quant-Ionen und FRAG-Angaben (`.lbr`, `.msl`, `.mdlib`). | `md_core::library` |
| **Library search** | Zuordnung von Components zu Library-Einträgen über Spektrum- und RI-Score, greedy eins zu eins. | `md_core::search` |
| **Retention index (RI)** | Zeitnormierung über Alkan-Marker; Kalibrierung RT zu RI. | `md_core::ri` |
| **Study** | Viele Runs einer Fragestellung, die gemeinsam quantifiziert werden: Library, Settings, Sample table, Gruppen. | `md_core::batch`; geplant: `run_study` (#363) |
| **Batch** | Die Quantifizierung einer Study: Substanz × Run Tabelle mit Quant-Ionen. | `md_core::batch::{scan_runs, complete_scan, table}` |
| **Rescan** | Zweiter Pass des Batch: fehlende Substanzen werden je Run an der erwarteten RT mit frischen Peaks nachgesucht. | `md_core::batch::rescan` |
| **Quant ion** | Die Masse, deren Peakfläche eine Substanz im Batch quantifiziert. | `QuantSettings` |
| **Normalisierung** | Umrechnung der Batch-Tabelle (interner Standard, Gesamtsignal, Pool, z-Score). | `md_core::normalize` |
| **MID** | Mass isotopomer distribution: Verteilung der Markierung je Fragment, korrigiert um natürliche Isotope. | `md_core::mid` |
| **Job** | Hintergrundarbeit der App mit Progress und Cancel; die UI pollt den Zustand. | `md_app::jobs::Job`; `md_core::job::{Cancel, Progress}` |
| **View model** | Fertig sortierte, gefilterte und formatierte Zeilen einer Tabelle, aus der Core-Tabelle plus Ansichtsoptionen gebaut und per Generation gecacht. | `md_core::view` (`BatchViewCache`) |
| **Legacy numeric** | Benannte Rundungen und Rechenwege, die 3.3 bit für bit nachbilden (C-Trunkierung, `round_sig`, glibc `pow`/`hypot`). | geplant: `md_core::legacy_numeric` (#365) |
| **Golden test** | Vergleich der neuen Kette gegen Ergebnisse des alten Binary 3.3 auf den Beispiel-Runs. | `crates/md-core/tests/golden_old_md.rs` |
| **Settings** | Die `AN_*`/`LS_*`-Parameter der Analyse; eine Quelle in md-core. | `md_core::settings::GcmsSettings`, `KEYS` (Key-Tabelle) |
