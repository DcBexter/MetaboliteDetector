//! CSV output of plotted series (chromatogram and spectrum export). Values are written with Rust's
//! shortest round-trip formatting, so a file holds exactly the numbers the plots show.

use std::io::{self, Write};

/// Writes `header` and then one row per index of the equally long `columns`.
pub fn write_columns(out: &mut impl Write, header: &[&str], columns: &[&[f64]]) -> io::Result<()> {
    let names: Vec<String> = header.iter().map(|h| quote(h)).collect();
    write_rows(out, &names, columns, ',')
}

/// The same table tab-separated, as spreadsheets paste it; tabs and line breaks in names become spaces.
pub fn write_tsv(out: &mut impl Write, header: &[&str], columns: &[&[f64]]) -> io::Result<()> {
    let names: Vec<String> = header
        .iter()
        .map(|h| h.replace(['\t', '\n', '\r'], " "))
        .collect();
    write_rows(out, &names, columns, '\t')
}

fn write_rows(
    out: &mut impl Write,
    names: &[String],
    columns: &[&[f64]],
    sep: char,
) -> io::Result<()> {
    let rows = columns.first().map_or(0, |c| c.len());
    if names.len() != columns.len() || columns.iter().any(|c| c.len() != rows) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "CSV columns differ in length",
        ));
    }
    writeln!(out, "{}", names.join(&sep.to_string()))?;
    for r in 0..rows {
        for (i, c) in columns.iter().enumerate() {
            if i > 0 {
                write!(out, "{sep}")?;
            }
            write!(out, "{}", c[r])?;
        }
        writeln!(out)?;
    }
    Ok(())
}

/// Rewrites a CSV of this crate for a spreadsheet with German locale, where `,` is the decimal mark and
/// `;` the separator: numeric cells get a decimal comma, cells with `;` get quoted. The header line
/// and the first `text_cols` columns (names, which must stay as they are) are not converted.
pub fn to_excel_de(csv: &str, text_cols: usize) -> String {
    let (mut out, mut cell) = (String::with_capacity(csv.len()), String::new());
    let (mut quoted, mut header, mut col) = (false, true, 0);
    let flush = |cell: &mut String, end: Option<char>, numeric: bool, out: &mut String| {
        let text = std::mem::take(cell);
        if text.starts_with('"') {
            out.push_str(&text);
        } else if numeric && text.parse::<f64>().is_ok() {
            out.push_str(&text.replace('.', ","));
        } else if text.contains(';') {
            out.push_str(&format!("\"{text}\""));
        } else {
            out.push_str(&text);
        }
        out.extend(end);
    };
    for c in csv.chars() {
        let numeric = !header && col >= text_cols;
        match c {
            '"' => {
                quoted = !quoted;
                cell.push(c);
            }
            ',' if !quoted => {
                flush(&mut cell, Some(';'), numeric, &mut out);
                col += 1;
            }
            '\n' if !quoted => {
                flush(&mut cell, Some('\n'), numeric, &mut out);
                (header, col) = (false, 0);
            }
            _ => cell.push(c),
        }
    }
    flush(&mut cell, None, !header && col >= text_cols, &mut out);
    out
}

/// `to_excel_de` applied to a file in place.
pub fn excel_de_file(path: &std::path::Path, text_cols: usize) -> io::Result<()> {
    std::fs::write(
        path,
        to_excel_de(&std::fs::read_to_string(path)?, text_cols),
    )
}

/// RFC 4180 quoting for names that contain a separator, quote or line break.
pub fn quote(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_owned()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tsv_uses_tabs() {
        let mut out = Vec::new();
        write_tsv(&mut out, &["a\tb", "c"], &[&[1.5], &[2.0]]).unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), "a b\tc\n1.5\t2\n");
    }

    #[test]
    fn excel_de_has_semicolons_and_decimal_commas() {
        let csv = "\u{feff}Substance,RT_s,Quant_ions,a.cdf\n\"A, \"\"x\"\"\",293.808,73 147,3643496.25\nB;1,1e3,,\n\"q\nr\",-0.5,1,inf\n";
        assert_eq!(
            to_excel_de(csv, 1),
            "\u{feff}Substance;RT_s;Quant_ions;a.cdf\n\"A, \"\"x\"\"\";293,808;73 147;3643496,25\n\"B;1\";1e3;;\n\"q\nr\";-0,5;1;inf\n"
        );
        // Names that look like numbers (a group `0.5`) stay as they are.
        assert_eq!(
            to_excel_de("Run,Group\nr1.cdf,0.5\n", 2),
            "Run;Group\nr1.cdf;0.5\n"
        );
    }

    #[test]
    fn rows_match_the_input_values() {
        let (t, y) = ([0.5, 1.0 / 3.0], [1.0e7, 0.25]);
        let mut out = Vec::new();
        write_columns(
            &mut out,
            &["Retention time (min)", "TIC, \"raw\""],
            &[&t, &y],
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();
        let mut lines = text.lines();
        assert_eq!(
            lines.next(),
            Some("Retention time (min),\"TIC, \"\"raw\"\"\"")
        );
        // Every value parses back to the exact f64 that was plotted.
        for (line, (t, y)) in lines.zip(t.iter().zip(&y)) {
            let (a, b) = line.split_once(',').unwrap();
            assert_eq!(
                (a.parse::<f64>().unwrap(), b.parse::<f64>().unwrap()),
                (*t, *y)
            );
        }
    }

    #[test]
    fn ragged_columns_are_rejected() {
        let mut out = Vec::new();
        assert!(write_columns(&mut out, &["a", "b"], &[&[1.0], &[1.0, 2.0]]).is_err());
        assert!(write_columns(&mut out, &["a"], &[&[1.0], &[2.0]]).is_err());
        assert!(out.is_empty());
    }
}
