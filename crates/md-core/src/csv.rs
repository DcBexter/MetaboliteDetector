//! CSV output of plotted series (chromatogram and spectrum export). Values are written with Rust's
//! shortest round-trip formatting, so a file holds exactly the numbers the plots show.

use std::io::{self, Write};

/// Writes `header` and then one row per index of the equally long `columns`.
pub fn write_columns(out: &mut impl Write, header: &[&str], columns: &[&[f64]]) -> io::Result<()> {
    let rows = columns.first().map_or(0, |c| c.len());
    if header.len() != columns.len() || columns.iter().any(|c| c.len() != rows) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "CSV columns differ in length",
        ));
    }
    let names: Vec<String> = header.iter().map(|h| quote(h)).collect();
    writeln!(out, "{}", names.join(","))?;
    for r in 0..rows {
        for (i, c) in columns.iter().enumerate() {
            if i > 0 {
                out.write_all(b",")?;
            }
            write!(out, "{}", c[r])?;
        }
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// RFC 4180 quoting for names that contain a separator, quote or line break.
fn quote(s: &str) -> String {
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
