//! Natural isotope patterns of molecular formulas (`M+0`, `M+1`, …), the old "Natural Isotope
//! Calculator" (`tools::NaturalMIDCalculator`, features F16/F17). The MID evaluation builds on it.
//!
//! The pattern is the product of one polynomial per atom, `Σ abundance · x^(nominal − lightest)`,
//! truncated after `len` terms; an element with `n` atoms is raised to the power `n` by squaring.
//! Index `i` of the result is the abundance of the ions that are `i` nominal masses heavier than
//! the all-lightest-isotope ion, like the old calculator. The old program prints the same values
//! for the same table (checked in the tests against its output).

use std::collections::{BTreeMap, HashMap};

/// Element symbol → number of atoms.
pub type Formula = BTreeMap<String, u32>;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Isotope {
    /// Nominal (integer) mass.
    pub nominal: u32,
    /// Exact mass in u.
    pub mass: f64,
    /// Share of this isotope in the element, 0..=1; the shares of an element sum to 1.
    pub abundance: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopePattern {
    /// Nominal mass of the all-lightest-isotope ion; the ion of `abundance[i]` is `nominal_mass + i`.
    pub nominal_mass: u64,
    /// Exact mass of the all-lightest-isotope ion (for most elements the monoisotopic mass).
    pub monoisotopic_mass: f64,
    /// Abundances of `M+0`, `M+1`, …, summing to 1 when `len` covers the whole pattern.
    pub abundance: Vec<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IsotopeTable(HashMap<String, Vec<Isotope>>);

impl IsotopeTable {
    /// Rows are (symbol, nominal mass, exact mass, abundance in percent). Abundances of an element
    /// are normalised to sum to 1; rows of one element may come in any order.
    pub fn new<'a>(rows: impl IntoIterator<Item = (&'a str, u32, f64, f64)>) -> Self {
        let mut map: HashMap<String, Vec<Isotope>> = HashMap::new();
        for (symbol, nominal, mass, percent) in rows {
            map.entry(symbol.to_owned()).or_default().push(Isotope {
                nominal,
                mass,
                abundance: percent,
            });
        }
        for isotopes in map.values_mut() {
            isotopes.sort_by_key(|i| i.nominal);
            let total: f64 = isotopes.iter().map(|i| i.abundance).sum();
            isotopes.iter_mut().for_each(|i| i.abundance /= total);
        }
        Self(map)
    }

    /// Natural abundances and masses of the elements common in GC-MS and metabolomics (IUPAC 2016
    /// isotopic compositions, atomic masses as in the NIST compilation).
    pub fn standard() -> Self {
        Self::new(STANDARD.iter().copied())
    }

    pub fn isotopes(&self, symbol: &str) -> Option<&[Isotope]> {
        self.0.get(symbol).map(Vec::as_slice)
    }

    /// Isotope pattern of `formula`, the first `len` ions `M+0` … `M+len-1`.
    pub fn pattern(&self, formula: &Formula, len: usize) -> Result<IsotopePattern, String> {
        let mut total = vec![0.0; len];
        if len > 0 {
            total[0] = 1.0;
        }
        let (mut nominal_mass, mut monoisotopic_mass) = (0u64, 0.0);
        for (symbol, &count) in formula {
            let isotopes = self
                .isotopes(symbol)
                .ok_or_else(|| format!("unknown element {symbol}"))?;
            let lightest = isotopes[0];
            nominal_mass += u64::from(count) * u64::from(lightest.nominal);
            monoisotopic_mass += f64::from(count) * lightest.mass;
            let mut atom = vec![0.0; len];
            for i in isotopes {
                if let Some(slot) = atom.get_mut((i.nominal - lightest.nominal) as usize) {
                    *slot += i.abundance;
                }
            }
            total = multiply(&total, &power(&atom, count, len), len);
        }
        Ok(IsotopePattern {
            nominal_mass,
            monoisotopic_mass,
            abundance: total,
        })
    }
}

/// `a · b` truncated to `len` terms.
fn multiply(a: &[f64], b: &[f64], len: usize) -> Vec<f64> {
    let mut out = vec![0.0; len];
    for (i, x) in a.iter().enumerate().filter(|(_, x)| **x != 0.0) {
        for (o, y) in out[i..].iter_mut().zip(b) {
            *o += x * y;
        }
    }
    out
}

/// `base^n` truncated to `len` terms, by squaring.
fn power(base: &[f64], mut n: u32, len: usize) -> Vec<f64> {
    let mut result = vec![0.0; len];
    if len > 0 {
        result[0] = 1.0;
    }
    let mut square = base.to_vec();
    while n > 0 {
        if n & 1 == 1 {
            result = multiply(&result, &square, len);
        }
        n >>= 1;
        if n > 0 {
            square = multiply(&square, &square, len);
        }
    }
    result
}

/// Parses a molecular formula such as `C6H12O6` or `Ca(OH)2`: element symbols (an upper-case
/// letter and at most one lower-case letter), optional counts, parentheses with an optional
/// multiplier. Repeated elements add up (`CH3CH2OH` = `C2H6O`); a count of 0 is allowed. The old
/// calculator has no parentheses, so every formula it accepts means the same here.
pub fn parse_formula(text: &str) -> Result<Formula, String> {
    let bytes = text.as_bytes();
    let mut stack = vec![Formula::new()];
    let mut i = 0;
    // Optional count at `i`, default 1.
    let count = |i: &mut usize| -> Result<u32, String> {
        let start = *i;
        while bytes.get(*i).is_some_and(u8::is_ascii_digit) {
            *i += 1;
        }
        if start == *i {
            return Ok(1);
        }
        text[start..*i]
            .parse()
            .map_err(|_| format!("count too large: {}", &text[start..*i]))
    };
    let add = |group: &mut Formula, symbol: String, n: u32| -> Result<(), String> {
        let slot = group.entry(symbol).or_default();
        *slot = slot.checked_add(n).ok_or("atom count too large")?;
        Ok(())
    };
    while i < bytes.len() {
        match bytes[i] {
            b'A'..=b'Z' => {
                let start = i;
                i += 1;
                if bytes.get(i).is_some_and(u8::is_ascii_lowercase) {
                    i += 1;
                }
                let symbol = text[start..i].to_owned();
                let n = count(&mut i)?;
                add(stack.last_mut().expect("never empty"), symbol, n)?;
            }
            b'(' => {
                stack.push(Formula::new());
                i += 1;
            }
            b')' => {
                i += 1;
                if stack.len() == 1 {
                    return Err("unmatched )".into());
                }
                let group = stack.pop().expect("checked");
                let n = count(&mut i)?;
                for (symbol, k) in group {
                    add(
                        stack.last_mut().expect("never empty"),
                        symbol,
                        k.checked_mul(n).ok_or("atom count too large")?,
                    )?;
                }
            }
            c => return Err(format!("unexpected {:?} at position {}", c as char, i + 1)),
        }
    }
    if stack.len() > 1 {
        return Err("unmatched (".into());
    }
    let formula = stack.pop().expect("never empty");
    if formula.is_empty() {
        return Err("empty formula".into());
    }
    Ok(formula)
}

/// (symbol, nominal mass, exact mass, abundance in percent).
const STANDARD: &[(&str, u32, f64, f64)] = &[
    ("H", 1, 1.00782503223, 99.9885),
    ("H", 2, 2.01410177812, 0.0115),
    ("Li", 6, 6.0151228874, 7.59),
    ("Li", 7, 7.0160034366, 92.41),
    ("B", 10, 10.01293695, 19.9),
    ("B", 11, 11.00930536, 80.1),
    ("C", 12, 12.0, 98.93),
    ("C", 13, 13.00335483507, 1.07),
    ("N", 14, 14.00307400443, 99.636),
    ("N", 15, 15.00010889888, 0.364),
    ("O", 16, 15.99491461957, 99.757),
    ("O", 17, 16.9991317565, 0.038),
    ("O", 18, 17.99915961286, 0.205),
    ("F", 19, 18.99840316273, 100.0),
    ("Na", 23, 22.989769282, 100.0),
    ("Mg", 24, 23.985041697, 78.99),
    ("Mg", 25, 24.985836976, 10.00),
    ("Mg", 26, 25.982592968, 11.01),
    ("Al", 27, 26.98153853, 100.0),
    ("Si", 28, 27.97692653465, 92.223),
    ("Si", 29, 28.9764946649, 4.685),
    ("Si", 30, 29.973770136, 3.092),
    ("P", 31, 30.97376199842, 100.0),
    ("S", 32, 31.9720711744, 94.99),
    ("S", 33, 32.9714589098, 0.75),
    ("S", 34, 33.967867004, 4.25),
    ("S", 36, 35.96708071, 0.01),
    ("Cl", 35, 34.968852682, 75.76),
    ("Cl", 37, 36.965902602, 24.24),
    ("K", 39, 38.9637064864, 93.2581),
    ("K", 40, 39.963998166, 0.0117),
    ("K", 41, 40.9618252579, 6.7302),
    ("Ca", 40, 39.962590863, 96.941),
    ("Ca", 42, 41.95861783, 0.647),
    ("Ca", 43, 42.95876644, 0.135),
    ("Ca", 44, 43.95548156, 2.086),
    ("Ca", 46, 45.9536893, 0.004),
    ("Ca", 48, 47.95252276, 0.187),
    ("Mn", 55, 54.93804391, 100.0),
    ("Fe", 54, 53.93960899, 5.845),
    ("Fe", 56, 55.93493633, 91.754),
    ("Fe", 57, 56.93539284, 2.119),
    ("Fe", 58, 57.93327443, 0.282),
    ("Co", 59, 58.93319429, 100.0),
    ("Cu", 63, 62.92959772, 69.15),
    ("Cu", 65, 64.9277897, 30.85),
    ("Zn", 64, 63.92914201, 49.17),
    ("Zn", 66, 65.92603381, 27.73),
    ("Zn", 67, 66.92712775, 4.04),
    ("Zn", 68, 67.92484455, 18.45),
    ("Zn", 70, 69.9253192, 0.61),
    ("Ge", 70, 69.92424875, 20.57),
    ("Ge", 72, 71.922075826, 27.45),
    ("Ge", 73, 72.923458956, 7.75),
    ("Ge", 74, 73.921177761, 36.50),
    ("Ge", 76, 75.921402726, 7.73),
    ("Se", 74, 73.922475934, 0.89),
    ("Se", 76, 75.919213704, 9.37),
    ("Se", 77, 76.919914154, 7.63),
    ("Se", 78, 77.91730928, 23.77),
    ("Se", 80, 79.9165218, 49.61),
    ("Se", 82, 81.9166995, 8.73),
    ("Br", 79, 78.9183376, 50.69),
    ("Br", 81, 80.9162897, 49.31),
    ("I", 127, 126.9044719, 100.0),
    ("Mo", 92, 91.90680796, 14.53),
    ("Mo", 94, 93.9050849, 9.15),
    ("Mo", 95, 94.90583877, 15.84),
    ("Mo", 96, 95.90467612, 16.67),
    ("Mo", 97, 96.90601812, 9.60),
    ("Mo", 98, 97.90540482, 24.39),
    ("Mo", 100, 99.9074718, 9.82),
    ("Ag", 107, 106.9050916, 51.839),
    ("Ag", 109, 108.9047553, 48.161),
    ("Hg", 196, 195.9658326, 0.15),
    ("Hg", 198, 197.9667686, 10.04),
    ("Hg", 199, 198.9682806, 16.94),
];

#[cfg(test)]
mod tests {
    use super::*;

    fn pattern(table: &IsotopeTable, formula: &str, len: usize) -> IsotopePattern {
        table
            .pattern(&parse_formula(formula).unwrap(), len)
            .unwrap()
    }

    fn assert_close(actual: &[f64], expected: &[f64]) {
        assert_eq!(actual.len(), expected.len());
        for (a, e) in actual.iter().zip(expected) {
            assert!(
                (a - e).abs() <= 1e-9 * e.abs().max(1e-12),
                "{actual:?} != {expected:?}"
            );
        }
    }

    /// The elements and values of the old program's `isotopes.csv` that the references below use.
    fn old_table() -> IsotopeTable {
        IsotopeTable::new([
            ("H", 1, 1.007825, 99.9844),
            ("H", 2, 2.0141, 0.0156),
            ("C", 12, 12.0, 98.9184),
            ("C", 13, 13.0034, 1.0816),
            ("N", 14, 14.0031, 99.6337),
            ("N", 15, 15.0001, 0.3663),
            ("O", 16, 15.9949, 99.758),
            ("O", 17, 16.9991, 0.038),
            ("O", 18, 17.9992, 0.204),
            ("Si", 28, 27.977, 92.22),
            ("Si", 29, 28.976, 4.69),
            ("Si", 30, 29.974, 3.09),
            ("Fe", 54, 53.94, 5.81),
            ("Fe", 56, 55.935, 91.64),
            ("Fe", 57, 56.935, 2.21),
            ("Fe", 58, 57.933, 0.34),
        ])
    }

    /// Output of `tools::NaturalMIDCalculator::calculateMID` of MetaboliteDetector 3.3 (called
    /// through `LD_PRELOAD` with the table above), length 6.
    #[test]
    fn matches_old_program() {
        let table = old_table();
        let cases: [(&str, [f64; 6]); 3] = [
            (
                "C6H12O6",
                [
                    0.921585933265,
                    0.0642929801922,
                    0.0132191399165,
                    0.00081601144871,
                    8.12756856304e-05,
                    4.36800507493e-06,
                ],
            ),
            (
                "C12H32O6Si4",
                [
                    0.622506230612,
                    0.212844889259,
                    0.123791184504,
                    0.0298982785304,
                    0.0089921905006,
                    0.00160565029305,
                ],
            ),
            (
                "CH3CH2OH",
                [
                    0.975203762877,
                    0.0226106823072,
                    0.00213962730461,
                    4.56472005392e-05,
                    2.80071227666e-07,
                    2.3930544513e-10,
                ],
            ),
        ];
        for (formula, expected) in cases {
            assert_close(&pattern(&table, formula, 6).abundance, &expected);
        }
        let c6 = pattern(&table, "C6H12O6", 6);
        assert_eq!(c6.nominal_mass, 180);
        assert!((c6.monoisotopic_mass - 180.0633).abs() < 1e-6);
        // Isotopes without a neighbour leave gaps, indexed by nominal mass from the lightest.
        let fe = pattern(&table, "Fe", 12).abundance;
        assert_close(&fe[..5], &[0.0581, 0.0, 0.9164, 0.0221, 0.0034]);
        assert!(fe[5..].iter().all(|v| *v == 0.0));
    }

    #[test]
    fn glucose_matches_textbook() {
        let p = pattern(&IsotopeTable::standard(), "C6H12O6", 4);
        assert_eq!(p.nominal_mass, 180);
        assert!((p.monoisotopic_mass - 180.063388).abs() < 1e-5);
        // Relative to M: M+1 ≈ 6.86 % (6 · 1.07/98.93 for C13 + 12 · 0.0115/99.99 for H2 + 6 · 0.038/99.76 for O17), M+2 ≈ 1.43 %.
        let rel: Vec<f64> = p
            .abundance
            .iter()
            .map(|a| 100.0 * a / p.abundance[0])
            .collect();
        assert!(
            (rel[1] - 6.86).abs() < 0.02 && (rel[2] - 1.43).abs() < 0.03,
            "{rel:?}"
        );
    }

    #[test]
    fn halogens_and_tms() {
        let table = IsotopeTable::standard();
        // Cl2 ≈ 100 : 64 : 10, Br2 ≈ 1 : 2 : 1.
        let cl2 = pattern(&table, "Cl2", 5).abundance;
        assert_close(
            &cl2,
            &[
                0.75760_f64.powi(2),
                0.0,
                2.0 * 0.7576 * 0.2424,
                0.0,
                0.2424_f64.powi(2),
            ],
        );
        let br2 = pattern(&table, "Br2", 5).abundance;
        assert!((br2[2] / br2[0] - 1.95).abs() < 0.01 && (br2[4] / br2[0] - 0.95).abs() < 0.01);
        // Each Si adds M+2 (Si30): four TMS groups give far more M+2 than one.
        let tms1 = pattern(&table, "C3H9Si", 3).abundance;
        let tms4 = pattern(&table, "C12H36Si4", 3).abundance;
        assert!(tms4[2] > 3.0 * tms1[2]);
        // The complete pattern sums to 1; truncation never exceeds that.
        let all: f64 = pattern(&table, "C12H32O6Si4", 80).abundance.iter().sum();
        assert!((all - 1.0).abs() < 1e-12);
    }

    #[test]
    fn large_counts_are_fast_and_exact() {
        // Binomial: M+1 of C1000 with p = 1.07 %.
        let table = IsotopeTable::standard();
        let a = pattern(&table, "C1000", 3).abundance;
        let (p, q) = (0.0107, 0.9893_f64);
        assert!((a[0] - q.powi(1000)).abs() < 1e-12);
        assert!((a[1] - 1000.0 * p * q.powi(999)).abs() < 1e-12);
        let _ = pattern(&table, "C4000000000H4000000000", 30);
    }

    #[test]
    fn parses_formulas() {
        let f = |s| parse_formula(s).unwrap().into_iter().collect::<Vec<_>>();
        let pair = |s: &str, n| (s.to_owned(), n);
        assert_eq!(f("C6H12O6"), [pair("C", 6), pair("H", 12), pair("O", 6)]);
        assert_eq!(f("CH3CH2OH"), [pair("C", 2), pair("H", 6), pair("O", 1)]);
        assert_eq!(f("Ca(OH)2"), [pair("Ca", 1), pair("H", 2), pair("O", 2)]);
        assert_eq!(
            f("((CH3)3Si)2O"),
            [pair("C", 6), pair("H", 18), pair("O", 1), pair("Si", 2)]
        );
        assert_eq!(f("C0"), [pair("C", 0)]);
        for bad in [
            "",
            "c6",
            "C6 H",
            "(CH3",
            "CH3)",
            "C6+",
            "C99999999999",
            "(C4000000000)2",
            "2C",
        ] {
            assert!(parse_formula(bad).is_err(), "{bad:?} must not parse");
        }
        let table = IsotopeTable::standard();
        assert!(table.pattern(&parse_formula("Xx2").unwrap(), 3).is_err());
    }
}
