//! UI language (S8, #257): English or German. Every UI text is an English/German pair at its call
//! site, [`t`] for fixed texts and [`tf!`] for formatted ones, so no text can miss a language.
//! CSV headers, file formats and saved keys stay plain English literals, so scripts keep working.
//!
//! The language is per thread: the UI thread sets it from the saved choice, [`crate::jobs::spawn`]
//! hands it to its worker, and parallel tests do not see each other's choice. A thread starts in English.

use std::cell::Cell;

use eframe::egui;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Lang {
    #[default]
    En,
    De,
}

thread_local! {
    static CURRENT: Cell<Lang> = const { Cell::new(Lang::En) };
}

impl Lang {
    /// The value saved under `language` (S4).
    pub fn key(self) -> &'static str {
        match self {
            Lang::En => "en",
            Lang::De => "de",
        }
    }

    pub fn from_key(v: &str) -> Option<Lang> {
        [Lang::En, Lang::De].into_iter().find(|l| l.key() == v)
    }

    /// The system's language if it is German, else English.
    pub fn system() -> Lang {
        Self::from_locale(sys_locale::get_locale().as_deref())
    }

    /// "de", "de-DE", "de_AT.UTF-8" are German; anything else or nothing is English.
    fn from_locale(locale: Option<&str>) -> Lang {
        let prefix = locale.and_then(|l| l.get(..2));
        if prefix.is_some_and(|p| p.eq_ignore_ascii_case("de")) {
            Lang::De
        } else {
            Lang::En
        }
    }
}

pub fn get() -> Lang {
    CURRENT.get()
}

pub fn set(lang: Lang) {
    CURRENT.set(lang);
}

pub fn de() -> bool {
    get() == Lang::De
}

/// The text in the current language.
pub fn t(en: &'static str, de: &'static str) -> &'static str {
    if self::de() { de } else { en }
}

/// `format!` in the current language: `tf!("Step {x} of {n}", "Schritt {x} von {n}")`.
macro_rules! tf {
    ($en:literal, $de:literal $(, $arg:expr)* $(,)?) => {
        if $crate::lang::de() {
            format!($de $(, $arg)*)
        } else {
            format!($en $(, $arg)*)
        }
    };
}
pub(crate) use tf;

/// The switch in Help and in Settings. Language names stay in their own language.
pub fn switch(ui: &mut egui::Ui) {
    let mut lang = get();
    ui.horizontal(|ui| {
        ui.label(t("Language", "Sprache"));
        ui.radio_value(&mut lang, Lang::En, "English");
        ui.radio_value(&mut lang, Lang::De, "Deutsch");
    })
    .response
    .on_hover_text(t(
        "The language of the window texts. CSV headers and files stay English.",
        "Die Sprache der Fenstertexte. CSV-Köpfe und Dateien bleiben englisch.",
    ));
    set(lang);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn locale_picks_german_only_for_de() {
        for (locale, lang) in [
            (Some("de-DE"), Lang::De),
            (Some("de_AT.UTF-8"), Lang::De),
            (Some("DE"), Lang::De),
            (Some("en-US"), Lang::En),
            (Some("fr"), Lang::En),
            (Some("d"), Lang::En),
            (None, Lang::En),
        ] {
            assert_eq!(Lang::from_locale(locale), lang, "{locale:?}");
        }
    }

    #[test]
    fn texts_follow_the_language() {
        let n = 3;
        set(Lang::De);
        assert_eq!(t("Open", "Öffnen"), "Öffnen");
        assert_eq!(tf!("{n} runs", "{n} Läufe"), "3 Läufe");
        set(Lang::En);
        assert_eq!(tf!("{n} runs", "{n} Läufe"), "3 runs");
    }
}
