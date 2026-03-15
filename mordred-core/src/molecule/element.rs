/// Chemical element data.
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Element {
    H,
    He,
    Li,
    Be,
    B,
    C,
    N,
    O,
    F,
    Ne,
    Na,
    Mg,
    Al,
    Si,
    P,
    S,
    Cl,
    Ar,
    K,
    Ca,
    Br,
    I,
    Fe,
    Cu,
    Zn,
    Se,
}

impl Element {
    /// Atomic number.
    pub fn atomic_number(self) -> u8 {
        match self {
            Self::H => 1,
            Self::He => 2,
            Self::Li => 3,
            Self::Be => 4,
            Self::B => 5,
            Self::C => 6,
            Self::N => 7,
            Self::O => 8,
            Self::F => 9,
            Self::Ne => 10,
            Self::Na => 11,
            Self::Mg => 12,
            Self::Al => 13,
            Self::Si => 14,
            Self::P => 15,
            Self::S => 16,
            Self::Cl => 17,
            Self::Ar => 18,
            Self::K => 19,
            Self::Ca => 20,
            Self::Br => 35,
            Self::I => 53,
            Self::Fe => 26,
            Self::Cu => 29,
            Self::Zn => 30,
            Self::Se => 34,
        }
    }

    /// Average atomic weight.
    pub fn atomic_weight(self) -> f64 {
        match self {
            Self::H => 1.008,
            Self::He => 4.003,
            Self::Li => 6.941,
            Self::Be => 9.012,
            Self::B => 10.81,
            Self::C => 12.011,
            Self::N => 14.007,
            Self::O => 15.999,
            Self::F => 18.998,
            Self::Ne => 20.180,
            Self::Na => 22.990,
            Self::Mg => 24.305,
            Self::Al => 26.982,
            Self::Si => 28.086,
            Self::P => 30.974,
            Self::S => 32.06,
            Self::Cl => 35.45,
            Self::Ar => 39.948,
            Self::K => 39.098,
            Self::Ca => 40.078,
            Self::Br => 79.904,
            Self::I => 126.904,
            Self::Fe => 55.845,
            Self::Cu => 63.546,
            Self::Zn => 65.38,
            Self::Se => 78.971,
        }
    }

    /// Default valence for implicit hydrogen calculation.
    pub fn default_valence(self) -> u8 {
        match self {
            Self::H => 1,
            Self::B => 3,
            Self::C => 4,
            Self::N => 3,
            Self::O => 2,
            Self::F => 1,
            Self::Si => 4,
            Self::P => 3,
            Self::S => 2,
            Self::Cl => 1,
            Self::Br => 1,
            Self::I => 1,
            Self::Se => 2,
            _ => 0,
        }
    }

    /// Parse element from symbol string.
    pub fn from_symbol(s: &str) -> Option<Self> {
        match s {
            "H" => Some(Self::H),
            "He" => Some(Self::He),
            "Li" => Some(Self::Li),
            "Be" => Some(Self::Be),
            "B" => Some(Self::B),
            "C" | "c" => Some(Self::C),
            "N" | "n" => Some(Self::N),
            "O" | "o" => Some(Self::O),
            "F" => Some(Self::F),
            "Ne" => Some(Self::Ne),
            "Na" => Some(Self::Na),
            "Mg" => Some(Self::Mg),
            "Al" => Some(Self::Al),
            "Si" => Some(Self::Si),
            "P" | "p" => Some(Self::P),
            "S" | "s" => Some(Self::S),
            "Cl" => Some(Self::Cl),
            "Ar" => Some(Self::Ar),
            "K" => Some(Self::K),
            "Ca" => Some(Self::Ca),
            "Br" => Some(Self::Br),
            "I" => Some(Self::I),
            "Fe" => Some(Self::Fe),
            "Cu" => Some(Self::Cu),
            "Zn" => Some(Self::Zn),
            "Se" => Some(Self::Se),
            _ => None,
        }
    }

    /// Monoisotopic (exact) mass for the most abundant isotope (NIST values).
    pub fn monoisotopic_mass(self) -> f64 {
        match self {
            Self::H => 1.00782503207,
            Self::He => 4.00260325413,
            Self::Li => 7.0160034366,
            Self::Be => 9.012183065,
            Self::B => 11.00930536,
            Self::C => 12.0,
            Self::N => 14.00307400443,
            Self::O => 15.99491461957,
            Self::F => 18.99840316273,
            Self::Ne => 19.9924401762,
            Self::Na => 22.9897692820,
            Self::Mg => 23.985041697,
            Self::Al => 26.98153853,
            Self::Si => 27.97692653465,
            Self::P => 30.97376199842,
            Self::S => 31.9720711744,
            Self::Cl => 34.96885268,
            Self::Ar => 39.9623831237,
            Self::K => 38.9637064864,
            Self::Ca => 39.962590863,
            Self::Br => 78.9183376,
            Self::I => 126.9044719,
            Self::Fe => 55.9349375,
            Self::Cu => 62.9295975,
            Self::Zn => 63.9291420,
            Self::Se => 79.9165218,
        }
    }

    /// Returns true if this is not hydrogen.
    pub fn is_heavy(self) -> bool {
        self != Self::H
    }
}
