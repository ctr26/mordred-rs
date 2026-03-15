use serde::{Deserialize, Serialize};

use super::element::Element;

/// An atom in a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    pub element: Element,
    pub charge: i8,
    pub is_aromatic: bool,
    pub implicit_h: u8,
    pub isotope: Option<u16>,
}

impl Atom {
    pub fn new(element: Element) -> Self {
        Self {
            element,
            charge: 0,
            is_aromatic: false,
            implicit_h: 0,
            isotope: None,
        }
    }

    /// Total degree = number of explicit bonds + implicit hydrogens.
    pub fn total_degree(&self, explicit_bonds: u8) -> u8 {
        explicit_bonds + self.implicit_h
    }

    /// Atomic weight including implicit hydrogens.
    pub fn mass(&self) -> f64 {
        self.element.atomic_weight() + f64::from(self.implicit_h) * Element::H.atomic_weight()
    }
}
