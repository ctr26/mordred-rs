use serde::{Deserialize, Serialize};

use super::element::Element;

/// An atom in a molecule.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Atom {
    /// Chemical element of this atom.
    pub element: Element,
    /// Formal charge (e.g. +1, -1).
    pub charge: i8,
    /// Whether this atom is part of an aromatic system.
    pub is_aromatic: bool,
    /// Number of implicit hydrogen atoms.
    pub implicit_h: u8,
    /// Isotope mass number, if specified.
    pub isotope: Option<u16>,
}

impl Atom {
    /// Creates a new atom with the given element and default properties.
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
