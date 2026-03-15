use serde::{Deserialize, Serialize};

/// Bond order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

impl BondOrder {
    /// Numeric bond order value.
    pub fn as_f64(self) -> f64 {
        match self {
            Self::Single => 1.0,
            Self::Double => 2.0,
            Self::Triple => 3.0,
            Self::Aromatic => 1.5,
        }
    }

    /// Integer contribution to valence.
    pub fn valence_contribution(self) -> u8 {
        match self {
            Self::Single | Self::Aromatic => 1,
            Self::Double => 2,
            Self::Triple => 3,
        }
    }
}

/// A bond between two atoms.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Bond {
    /// Bond order (single, double, triple, or aromatic).
    pub order: BondOrder,
    /// Whether this bond is part of an aromatic system.
    pub is_aromatic: bool,
    /// Whether this bond belongs to a ring.
    pub is_ring: bool,
}

impl Bond {
    /// Creates a new bond with the given order.
    pub fn new(order: BondOrder) -> Self {
        Self {
            is_aromatic: order == BondOrder::Aromatic,
            order,
            is_ring: false,
        }
    }
}
