//! Molecule representation and SMILES parsing.
//!
//! This module provides the core data structures for representing molecules
//! as undirected graphs of [`Atom`] nodes and [`Bond`] edges, built on
//! [`petgraph`]. The [`smiles`] submodule converts SMILES strings into
//! [`Molecule`] instances.

pub mod aromaticity;
pub mod atom;
pub mod bond;
pub mod element;
pub mod mol;
pub mod rings;
pub mod smiles;

pub use aromaticity::perceive_aromaticity;
pub use atom::Atom;
pub use bond::{Bond, BondOrder};
pub use element::Element;
pub use mol::{Molecule, MolecularProperties};
pub use rings::RingInfo;
pub use smiles::parse_smiles;
