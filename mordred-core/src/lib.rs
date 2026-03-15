//! mordred-core: molecular descriptor calculator.
//!
//! This crate provides a pure-Rust implementation of the
//! [mordred](https://github.com/mordred-descriptor/mordred) molecular
//! descriptor library. It includes a SMILES parser, molecule graph
//! representation, and a growing set of 2D descriptor calculators.
//!
//! # Quick start
//!
//! ```
//! use mordred_core::{parse_smiles, DescriptorSet};
//!
//! let mol = parse_smiles("c1ccccc1").unwrap();
//! let set = DescriptorSet::all();
//! for (name, value) in set.calculate(&mol) {
//!     println!("{}: {:?}", name, value);
//! }
//! ```

pub mod descriptor;
pub mod error;
pub mod molecule;

pub use descriptor::{Descriptor, DescriptorSet};
pub use error::MordredError;
pub use molecule::{Molecule, parse_smiles};
