/// Rotatable-bond descriptors: count and ratio.
///
/// A bond is rotatable iff it is single, not in a ring, has both endpoints
/// with ≥ 2 heavy neighbors (terminal-bond rotation is degenerate), and is
/// not an amide C-N bond. Aromatic bonds aren't single, so they're excluded
/// by construction. Reference: upstream `mordred/RotatableBond.py` (PaDEL-
/// descriptor formulation, which follows RDKit's
/// `Lipinski.NumRotatableBonds`).
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{BondOrder, Element, Molecule};

/// Returns true if `bond (a, b)` is the C-N single bond of an amide group:
/// the C has at least one neighbor connected via a double bond to an O.
fn is_amide_bond(mol: &Molecule, a: NodeIndex, b: NodeIndex) -> bool {
    let (carbon, _nitrogen) = match (mol.atom(a).element, mol.atom(b).element) {
        (Element::C, Element::N) => (a, b),
        (Element::N, Element::C) => (b, a),
        _ => return false,
    };
    mol.graph.edges(carbon).any(|e| {
        let other = if e.source() == carbon {
            e.target()
        } else {
            e.source()
        };
        e.weight().order == BondOrder::Double && mol.atom(other).element == Element::O
    })
}

fn count_rotatable_bonds(mol: &Molecule) -> u32 {
    let ring_info = mol.ring_info();
    mol.bonds()
        .filter(|(_, _, b)| b.order == BondOrder::Single)
        .filter(|(a, b, _)| !ring_info.is_ring_bond(*a, *b))
        .filter(|(a, b, _)| {
            mol.graph.neighbors(*a).count() >= 2 && mol.graph.neighbors(*b).count() >= 2
        })
        .filter(|(a, b, _)| !is_amide_bond(mol, *a, *b))
        .count() as u32
}

pub struct RotatableBondCount;

impl Descriptor for RotatableBondCount {
    fn name(&self) -> &str {
        "nRot"
    }
    fn description(&self) -> &str {
        "Number of rotatable bonds"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(f64::from(count_rotatable_bonds(mol)))
    }
}

/// `RotRatio = nRot / nHeavyBonds`. NaN when no heavy bonds.
pub struct RotatableBondRatio;

impl Descriptor for RotatableBondRatio {
    fn name(&self) -> &str {
        "RotRatio"
    }
    fn description(&self) -> &str {
        "Rotatable bond fraction (nRot / nHeavyBonds)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let n_heavy = mol.bond_count();
        if n_heavy == 0 {
            return Ok(f64::NAN);
        }
        Ok(f64::from(count_rotatable_bonds(mol)) / n_heavy as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_hexane_three_rotatable() {
        // CCCCCC: terminal Cs have degree 1; only the middle 3 C-C bonds qualify.
        let mol = parse_smiles("CCCCCC").unwrap();
        assert_eq!(RotatableBondCount.calculate(&mol).unwrap(), 3.0);
        assert!((RotatableBondRatio.calculate(&mol).unwrap() - 0.6).abs() < 1e-10);
    }

    #[test]
    fn test_benzene_no_rotatable() {
        // All ring bonds are aromatic (or in-ring) → 0 rotatable.
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(RotatableBondCount.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_methylcyclopropane_no_rotatable() {
        // Methylcyclopropane CC1CC1: only C-C single bonds are in the ring
        // or terminal (methyl). 0 rotatable.
        let mol = parse_smiles("CC1CC1").unwrap();
        assert_eq!(RotatableBondCount.calculate(&mol).unwrap(), 0.0);
    }
}
