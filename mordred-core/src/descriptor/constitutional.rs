/// Constitutional descriptors — basic molecular composition.
///
/// Reference: mordred Python `_constitutional.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Total number of atoms including implicit hydrogens.
pub struct AtomCount;

impl Descriptor for AtomCount {
    fn name(&self) -> &str {
        "nAtom"
    }

    fn description(&self) -> &str {
        "Total atom count (including implicit H)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.total_atom_count() as f64)
    }
}

/// Number of heavy (non-hydrogen) atoms.
pub struct HeavyAtomCount;

impl Descriptor for HeavyAtomCount {
    fn name(&self) -> &str {
        "nHeavyAtom"
    }

    fn description(&self) -> &str {
        "Heavy atom count (non-hydrogen)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.heavy_atom_count() as f64)
    }
}

/// Total number of bonds including implicit hydrogen bonds.
pub struct BondCount;

impl Descriptor for BondCount {
    fn name(&self) -> &str {
        "nBonds"
    }

    fn description(&self) -> &str {
        "Bond count (including implicit H bonds)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let implicit_h_bonds: usize = mol.atoms().map(|(_, a)| a.implicit_h as usize).sum();
        Ok((mol.bond_count() + implicit_h_bonds) as f64)
    }
}

/// Molecular weight (monoisotopic exact mass, including implicit H).
pub struct MolecularWeight;

impl Descriptor for MolecularWeight {
    fn name(&self) -> &str {
        "MW"
    }

    fn description(&self) -> &str {
        "Molecular weight (exact mass, monoisotopic)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.molecular_weight())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_methane_atom_count() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(AtomCount.calculate(&mol).unwrap(), 5.0); // C + 4H
    }

    #[test]
    fn test_methane_heavy_atom_count() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_ethanol_bond_count() {
        let mol = parse_smiles("CCO").unwrap();
        // 2 heavy-atom bonds + 6 implicit H bonds (3+2+1) = 8
        assert_eq!(BondCount.calculate(&mol).unwrap(), 8.0);
    }

    #[test]
    fn test_methane_mw() {
        let mol = parse_smiles("C").unwrap();
        let mw = MolecularWeight.calculate(&mol).unwrap();
        // CH4 = 12.0 + 4*1.00783 ≈ 16.031
        assert!((mw - 16.031).abs() < 0.01);
    }

    #[test]
    fn test_water_mw() {
        let mol = parse_smiles("O").unwrap();
        let mw = MolecularWeight.calculate(&mol).unwrap();
        // H2O = 15.995 + 2*1.00783 ≈ 18.011
        assert!((mw - 18.011).abs() < 0.01);
    }
}
