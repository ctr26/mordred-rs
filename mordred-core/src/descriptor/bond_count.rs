/// Bond type count descriptors.
///
/// Reference: mordred Python `BondCount.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Number of single bonds including implicit hydrogen bonds.
pub struct SingleBondCount;

impl Descriptor for SingleBondCount {
    fn name(&self) -> &str {
        "nBondsS"
    }
    fn description(&self) -> &str {
        "Number of single bonds (including implicit H)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let props = mol.properties();
        Ok((props.single_bond_count + props.implicit_h_sum) as f64)
    }
}

/// Number of double bonds.
pub struct DoubleBondCount;

impl Descriptor for DoubleBondCount {
    fn name(&self) -> &str {
        "nBondsD"
    }
    fn description(&self) -> &str {
        "Number of double bonds"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.properties().double_bond_count as f64)
    }
}

/// Number of triple bonds.
pub struct TripleBondCount;

impl Descriptor for TripleBondCount {
    fn name(&self) -> &str {
        "nBondsT"
    }
    fn description(&self) -> &str {
        "Number of triple bonds"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.properties().triple_bond_count as f64)
    }
}

/// Number of aromatic bonds.
pub struct AromaticBondCount;

impl Descriptor for AromaticBondCount {
    fn name(&self) -> &str {
        "nBondsA"
    }
    fn description(&self) -> &str {
        "Number of aromatic bonds"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.properties().aromatic_bond_count as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_single_bond_count_ethane() {
        let mol = parse_smiles("CC").unwrap();
        // 1 explicit single + 6 implicit H = 7
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 7.0);
    }

    #[test]
    fn test_double_bond_count_formaldehyde() {
        // Formaldehyde: C=O (CH2O: C has 2 implicit H)
        let mol = parse_smiles("C=O").unwrap();
        assert_eq!(DoubleBondCount.calculate(&mol).unwrap(), 1.0);
        // 0 explicit single + 2 implicit H = 2
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_triple_bond_count_acetylene() {
        // Acetylene: C#C
        let mol = parse_smiles("C#C").unwrap();
        assert_eq!(TripleBondCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_aromatic_bond_count_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(AromaticBondCount.calculate(&mol).unwrap(), 6.0);
    }

    #[test]
    fn test_mixed_bonds_acetic_acid() {
        // Acetic acid: CC(=O)O
        let mol = parse_smiles("CC(=O)O").unwrap();
        // 2 explicit single + 4 implicit H (3+0+0+1) = 6
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 6.0);
        assert_eq!(DoubleBondCount.calculate(&mol).unwrap(), 1.0);
        assert_eq!(TripleBondCount.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_no_aromatic_bonds_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(AromaticBondCount.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_triple_bond_hydrogen_cyanide() {
        // HCN: C#N
        let mol = parse_smiles("C#N").unwrap();
        assert_eq!(TripleBondCount.calculate(&mol).unwrap(), 1.0);
        // 0 explicit single + 1 implicit H (H on C) = 1
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 1.0);
    }
}
