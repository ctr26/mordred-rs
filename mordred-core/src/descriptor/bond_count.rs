/// Bond type count descriptors.
///
/// Reference: mordred Python `BondCount.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{BondOrder, Molecule};

macro_rules! bond_count_descriptor {
    ($name:ident, $display:expr, $desc:expr, $order:expr) => {
        pub struct $name;

        impl Descriptor for $name {
            fn name(&self) -> &str {
                $display
            }
            fn description(&self) -> &str {
                $desc
            }
            fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
                let count = mol.bonds().filter(|(_, _, b)| b.order == $order).count();
                Ok(count as f64)
            }
        }
    };
}

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
        let explicit_single = mol.bonds().filter(|(_, _, b)| b.order == BondOrder::Single).count();
        let implicit_h_bonds: usize = mol.atoms().map(|(_, a)| a.implicit_h as usize).sum();
        Ok((explicit_single + implicit_h_bonds) as f64)
    }
}

bond_count_descriptor!(
    DoubleBondCount,
    "nBondsD",
    "Number of double bonds",
    BondOrder::Double
);
bond_count_descriptor!(
    TripleBondCount,
    "nBondsT",
    "Number of triple bonds",
    BondOrder::Triple
);
bond_count_descriptor!(
    AromaticBondCount,
    "nBondsA",
    "Number of aromatic bonds",
    BondOrder::Aromatic
);

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
