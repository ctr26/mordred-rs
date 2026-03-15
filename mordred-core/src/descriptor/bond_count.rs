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

bond_count_descriptor!(
    SingleBondCount,
    "nSingle",
    "Number of single bonds",
    BondOrder::Single
);
bond_count_descriptor!(
    DoubleBondCount,
    "nDouble",
    "Number of double bonds",
    BondOrder::Double
);
bond_count_descriptor!(
    TripleBondCount,
    "nTriple",
    "Number of triple bonds",
    BondOrder::Triple
);
bond_count_descriptor!(
    AromaticBondCount,
    "nAromatic",
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
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_double_bond_count_formaldehyde() {
        // Formaldehyde: C=O
        let mol = parse_smiles("C=O").unwrap();
        assert_eq!(DoubleBondCount.calculate(&mol).unwrap(), 1.0);
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 0.0);
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
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 2.0);
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
        assert_eq!(SingleBondCount.calculate(&mol).unwrap(), 0.0);
    }
}
