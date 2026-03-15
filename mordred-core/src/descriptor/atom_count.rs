/// Element-specific atom count descriptors.
///
/// Reference: mordred Python `AtomCount.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{Element, Molecule};

/// Macro to reduce boilerplate for element count descriptors.
macro_rules! element_count_descriptor {
    ($name:ident, $display:expr, $desc:expr, $element:expr) => {
        pub struct $name;

        impl Descriptor for $name {
            fn name(&self) -> &str {
                $display
            }
            fn description(&self) -> &str {
                $desc
            }
            fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
                Ok(mol.count_element($element) as f64)
            }
        }
    };
}

element_count_descriptor!(CarbonCount, "nC", "Number of carbon atoms", Element::C);
element_count_descriptor!(NitrogenCount, "nN", "Number of nitrogen atoms", Element::N);
element_count_descriptor!(OxygenCount, "nO", "Number of oxygen atoms", Element::O);
element_count_descriptor!(SulfurCount, "nS", "Number of sulfur atoms", Element::S);
element_count_descriptor!(
    PhosphorusCount,
    "nP",
    "Number of phosphorus atoms",
    Element::P
);

/// Number of halogen atoms (F, Cl, Br, I).
pub struct HalogenCount;

impl Descriptor for HalogenCount {
    fn name(&self) -> &str {
        "nX"
    }
    fn description(&self) -> &str {
        "Number of halogen atoms"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let count = mol
            .atoms()
            .filter(|(_, a)| {
                matches!(
                    a.element,
                    Element::F | Element::Cl | Element::Br | Element::I
                )
            })
            .count();
        Ok(count as f64)
    }
}

/// Number of heteroatoms (non-C, non-H).
pub struct HeteroatomCount;

impl Descriptor for HeteroatomCount {
    fn name(&self) -> &str {
        "nHetero"
    }
    fn description(&self) -> &str {
        "Number of heteroatoms"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let count = mol
            .atoms()
            .filter(|(_, a)| a.element != Element::C && a.element != Element::H)
            .count();
        Ok(count as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_carbon_count_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(CarbonCount.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_nitrogen_count_aniline() {
        let mol = parse_smiles("c1ccccc1N").unwrap();
        assert_eq!(NitrogenCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_oxygen_count_acetic_acid() {
        let mol = parse_smiles("CC(=O)O").unwrap();
        assert_eq!(OxygenCount.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_sulfur_count_dmso() {
        // DMSO: CS(=O)C
        let mol = parse_smiles("CS(=O)C").unwrap();
        assert_eq!(SulfurCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_phosphorus_count() {
        // Trimethylphosphine: CP(C)C
        let mol = parse_smiles("CP(C)C").unwrap();
        assert_eq!(PhosphorusCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_halogen_count() {
        // Chloroform: ClC(Cl)Cl
        let mol = parse_smiles("ClC(Cl)Cl").unwrap();
        assert_eq!(HalogenCount.calculate(&mol).unwrap(), 3.0);
    }

    #[test]
    fn test_halogen_count_mixed() {
        // Mixed halogens: FC(Cl)(Br)I
        let mol = parse_smiles("FC(Cl)(Br)I").unwrap();
        assert_eq!(HalogenCount.calculate(&mol).unwrap(), 4.0);
    }

    #[test]
    fn test_heteroatom_count_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(HeteroatomCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_heteroatom_count_glycine() {
        // Glycine: NCC(=O)O
        let mol = parse_smiles("NCC(=O)O").unwrap();
        assert_eq!(HeteroatomCount.calculate(&mol).unwrap(), 3.0);
    }

    #[test]
    fn test_carbon_count_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(CarbonCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_no_nitrogen() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(NitrogenCount.calculate(&mol).unwrap(), 0.0);
    }
}
