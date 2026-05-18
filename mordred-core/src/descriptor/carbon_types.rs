/// Carbon-type counts: `C{n}SP{m}` for n carbon-neighbors and sp{m} hybridization.
///
/// 9 descriptors total — `C1SP1, C2SP1, C1SP2, C2SP2, C3SP2, C1SP3, C2SP3, C3SP3,
/// C4SP3` — plus `HybRatio`.
///
/// C atoms with **zero** carbon neighbors (e.g. caffeine's N-methyls) are
/// excluded from CarbonTypes counters but still feed into the `HybRatio`
/// denominator. Reference: upstream `mordred/CarbonTypes.py` (PaDEL-descriptor
/// formulation).
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{Atom, BondOrder, Element, Molecule};

/// Returns sp¹ (1), sp² (2), or sp³ (3) for a carbon atom based on its bond
/// environment. Aromatic atoms are sp². Triple bond → sp¹. Double bond → sp².
fn carbon_hybridization(mol: &Molecule, idx: NodeIndex, atom: &Atom) -> u8 {
    if atom.is_aromatic {
        return 2;
    }
    let mut has_double = false;
    for edge in mol.graph.edges(idx) {
        match edge.weight().order {
            BondOrder::Triple => return 1,
            BondOrder::Double => has_double = true,
            BondOrder::Aromatic => return 2,
            BondOrder::Single => {}
        }
    }
    if has_double {
        2
    } else {
        3
    }
}

fn count_carbon_neighbors(mol: &Molecule, idx: NodeIndex) -> u8 {
    mol.graph
        .neighbors(idx)
        .filter(|n| mol.atom(*n).element == Element::C)
        .count() as u8
}

/// `C{n}SP{m}`: count of carbons with exactly `n` carbon neighbors and `sp{m}`
/// hybridization. Valid (n, m) pairs: (1,1)/(2,1) for sp¹; (1,2)/(2,2)/(3,2) for
/// sp²; (1,3)/(2,3)/(3,3)/(4,3) for sp³.
pub struct CarbonTypeCount {
    n_c_neighbors: u8,
    sp_class: u8,
    name: String,
    description: String,
}

impl CarbonTypeCount {
    pub fn new(n_c_neighbors: u8, sp_class: u8) -> Self {
        let name = format!("C{n_c_neighbors}SP{sp_class}");
        let description = format!(
            "Number of sp{sp_class} carbons with {n_c_neighbors} carbon neighbors"
        );
        Self {
            n_c_neighbors,
            sp_class,
            name,
            description,
        }
    }
}

impl Descriptor for CarbonTypeCount {
    fn name(&self) -> &str {
        &self.name
    }
    fn description(&self) -> &str {
        &self.description
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let count = mol
            .atoms()
            .filter(|(_, a)| a.element == Element::C)
            .filter(|(idx, a)| carbon_hybridization(mol, *idx, a) == self.sp_class)
            .filter(|(idx, _)| count_carbon_neighbors(mol, *idx) == self.n_c_neighbors)
            .count();
        Ok(count as f64)
    }
}

/// `HybRatio = n_sp3 / (n_sp2 + n_sp3)` over all carbon atoms (including those
/// with zero carbon neighbors). Returns `NaN` when no sp²/sp³ carbons exist.
pub struct HybRatio;

impl Descriptor for HybRatio {
    fn name(&self) -> &str {
        "HybRatio"
    }
    fn description(&self) -> &str {
        "sp3 / (sp2 + sp3) carbon hybridization ratio"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let mut n_sp2 = 0u32;
        let mut n_sp3 = 0u32;
        for (idx, atom) in mol.atoms() {
            if atom.element != Element::C {
                continue;
            }
            match carbon_hybridization(mol, idx, atom) {
                2 => n_sp2 += 1,
                3 => n_sp3 += 1,
                _ => {}
            }
        }
        let denom = n_sp2 + n_sp3;
        if denom == 0 {
            Ok(f64::NAN)
        } else {
            Ok(n_sp3 as f64 / denom as f64)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_hexane_all_sp3() {
        // Hexane: 6 sp³ C's. 2 terminal (1 C neighbor) + 4 internal (2 C neighbors).
        let mol = parse_smiles("CCCCCC").unwrap();
        assert_eq!(CarbonTypeCount::new(1, 3).calculate(&mol).unwrap(), 2.0);
        assert_eq!(CarbonTypeCount::new(2, 3).calculate(&mol).unwrap(), 4.0);
        assert_eq!(HybRatio.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_benzene_all_sp2() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        // Each aromatic C has 2 ring-C neighbors → C2SP2.
        assert_eq!(CarbonTypeCount::new(2, 2).calculate(&mol).unwrap(), 6.0);
        assert_eq!(HybRatio.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_methyl_attached_to_n_excluded() {
        // Methylamine CH3-NH2: 1 C with 0 C neighbors → not counted in CarbonTypes,
        // but still in HybRatio denominator.
        let mol = parse_smiles("CN").unwrap();
        assert_eq!(CarbonTypeCount::new(1, 3).calculate(&mol).unwrap(), 0.0);
        assert_eq!(HybRatio.calculate(&mol).unwrap(), 1.0); // 1 sp3 / (1 sp3 + 0 sp2)
    }

    #[test]
    fn test_allicin_mixed_hybridization() {
        // Allicin C=CCSS(=O)CC=C: 4 sp² + 2 sp³ Cs. HybRatio = 2/6 ≈ 0.333.
        let mol = parse_smiles("C=CCSS(=O)CC=C").unwrap();
        let ratio = HybRatio.calculate(&mol).unwrap();
        assert!((ratio - 1.0 / 3.0).abs() < 1e-6, "got {ratio}");
    }

    #[test]
    fn test_acetylene_sp1() {
        // HC≡CH: 2 sp¹ carbons, each with 1 C neighbor.
        let mol = parse_smiles("C#C").unwrap();
        assert_eq!(CarbonTypeCount::new(1, 1).calculate(&mol).unwrap(), 2.0);
    }
}
