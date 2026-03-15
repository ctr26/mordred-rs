/// Connectivity indices (Randic-type chi indices).
///
/// Reference: mordred Python `_chi.py`
///
/// Chi indices use vertex degrees (including implicit hydrogens) to compute
/// graph-connectivity descriptors.
///
/// Chi0 = Σ_v  1/sqrt(delta_v)           (sum over all vertices)
/// Chi1 = Σ_e  1/sqrt(delta_u * delta_v) (sum over all edges)
///
/// where delta_v = total degree of vertex v (explicit bonds + implicit H).
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Zero-order connectivity index (Chi0).
///
/// Chi0 = Σ 1/sqrt(δᵢ) for each atom i, where δᵢ is the total degree.
pub struct Chi0;

impl Descriptor for Chi0 {
    fn name(&self) -> &str {
        "Chi0"
    }

    fn description(&self) -> &str {
        "Zero-order Randic connectivity index"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let mut sum = 0.0f64;
        for idx in mol.graph.node_indices() {
            let delta = mol.total_degree(idx);
            if delta > 0 {
                sum += 1.0 / (delta as f64).sqrt();
            }
        }
        Ok(sum)
    }
}

/// First-order connectivity index (Chi1).
///
/// Chi1 = Σ 1/sqrt(δᵢ · δⱼ) for each edge (i,j).
pub struct Chi1;

impl Descriptor for Chi1 {
    fn name(&self) -> &str {
        "Chi1"
    }

    fn description(&self) -> &str {
        "First-order Randic connectivity index"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let mut sum = 0.0f64;
        for (a, b, _) in mol.bonds() {
            let da = mol.total_degree(a);
            let db = mol.total_degree(b);
            if da > 0 && db > 0 {
                sum += 1.0 / ((da * db) as f64).sqrt();
            }
        }
        Ok(sum)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_chi0_methane() {
        // CH4: one atom with total degree 4
        // Chi0 = 1/sqrt(4) = 0.5
        let mol = parse_smiles("C").unwrap();
        let val = Chi0.calculate(&mol).unwrap();
        assert!((val - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_chi0_ethane() {
        // CC: two C atoms, each with total degree 4 (1 bond + 3 implicit H)
        // Chi0 = 2 * 1/sqrt(4) = 1.0
        let mol = parse_smiles("CC").unwrap();
        let val = Chi0.calculate(&mol).unwrap();
        assert!((val - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_chi1_ethane() {
        // CC: one edge, both atoms have total degree 4
        // Chi1 = 1/sqrt(4*4) = 1/4 = 0.25
        let mol = parse_smiles("CC").unwrap();
        let val = Chi1.calculate(&mol).unwrap();
        assert!((val - 0.25).abs() < 1e-6);
    }

    #[test]
    fn test_chi1_propane() {
        // CCC: two edges
        // C0: degree 4 (1 bond + 3H), C1: degree 4 (2 bonds + 2H), C2: degree 4 (1 bond + 3H)
        // Edge(0,1): 1/sqrt(4*4)=0.25, Edge(1,2): 1/sqrt(4*4)=0.25
        // Chi1 = 0.5
        let mol = parse_smiles("CCC").unwrap();
        let val = Chi1.calculate(&mol).unwrap();
        assert!((val - 0.5).abs() < 1e-6);
    }
}
