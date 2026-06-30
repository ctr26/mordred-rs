/// Vertex adjacency descriptor: log₂ of the magnitude of the vertex-adjacency matrix.
///
/// `VAdjMat = log₂(2 · nHeavyBonds)`. Returns NaN for molecules with no
/// explicit (heavy-atom) bonds. Reference: upstream `mordred/VertexAdjacencyInformation.py`.
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

pub struct VertexAdjacencyMagnitude;

impl Descriptor for VertexAdjacencyMagnitude {
    fn name(&self) -> &str {
        "VAdjMat"
    }
    fn description(&self) -> &str {
        "log2(2 × number of heavy-atom bonds)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let n = mol.bond_count();
        if n == 0 {
            return Ok(f64::NAN);
        }
        Ok((2.0 * n as f64).log2())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_hexane_vadjmat() {
        // 5 heavy bonds → log2(10) ≈ 3.321928.
        let mol = parse_smiles("CCCCCC").unwrap();
        let v = VertexAdjacencyMagnitude.calculate(&mol).unwrap();
        assert!((v - 3.321928).abs() < 1e-5);
    }

    #[test]
    fn test_benzene_vadjmat() {
        // 6 heavy bonds → log2(12) ≈ 3.5849626.
        let mol = parse_smiles("c1ccccc1").unwrap();
        let v = VertexAdjacencyMagnitude.calculate(&mol).unwrap();
        assert!((v - 3.5849626).abs() < 1e-5);
    }

    #[test]
    fn test_acetonitrile_vadjmat() {
        // CH3C#N: 2 heavy bonds (C-C single + C≡N triple) → log2(4) = 2.
        let mol = parse_smiles("CC#N").unwrap();
        let v = VertexAdjacencyMagnitude.calculate(&mol).unwrap();
        assert!((v - 2.0).abs() < 1e-10);
    }
}
