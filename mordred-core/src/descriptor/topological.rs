/// Topological descriptors — graph-distance based.
///
/// Reference: mordred Python `_topology.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Wiener index: sum of all shortest-path distances between pairs of atoms.
///
/// W = Σ_{i<j} d(i,j)
pub struct WienerIndex;

impl Descriptor for WienerIndex {
    fn name(&self) -> &str {
        "WienerIndex"
    }

    fn description(&self) -> &str {
        "Wiener index (sum of shortest path distances)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        if mol.atom_count() < 2 {
            return Ok(0.0);
        }
        let dist = mol.distance_matrix();
        let nodes: Vec<_> = mol.graph.node_indices().collect();
        let mut sum = 0i64;
        for i in 0..nodes.len() {
            for j in (i + 1)..nodes.len() {
                if let Some(&d) = dist.get(&(nodes[i], nodes[j])) {
                    sum += d;
                }
            }
        }
        Ok(sum as f64)
    }
}

/// First Zagreb index: M1 = Σ deg(v)^2 for all vertices v.
pub struct ZagrebIndex1;

impl Descriptor for ZagrebIndex1 {
    fn name(&self) -> &str {
        "Zagreb1"
    }

    fn description(&self) -> &str {
        "First Zagreb index (sum of squared degrees)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let sum: usize = mol
            .graph
            .node_indices()
            .map(|i| {
                let d = mol.degree(i);
                d * d
            })
            .sum();
        Ok(sum as f64)
    }
}

/// Second Zagreb index: M2 = Σ deg(u)*deg(v) for all edges (u,v).
pub struct ZagrebIndex2;

impl Descriptor for ZagrebIndex2 {
    fn name(&self) -> &str {
        "Zagreb2"
    }

    fn description(&self) -> &str {
        "Second Zagreb index (sum of degree products over edges)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let sum: usize = mol
            .bonds()
            .map(|(a, b, _)| mol.degree(a) * mol.degree(b))
            .sum();
        Ok(sum as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_wiener_ethane() {
        // C-C: distance 1, so W = 1
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(WienerIndex.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_wiener_propane() {
        // C-C-C: d(0,1)=1, d(0,2)=2, d(1,2)=1, W=4
        let mol = parse_smiles("CCC").unwrap();
        assert_eq!(WienerIndex.calculate(&mol).unwrap(), 4.0);
    }

    #[test]
    fn test_wiener_butane() {
        // C-C-C-C: d(0,1)=1, d(0,2)=2, d(0,3)=3, d(1,2)=1, d(1,3)=2, d(2,3)=1 => W=10
        let mol = parse_smiles("CCCC").unwrap();
        assert_eq!(WienerIndex.calculate(&mol).unwrap(), 10.0);
    }

    #[test]
    fn test_zagreb1_ethane() {
        // Two atoms each with degree 1: 1^2 + 1^2 = 2
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(ZagrebIndex1.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_zagreb1_propane() {
        // Degrees: 1, 2, 1 => 1 + 4 + 1 = 6
        let mol = parse_smiles("CCC").unwrap();
        assert_eq!(ZagrebIndex1.calculate(&mol).unwrap(), 6.0);
    }

    #[test]
    fn test_zagreb2_propane() {
        // Edges: (0,1) with deg 1*2=2, (1,2) with deg 2*1=2 => M2=4
        let mol = parse_smiles("CCC").unwrap();
        assert_eq!(ZagrebIndex2.calculate(&mol).unwrap(), 4.0);
    }
}
