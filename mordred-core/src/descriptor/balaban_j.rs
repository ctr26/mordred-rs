/// Balaban's J topological index.
///
/// Reference: A. T. Balaban, "Highly discriminating distance-based topological
/// index", Chem. Phys. Lett. 89 (1982) 399-404.
///
/// Upstream Python wrapper: `mordred/BalabanJ.py` (delegates to
/// `rdkit.Chem.GraphDescriptors.BalabanJ`).
///
/// Definition (heavy-atom graph only; no explicit hydrogens):
///
/// ```text
///     J = M / (mu + 1) * Σ_{(i,j) ∈ E} 1 / sqrt(s_i · s_j)
/// ```
///
/// where
///   * `M`  = number of edges (bonds) in the heavy-atom graph,
///   * `N`  = number of heavy atoms,
///   * `mu` = cyclomatic number = `M - N + 1` (for a connected graph),
///   * `E`  = set of heavy-atom bonds,
///   * `s_i` = row sum of the topological distance matrix, i.e.
///     `Σ_j d(i, j)` over all heavy atoms `j != i`.
///
/// Edge cases:
///   * For graphs with fewer than 2 atoms or no bonds, J is defined to be 0.
///   * For disconnected molecules the upstream RDKit implementation
///     restricts to one connected component; here we follow the same
///     convention by skipping atom pairs whose shortest-path distance is
///     not represented in the cached distance matrix.
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;
use petgraph::graph::NodeIndex;
use std::collections::HashMap;

/// Balaban's J index (graph-distance based connectivity descriptor).
pub struct BalabanJ;

impl Descriptor for BalabanJ {
    fn name(&self) -> &str {
        "BalabanJ"
    }

    fn description(&self) -> &str {
        "Balaban's J index"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let n_atoms = mol.atom_count();
        let n_bonds = mol.bond_count();
        if n_atoms < 2 || n_bonds == 0 {
            return Ok(0.0);
        }

        let dist: &HashMap<(NodeIndex, NodeIndex), i64> = mol.distance_matrix();
        let nodes: Vec<NodeIndex> = mol.graph.node_indices().collect();

        // Row sums s_i = Σ_j d(i, j), j != i. Floyd-Warshall stores both
        // ordered pairs (a, b) and (b, a) plus the diagonal (a, a) = 0,
        // so we can sum every (i, j) entry directly.
        let mut row_sums: HashMap<NodeIndex, f64> = HashMap::with_capacity(nodes.len());
        for &i in &nodes {
            let mut s: i64 = 0;
            for &j in &nodes {
                if i == j {
                    continue;
                }
                if let Some(&d) = dist.get(&(i, j)) {
                    s += d;
                }
            }
            row_sums.insert(i, s as f64);
        }

        // Σ_{(i,j) ∈ E} 1 / sqrt(s_i · s_j)
        let mut edge_sum = 0.0f64;
        for (a, b, _) in mol.bonds() {
            let sa = *row_sums.get(&a).unwrap_or(&0.0);
            let sb = *row_sums.get(&b).unwrap_or(&0.0);
            let prod = sa * sb;
            if prod > 0.0 {
                edge_sum += 1.0 / prod.sqrt();
            }
        }

        // mu = M - N + 1 for a connected graph; matches upstream/RDKit.
        let m = n_bonds as i64;
        let n = n_atoms as i64;
        let mu = m - n + 1;
        let denom = (mu + 1) as f64;
        if denom == 0.0 {
            return Ok(0.0);
        }
        let j = (m as f64 / denom) * edge_sum;
        Ok(j)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn test_balaban_j_ethane() {
        // CC: N=2, M=1, mu=0, M/(mu+1)=1
        // d(0,1)=1 -> s_0=s_1=1; edge sum = 1/sqrt(1*1) = 1
        // J = 1 * 1 = 1.0
        let mol = parse_smiles("CC").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        assert!(approx(v, 1.0, 1e-9), "ethane J = {}", v);
    }

    #[test]
    fn test_balaban_j_propane() {
        // CCC: N=3, M=2, mu=0, M/(mu+1)=2
        // s_0 = 1+2 = 3, s_1 = 1+1 = 2, s_2 = 2+1 = 3
        // edges: (0,1) 1/sqrt(6), (1,2) 1/sqrt(6) -> sum = 2/sqrt(6)
        // J = 2 * 2/sqrt(6) = 4/sqrt(6) ≈ 1.6329932
        let mol = parse_smiles("CCC").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        let expected = 4.0 / 6.0f64.sqrt();
        assert!(
            approx(v, expected, 1e-9),
            "propane J = {} expected {}",
            v,
            expected
        );
    }

    #[test]
    fn test_balaban_j_butane() {
        // CCCC: N=4, M=3, mu=0, M/(mu+1)=3
        // distances:
        //   d(0,1)=1, d(0,2)=2, d(0,3)=3 -> s_0=6
        //   d(1,0)=1, d(1,2)=1, d(1,3)=2 -> s_1=4
        //   s_2 = 4 (by symmetry), s_3 = 6
        // edges: (0,1) 1/sqrt(24), (1,2) 1/sqrt(16)=1/4, (2,3) 1/sqrt(24)
        //   sum = 2/sqrt(24) + 1/4
        // J = 3 * (2/sqrt(24) + 1/4)
        let mol = parse_smiles("CCCC").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        let expected = 3.0 * (2.0 / 24.0f64.sqrt() + 0.25);
        assert!(
            approx(v, expected, 1e-9),
            "butane J = {} expected {}",
            v,
            expected
        );
    }

    #[test]
    fn test_balaban_j_hexane() {
        // CCCCCC: N=6, M=5, mu=0, M/(mu+1)=5
        // s_0=1+2+3+4+5=15, s_1=1+1+2+3+4=11, s_2=2+1+1+2+3=9,
        // s_3=9, s_4=11, s_5=15
        // edges: (0,1) 1/sqrt(15*11), (1,2) 1/sqrt(11*9), (2,3) 1/sqrt(9*9)=1/9,
        //         (3,4) 1/sqrt(9*11), (4,5) 1/sqrt(11*15)
        // sum = 2/sqrt(165) + 2/sqrt(99) + 1/9
        // J ≈ 5 * 0.46778 ≈ 2.33889  (Balaban's classic hexane value)
        let mol = parse_smiles("CCCCCC").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        let s = 2.0 / 165.0f64.sqrt() + 2.0 / 99.0f64.sqrt() + 1.0 / 9.0;
        let expected = 5.0 * s;
        assert!(
            approx(v, expected, 1e-9),
            "hexane J = {} expected {}",
            v,
            expected
        );
        // Sanity: textbook value ≈ 2.339
        assert!(approx(v, 2.3389, 1e-3), "hexane J = {} (textbook 2.339)", v);
    }

    #[test]
    fn test_balaban_j_benzene() {
        // c1ccccc1: N=6, M=6, mu=1, M/(mu+1)=3
        // s_i = 1+1+2+2+3 = 9 for every atom (regular ring)
        // edge sum = 6 * 1/sqrt(81) = 6/9 = 2/3
        // J = 3 * 2/3 = 2.0
        let mol = parse_smiles("c1ccccc1").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        assert!(approx(v, 2.0, 1e-9), "benzene J = {}", v);
    }

    #[test]
    fn test_balaban_j_cyclohexane() {
        // C1CCCCC1: same heavy-atom graph as benzene
        let mol = parse_smiles("C1CCCCC1").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        assert!(approx(v, 2.0, 1e-9), "cyclohexane J = {}", v);
    }

    #[test]
    fn test_balaban_j_single_atom() {
        // Single atom: defined to be 0.
        let mol = parse_smiles("C").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        assert_eq!(v, 0.0);
    }

    #[test]
    fn test_balaban_j_isobutane() {
        // CC(C)C: N=4, M=3, mu=0, M/(mu+1)=3
        // central C connected to three terminal C; star graph K_{1,3}.
        // s_center = 1+1+1 = 3, s_term = 1+2+2 = 5 (each terminal)
        // 3 edges, each 1/sqrt(3*5) = 1/sqrt(15)
        // sum = 3/sqrt(15);  J = 3 * 3/sqrt(15) = 9/sqrt(15)
        let mol = parse_smiles("CC(C)C").unwrap();
        let v = BalabanJ.calculate(&mol).unwrap();
        let expected = 9.0 / 15.0f64.sqrt();
        assert!(
            approx(v, expected, 1e-9),
            "isobutane J = {} expected {}",
            v,
            expected
        );
    }
}
