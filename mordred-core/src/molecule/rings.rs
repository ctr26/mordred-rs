/// Ring detection for molecular graphs.
///
/// Implements SSSR (Smallest Set of Smallest Rings) detection using
/// a BFS-based approach suitable for molecular graphs.
use petgraph::graph::{NodeIndex, UnGraph};
use petgraph::visit::NodeIndexable;
use std::collections::{BTreeSet, HashSet, VecDeque};

use super::atom::Atom;
use super::bond::Bond;

/// Information about rings in a molecule.
#[derive(Debug, Clone)]
pub struct RingInfo {
    /// The SSSR rings, each as a vector of atom indices.
    rings: Vec<Vec<NodeIndex>>,
    /// Set of atoms that are in at least one ring.
    ring_atoms: HashSet<NodeIndex>,
    /// Set of bonds (as sorted pairs) that are in at least one ring.
    ring_bonds: HashSet<(NodeIndex, NodeIndex)>,
}

impl RingInfo {
    /// Detect rings in a molecular graph.
    pub fn detect(graph: &UnGraph<Atom, Bond>) -> Self {
        let rings = find_sssr(graph);
        let mut ring_atoms = HashSet::new();
        let mut ring_bonds = HashSet::new();

        for ring in &rings {
            for &atom in ring {
                ring_atoms.insert(atom);
            }
            for i in 0..ring.len() {
                let a = ring[i];
                let b = ring[(i + 1) % ring.len()];
                let bond = if a < b { (a, b) } else { (b, a) };
                ring_bonds.insert(bond);
            }
        }

        Self {
            rings,
            ring_atoms,
            ring_bonds,
        }
    }

    /// Get all rings.
    pub fn rings(&self) -> &[Vec<NodeIndex>] {
        &self.rings
    }

    /// Number of rings in the SSSR.
    pub fn num_rings(&self) -> usize {
        self.rings.len()
    }

    /// Check if an atom is in any ring.
    pub fn is_in_ring(&self, atom: NodeIndex) -> bool {
        self.ring_atoms.contains(&atom)
    }

    /// Check if a bond is in any ring.
    pub fn is_ring_bond(&self, a: NodeIndex, b: NodeIndex) -> bool {
        let bond = if a < b { (a, b) } else { (b, a) };
        self.ring_bonds.contains(&bond)
    }

    /// Get the smallest ring size containing an atom (None if not in a ring).
    pub fn smallest_ring_size(&self, atom: NodeIndex) -> Option<usize> {
        self.rings
            .iter()
            .filter(|ring| ring.contains(&atom))
            .map(|ring| ring.len())
            .min()
    }

    /// Get rings of a specific size.
    pub fn rings_of_size(&self, size: usize) -> Vec<&Vec<NodeIndex>> {
        self.rings.iter().filter(|r| r.len() == size).collect()
    }
}

/// Find the SSSR using a BFS-based approach.
fn find_sssr(graph: &UnGraph<Atom, Bond>) -> Vec<Vec<NodeIndex>> {
    let n_vertices = graph.node_count();
    let n_edges = graph.edge_count();
    if n_vertices == 0 || n_edges == 0 {
        return Vec::new();
    }

    // Circuit rank = |E| - |V| + number of connected components
    let n_components = count_components(graph);
    // Use saturating arithmetic to avoid underflow for acyclic graphs
    let circuit_rank = (n_edges + n_components).saturating_sub(n_vertices);

    if circuit_rank == 0 {
        return Vec::new();
    }

    // Find candidate rings using BFS from each node
    let mut all_rings: Vec<Vec<NodeIndex>> = Vec::new();

    let nodes: Vec<NodeIndex> = graph.node_indices().collect();

    for &start in &nodes {
        let rings = find_rings_from_node(graph, start);
        for ring in rings {
            if !is_duplicate_ring(&all_rings, &ring) {
                all_rings.push(ring);
            }
        }
    }

    // Sort by size to prefer smallest rings
    all_rings.sort_by_key(|r| r.len());

    // Select linearly independent rings up to the circuit rank
    let mut sssr = Vec::new();
    for ring in all_rings {
        if sssr.len() >= circuit_rank {
            break;
        }
        if is_independent(&sssr, &ring, graph) {
            sssr.push(ring);
        }
    }

    sssr
}

/// Find rings starting from a given node using BFS.
///
/// When two shortest paths from the start node meet at a node, the
/// combination of those paths forms a ring candidate.
fn find_rings_from_node(graph: &UnGraph<Atom, Bond>, start: NodeIndex) -> Vec<Vec<NodeIndex>> {
    let n = graph.node_bound();
    let mut dist: Vec<Option<usize>> = vec![None; n];
    let mut parent: Vec<Vec<NodeIndex>> = vec![Vec::new(); n];
    let mut rings = Vec::new();

    dist[start.index()] = Some(0);
    let mut queue = VecDeque::new();
    queue.push_back(start);

    while let Some(current) = queue.pop_front() {
        let current_dist = dist[current.index()].unwrap();

        for neighbor in graph.neighbors(current) {
            let ni = neighbor.index();
            match dist[ni] {
                None => {
                    // Unvisited: set distance and parent
                    dist[ni] = Some(current_dist + 1);
                    parent[ni] = vec![current];
                    queue.push_back(neighbor);
                }
                Some(d) if d == current_dist + 1 => {
                    // Same level: another shortest path parent
                    parent[ni].push(current);
                }
                Some(d) if d == current_dist => {
                    // Same distance as current — odd-length ring
                    // Build ring from the two paths back to start
                    let path1 = trace_path(&parent, start, current);
                    let path2 = trace_path(&parent, start, neighbor);
                    if let Some(ring) = build_ring(&path1, &path2) {
                        if ring.len() >= 3 {
                            rings.push(ring);
                        }
                    }
                }
                _ => {
                    // Already visited at shorter distance — back edge
                    // This creates an even-length ring
                }
            }
        }
    }

    // Also find even-length rings: when a node has multiple parents at the
    // same distance, the two paths form a ring.
    for ni in 0..n {
        if parent[ni].len() > 1 {
            let node = NodeIndex::new(ni);
            for i in 0..parent[ni].len() {
                for j in (i + 1)..parent[ni].len() {
                    let path1 = trace_path(&parent, start, parent[ni][i]);
                    let path2 = trace_path(&parent, start, parent[ni][j]);
                    // Build ring: path1 -> node -> reverse(path2)
                    let mut ring_path = path1;
                    ring_path.push(node);
                    let mut p2 = path2;
                    p2.reverse();
                    // Remove the start node from p2 since path1 already starts with it
                    if !p2.is_empty() {
                        p2.pop(); // remove start from the reversed path
                    }
                    ring_path.extend(p2);
                    if ring_path.len() >= 3 {
                        rings.push(ring_path);
                    }
                }
            }
        }
    }

    rings
}

/// Trace the shortest path from start to target using the parent array.
fn trace_path(parent: &[Vec<NodeIndex>], start: NodeIndex, target: NodeIndex) -> Vec<NodeIndex> {
    let mut path = vec![target];
    let mut current = target;
    while current != start {
        if parent[current.index()].is_empty() {
            return path; // shouldn't happen in a connected component
        }
        current = parent[current.index()][0];
        path.push(current);
    }
    path.reverse();
    path
}

/// Build a ring from two paths that share a common start.
fn build_ring(path1: &[NodeIndex], path2: &[NodeIndex]) -> Option<Vec<NodeIndex>> {
    if path1.is_empty() || path2.is_empty() {
        return None;
    }
    // Both paths start at the same node.
    // Ring = path1[1..] + reverse(path2[1..])
    // (skip the shared start, which closes the ring implicitly)
    let mut ring = Vec::new();
    ring.extend_from_slice(path1);
    let mut p2_rev: Vec<NodeIndex> = path2[1..].to_vec();
    p2_rev.reverse();
    ring.extend(p2_rev);

    // Deduplicate: if path1 and path2 share intermediate nodes, skip
    let unique: HashSet<NodeIndex> = ring.iter().copied().collect();
    if unique.len() != ring.len() {
        return None;
    }

    Some(ring)
}

/// Check if a ring (as a set of atoms) already exists in the list.
fn is_duplicate_ring(existing: &[Vec<NodeIndex>], candidate: &[NodeIndex]) -> bool {
    let candidate_set: BTreeSet<NodeIndex> = candidate.iter().copied().collect();
    existing.iter().any(|r| {
        let r_set: BTreeSet<NodeIndex> = r.iter().copied().collect();
        r_set == candidate_set
    })
}

/// Convert a ring to its edge set representation (as sorted pairs).
fn ring_to_edge_set(ring: &[NodeIndex], graph: &UnGraph<Atom, Bond>) -> HashSet<(usize, usize)> {
    let mut edges = HashSet::new();
    for i in 0..ring.len() {
        let a = ring[i].index();
        let b = ring[(i + 1) % ring.len()].index();
        // Only include if this edge actually exists in the graph
        let na = NodeIndex::new(a);
        let nb = NodeIndex::new(b);
        if graph.find_edge(na, nb).is_some() {
            let edge = if a < b { (a, b) } else { (b, a) };
            edges.insert(edge);
        }
    }
    edges
}

/// Check if a ring is linearly independent from existing SSSR rings.
///
/// A ring is independent if its edge set cannot be expressed as the
/// symmetric difference (XOR) of any combination of existing ring edge sets.
/// For efficiency, we check if the candidate has at least one edge not
/// present in the union of existing rings, or if it introduces a new
/// combination.
fn is_independent(
    sssr: &[Vec<NodeIndex>],
    candidate: &[NodeIndex],
    graph: &UnGraph<Atom, Bond>,
) -> bool {
    if sssr.is_empty() {
        return true;
    }

    let candidate_edges = ring_to_edge_set(candidate, graph);
    if candidate_edges.is_empty() {
        return false;
    }

    // Use Gaussian elimination over GF(2) on edge vectors.
    // Collect all unique edges across all rings + candidate for indexing.
    let mut all_edges_set = BTreeSet::new();
    let mut edge_sets: Vec<HashSet<(usize, usize)>> = Vec::new();
    for r in sssr {
        let es = ring_to_edge_set(r, graph);
        for e in &es {
            all_edges_set.insert(*e);
        }
        edge_sets.push(es);
    }
    for e in &candidate_edges {
        all_edges_set.insert(*e);
    }

    let all_edges: Vec<(usize, usize)> = all_edges_set.into_iter().collect();
    let edge_to_idx: std::collections::HashMap<(usize, usize), usize> =
        all_edges.iter().enumerate().map(|(i, &e)| (e, i)).collect();
    let n_edges = all_edges.len();

    // Build GF(2) basis from existing SSSR rings
    let mut basis: Vec<Vec<bool>> = Vec::new();
    for es in &edge_sets {
        let mut row = vec![false; n_edges];
        for e in es {
            if let Some(&idx) = edge_to_idx.get(e) {
                row[idx] = true;
            }
        }
        reduce_and_insert(&mut basis, row, n_edges);
    }

    // Try to reduce the candidate row against the basis
    let mut row = vec![false; n_edges];
    for e in &candidate_edges {
        if let Some(&idx) = edge_to_idx.get(e) {
            row[idx] = true;
        }
    }
    for basis_row in &basis {
        if let Some(pivot) = basis_row.iter().position(|&b| b) {
            if row[pivot] {
                for j in 0..n_edges {
                    row[j] ^= basis_row[j];
                }
            }
        }
    }

    // If row is not all zeros, the candidate is independent
    row.iter().any(|&b| b)
}

/// Insert a row into the GF(2) basis using Gaussian elimination.
fn reduce_and_insert(basis: &mut Vec<Vec<bool>>, mut row: Vec<bool>, n: usize) {
    for basis_row in basis.iter() {
        if let Some(pivot) = basis_row.iter().position(|&b| b) {
            if row[pivot] {
                for j in 0..n {
                    row[j] ^= basis_row[j];
                }
            }
        }
    }
    // If row is non-zero, add to basis
    if row.iter().any(|&b| b) {
        basis.push(row);
    }
}

/// Count connected components in the graph using BFS.
fn count_components(graph: &UnGraph<Atom, Bond>) -> usize {
    let mut visited = HashSet::new();
    let mut components = 0;

    for node in graph.node_indices() {
        if visited.contains(&node) {
            continue;
        }
        components += 1;
        let mut queue = VecDeque::new();
        queue.push_back(node);
        visited.insert(node);
        while let Some(current) = queue.pop_front() {
            for neighbor in graph.neighbors(current) {
                if visited.insert(neighbor) {
                    queue.push_back(neighbor);
                }
            }
        }
    }

    components
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_benzene_one_ring() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 1);
        assert_eq!(info.rings()[0].len(), 6);
    }

    #[test]
    fn test_naphthalene_two_rings() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 2);
        for ring in info.rings() {
            assert_eq!(ring.len(), 6);
        }
    }

    #[test]
    fn test_cyclohexane_one_ring() {
        let mol = parse_smiles("C1CCCCC1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 1);
        assert_eq!(info.rings()[0].len(), 6);
    }

    #[test]
    fn test_cyclopropane_one_ring_size_3() {
        let mol = parse_smiles("C1CC1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 1);
        assert_eq!(info.rings()[0].len(), 3);
    }

    #[test]
    fn test_ethane_no_rings() {
        let mol = parse_smiles("CC").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 0);
    }

    #[test]
    fn test_methane_no_rings() {
        let mol = parse_smiles("C").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.num_rings(), 0);
    }

    #[test]
    fn test_is_in_ring() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        for node in mol.graph.node_indices() {
            assert!(info.is_in_ring(node));
        }
    }

    #[test]
    fn test_not_in_ring() {
        let mol = parse_smiles("CC").unwrap();
        let info = RingInfo::detect(&mol.graph);
        for node in mol.graph.node_indices() {
            assert!(!info.is_in_ring(node));
        }
    }

    #[test]
    fn test_smallest_ring_size() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        let node = mol.graph.node_indices().next().unwrap();
        assert_eq!(info.smallest_ring_size(node), Some(6));
    }

    #[test]
    fn test_rings_of_size() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let info = RingInfo::detect(&mol.graph);
        assert_eq!(info.rings_of_size(6).len(), 2);
        assert_eq!(info.rings_of_size(5).len(), 0);
    }
}
