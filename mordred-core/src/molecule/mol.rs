use std::cell::OnceCell;

use petgraph::algo::floyd_warshall;
use petgraph::graph::{NodeIndex, UnGraph};

use super::atom::Atom;
use super::bond::Bond;
use super::element::Element;
use super::rings::RingInfo;

/// A molecular graph backed by an undirected [`petgraph::graph::UnGraph`].
///
/// Molecules are constructed via [`parse_smiles`](super::smiles::parse_smiles)
/// or built manually with [`add_atom`](Self::add_atom) and
/// [`add_bond`](Self::add_bond).
pub struct Molecule {
    /// The underlying atom/bond graph.
    pub graph: UnGraph<Atom, Bond>,
    /// Lazily computed ring information (SSSR).
    ring_info: OnceCell<RingInfo>,
}

impl std::fmt::Debug for Molecule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Molecule")
            .field("graph", &self.graph)
            .finish()
    }
}

impl Clone for Molecule {
    fn clone(&self) -> Self {
        Self {
            graph: self.graph.clone(),
            ring_info: OnceCell::new(),
        }
    }
}

impl Molecule {
    /// Creates a new empty molecule.
    pub fn new() -> Self {
        Self {
            graph: UnGraph::new_undirected(),
            ring_info: OnceCell::new(),
        }
    }

    /// Get ring information (lazily computed).
    pub fn ring_info(&self) -> &RingInfo {
        self.ring_info.get_or_init(|| RingInfo::detect(&self.graph))
    }

    /// Check if an atom is in a ring.
    pub fn is_in_ring(&self, atom: NodeIndex) -> bool {
        self.ring_info().is_in_ring(atom)
    }

    /// Number of SSSR rings.
    pub fn num_rings(&self) -> usize {
        self.ring_info().num_rings()
    }

    /// Add an atom and return its index.
    pub fn add_atom(&mut self, atom: Atom) -> NodeIndex {
        self.graph.add_node(atom)
    }

    /// Add a bond between two atoms.
    pub fn add_bond(&mut self, a: NodeIndex, b: NodeIndex, bond: Bond) {
        self.graph.add_edge(a, b, bond);
    }

    /// Number of atoms (heavy + explicit H, not implicit H).
    pub fn atom_count(&self) -> usize {
        self.graph.node_count()
    }

    /// Number of bonds.
    pub fn bond_count(&self) -> usize {
        self.graph.edge_count()
    }

    /// Number of heavy (non-hydrogen) atoms.
    pub fn heavy_atom_count(&self) -> usize {
        self.graph
            .node_weights()
            .filter(|a| a.element.is_heavy())
            .count()
    }

    /// Total molecular weight (monoisotopic exact mass) including implicit hydrogens.
    pub fn molecular_weight(&self) -> f64 {
        self.graph.node_weights().map(|a| a.exact_mass()).sum()
    }

    /// Iterator over all atoms with their indices.
    pub fn atoms(&self) -> impl Iterator<Item = (NodeIndex, &Atom)> {
        self.graph.node_indices().map(move |i| (i, &self.graph[i]))
    }

    /// Iterator over all bonds with endpoint indices.
    pub fn bonds(&self) -> impl Iterator<Item = (NodeIndex, NodeIndex, &Bond)> {
        self.graph.edge_indices().map(move |e| {
            let (a, b) = self.graph.edge_endpoints(e).unwrap();
            (a, b, &self.graph[e])
        })
    }

    /// Degree of an atom (number of explicit bonds).
    pub fn degree(&self, idx: NodeIndex) -> usize {
        self.graph.neighbors(idx).count()
    }

    /// Total degree of an atom including implicit hydrogens.
    pub fn total_degree(&self, idx: NodeIndex) -> usize {
        self.degree(idx) + self.graph[idx].implicit_h as usize
    }

    /// Shortest-path distance matrix using Floyd-Warshall.
    /// Returns a map from (NodeIndex, NodeIndex) -> distance.
    pub fn distance_matrix(&self) -> std::collections::HashMap<(NodeIndex, NodeIndex), i64> {
        // Build a weighted copy where each edge has weight 1
        let weighted: UnGraph<(), i64> = self.graph.map(|_, _| (), |_, _| 1i64);
        floyd_warshall(&weighted, |e| *e.weight()).expect("no negative cycles in molecular graph")
    }

    /// Total number of atoms including implicit hydrogens.
    pub fn total_atom_count(&self) -> usize {
        self.atom_count()
            + self
                .graph
                .node_weights()
                .map(|a| a.implicit_h as usize)
                .sum::<usize>()
    }

    /// Get atom by index.
    pub fn atom(&self, idx: NodeIndex) -> &Atom {
        &self.graph[idx]
    }

    /// Count atoms of a given element.
    pub fn count_element(&self, element: Element) -> usize {
        self.graph
            .node_weights()
            .filter(|a| a.element == element)
            .count()
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}
