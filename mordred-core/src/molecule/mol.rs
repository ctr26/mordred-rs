use std::cell::OnceCell;
use std::collections::HashMap;

use petgraph::algo::floyd_warshall;
use petgraph::graph::{NodeIndex, UnGraph};

use super::atom::Atom;
use super::bond::{Bond, BondOrder};
use super::element::Element;
use super::rings::RingInfo;

/// Pre-computed molecular properties from a single pass over atoms and bonds.
pub struct MolecularProperties {
    pub element_counts: [u32; 26],
    pub heavy_atom_count: u32,
    pub implicit_h_sum: u32,
    pub total_atom_count: u32,
    pub molecular_weight: f64,
    pub single_bond_count: u32,
    pub double_bond_count: u32,
    pub triple_bond_count: u32,
    pub aromatic_bond_count: u32,
    pub total_bond_count: u32,
    pub degrees: Vec<u32>,
    pub halogen_count: u32,
    pub heteroatom_count: u32,
    pub hydrogen_count: u32,
}

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
    /// Lazily computed distance matrix (Floyd-Warshall).
    dist_cache: OnceCell<HashMap<(NodeIndex, NodeIndex), i64>>,
    /// Lazily computed molecular properties (fused counting).
    props_cache: OnceCell<MolecularProperties>,
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
            dist_cache: OnceCell::new(),
            props_cache: OnceCell::new(),
        }
    }
}

impl Molecule {
    /// Creates a new empty molecule.
    pub fn new() -> Self {
        Self {
            graph: UnGraph::new_undirected(),
            ring_info: OnceCell::new(),
            dist_cache: OnceCell::new(),
            props_cache: OnceCell::new(),
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

    /// Shortest-path distance matrix using Floyd-Warshall (lazily cached).
    pub fn distance_matrix(&self) -> &HashMap<(NodeIndex, NodeIndex), i64> {
        self.dist_cache.get_or_init(|| {
            let weighted: UnGraph<(), i64> = self.graph.map(|_, _| (), |_, _| 1i64);
            floyd_warshall(&weighted, |e| *e.weight())
                .expect("no negative cycles in molecular graph")
        })
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
        self.properties().element_counts[element.discriminant_index()] as usize
    }

    /// Cached molecular properties computed via SOA extraction + vectorized reductions.
    ///
    /// Extracts element discriminants, implicit H values, and masses into
    /// contiguous arrays, then processes them with tight loops that LLVM
    /// can auto-vectorize.
    pub fn properties(&self) -> &MolecularProperties {
        self.props_cache.get_or_init(|| {
            let node_count = self.graph.node_count();

            // SOA extraction: pull fields into contiguous arrays for auto-vectorization
            let mut discriminants = Vec::with_capacity(node_count);
            let mut implicit_h_values = Vec::with_capacity(node_count);
            let mut masses = Vec::with_capacity(node_count);
            let mut degrees = vec![0u32; node_count];

            for idx in self.graph.node_indices() {
                let atom = &self.graph[idx];
                discriminants.push(atom.element.discriminant_index() as u8);
                implicit_h_values.push(atom.implicit_h);
                masses.push(atom.exact_mass());
                degrees[idx.index()] = self.graph.neighbors(idx).count() as u32;
            }

            // Vectorized reductions over contiguous arrays
            let element_counts = super::simd::element_histogram(&discriminants);
            let implicit_h_sum = super::simd::sum_implicit_h(&implicit_h_values);
            let molecular_weight = super::simd::sum_molecular_weight(&masses);
            let heavy_atom_count = super::simd::count_heavy(&discriminants);
            let halogen_count = super::simd::count_halogens(&discriminants);
            let heteroatom_count = super::simd::count_heteroatoms(&discriminants);

            let hydrogen_count = element_counts[Element::H.discriminant_index()]
                + implicit_h_sum;
            let total_atom_count = node_count as u32 + implicit_h_sum;

            // Bond pass
            let mut single_bond_count = 0u32;
            let mut double_bond_count = 0u32;
            let mut triple_bond_count = 0u32;
            let mut aromatic_bond_count = 0u32;
            for edge in self.graph.edge_indices() {
                match self.graph[edge].order {
                    BondOrder::Single => single_bond_count += 1,
                    BondOrder::Double => double_bond_count += 1,
                    BondOrder::Triple => triple_bond_count += 1,
                    BondOrder::Aromatic => aromatic_bond_count += 1,
                }
            }
            let total_bond_count = self.graph.edge_count() as u32 + implicit_h_sum;

            MolecularProperties {
                element_counts,
                heavy_atom_count,
                implicit_h_sum,
                total_atom_count,
                molecular_weight,
                single_bond_count,
                double_bond_count,
                triple_bond_count,
                aromatic_bond_count,
                total_bond_count,
                degrees,
                halogen_count,
                heteroatom_count,
                hydrogen_count,
            }
        })
    }
}

impl Default for Molecule {
    fn default() -> Self {
        Self::new()
    }
}
