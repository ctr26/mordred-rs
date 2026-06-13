/// Aromatic perception for Kekulized SMILES input.
///
/// When a molecule is parsed from Kekulized SMILES (e.g. `C1=CC=CC=C1`),
/// aromaticity is not explicitly encoded. This module detects aromatic rings
/// from the bond topology and marks atoms/bonds accordingly, so that
/// Kekulized and aromatic SMILES produce identical descriptor values.
use petgraph::graph::NodeIndex;
use petgraph::visit::EdgeRef;

use super::bond::BondOrder;
use super::element::Element;
use super::mol::Molecule;
use super::rings::RingInfo;

/// Perceive aromaticity from ring structure and bond alternation.
///
/// For each SSSR ring of size 5 or 6, check whether the ring has alternating
/// single/double bonds and all ring atoms are sp2-capable (C, N, O, S).
/// If so, mark the ring atoms as aromatic and set ring bond orders to
/// [`BondOrder::Aromatic`].
///
/// This function should be called after the SMILES parser has built the
/// molecule and filled implicit hydrogens. It will recalculate implicit
/// hydrogens for any atoms whose aromaticity status changes.
pub fn perceive_aromaticity(mol: &mut Molecule) {
    // If the molecule already has aromatic atoms (from lowercase SMILES),
    // skip perception — aromaticity was already explicit.
    let has_aromatic = mol.graph.node_weights().any(|a| a.is_aromatic);
    if has_aromatic {
        return;
    }

    let ring_info = RingInfo::detect(&mol.graph);
    let rings: Vec<Vec<NodeIndex>> = ring_info
        .rings()
        .iter()
        .filter(|r| r.len() == 5 || r.len() == 6)
        .cloned()
        .collect();

    // Iterate to fixpoint: a ring may become aromatic only after a neighboring
    // fused ring has been marked (its bonds become Aromatic, which then count
    // as in-ring π contributors). EllagicAcid's lactones and the chromenylium
    // core of Cyanidin need this propagation.
    let mut aromatic_atoms: std::collections::HashSet<NodeIndex> = std::collections::HashSet::new();
    let mut newly_marked = true;
    while newly_marked {
        newly_marked = false;
        for ring in &rings {
            // Skip rings whose every atom is already aromatic — they have no
            // new state to discover and their bonds are already aromatic.
            if ring.iter().all(|a| mol.graph[*a].is_aromatic) {
                continue;
            }
            if !is_aromatic_ring(mol, ring) {
                continue;
            }
            for &atom in ring {
                if aromatic_atoms.insert(atom) {
                    mol.graph[atom].is_aromatic = true;
                }
            }
            for i in 0..ring.len() {
                if let Some(edge) = mol.graph.find_edge(ring[i], ring[(i + 1) % ring.len()]) {
                    let bond = &mut mol.graph[edge];
                    bond.order = BondOrder::Aromatic;
                    bond.is_aromatic = true;
                }
            }
            newly_marked = true;
        }
    }

    if !aromatic_atoms.is_empty() {
        recalculate_implicit_h(mol, &aromatic_atoms);
    }
}

/// Check whether a ring is aromatic via Hückel's rule on per-atom π contributions.
///
/// Each ring atom contributes π electrons based on its local context:
/// - In-ring double (or aromatic) bond: 1 electron (pyridine-like sp²).
/// - Exocyclic double bond (=O, =N, =C): 0 to the ring, but sp² maintained
///   (handles carbonyl carbons in uracil/caffeine-style pyrimidinediones).
/// - In-ring singles only + N/O/S lone-pair donor: 2 electrons (pyrrole-like).
/// - Otherwise: not sp² → ring is not aromatic.
///
/// The ring is aromatic iff the π-electron sum satisfies Hückel's rule
/// (`π = 4n + 2`, n ≥ 0 → π ∈ {2, 6, 10, …}) and all atoms are sp²-capable.
fn is_aromatic_ring(mol: &Molecule, ring: &[NodeIndex]) -> bool {
    let n = ring.len();

    for &idx in ring {
        if !is_sp2_capable(mol.graph[idx].element) {
            return false;
        }
    }

    let mut bond_orders: Vec<BondOrder> = Vec::with_capacity(n);
    for i in 0..n {
        match mol.graph.find_edge(ring[i], ring[(i + 1) % n]) {
            Some(e) => bond_orders.push(mol.graph[e].order),
            None => return false,
        }
    }

    let ring_set: std::collections::HashSet<NodeIndex> = ring.iter().copied().collect();
    let mut pi: u32 = 0;
    for (i, &atom) in ring.iter().enumerate().take(n) {
        match ring_pi_contribution(mol, atom, &ring_set, &bond_orders, i) {
            Some(p) => pi += p,
            None => return false,
        }
    }

    pi >= 2 && (pi - 2) % 4 == 0
}

/// π-electron contribution for one ring atom, or `None` if it can't be sp².
fn ring_pi_contribution(
    mol: &Molecule,
    idx: NodeIndex,
    ring_set: &std::collections::HashSet<NodeIndex>,
    bond_orders: &[BondOrder],
    i: usize,
) -> Option<u32> {
    let n = bond_orders.len();
    let prev = bond_orders[(i + n - 1) % n];
    let next = bond_orders[i];

    let in_ring_pi = matches!(
        prev,
        BondOrder::Double | BondOrder::Aromatic | BondOrder::Triple
    ) || matches!(
        next,
        BondOrder::Double | BondOrder::Aromatic | BondOrder::Triple
    );
    if in_ring_pi {
        return Some(1);
    }

    let has_exocyclic_double = mol.graph.edges(idx).any(|e| {
        let other = if e.source() == idx {
            e.target()
        } else {
            e.source()
        };
        !ring_set.contains(&other) && e.weight().order == BondOrder::Double
    });
    if has_exocyclic_double {
        return Some(0);
    }

    match mol.graph[idx].element {
        Element::N | Element::O | Element::S | Element::Se => Some(2),
        _ => None,
    }
}

/// Check if an element is capable of sp2 hybridization.
fn is_sp2_capable(elem: Element) -> bool {
    matches!(
        elem,
        Element::C | Element::N | Element::O | Element::S | Element::Se
    )
}

/// Recalculate implicit hydrogens for atoms that became aromatic.
///
/// When a bond changes from Single/Double to Aromatic, the valence
/// contribution changes. Aromatic bonds contribute 1 to valence
/// (like single bonds), but aromatic atoms also have an implicit pi
/// contribution that uses up 1 valence unit.
fn recalculate_implicit_h(
    mol: &mut Molecule,
    aromatic_atoms: &std::collections::HashSet<NodeIndex>,
) {
    for &idx in aromatic_atoms {
        let valence = mol.graph[idx].element.default_valence();
        if valence == 0 {
            continue;
        }

        let bond_valence: u8 = mol
            .graph
            .edges(idx)
            .map(|e| e.weight().order.valence_contribution())
            .sum();

        let atom = &mut mol.graph[idx];
        // Aromatic atoms have an extra pi electron contribution.
        // e.g. aromatic C: valence 4, 2 aromatic bonds (contribute 1 each) + 1 pi = 3
        // so implicit_h = 4 - 3 = 1 for unsubstituted aromatic C.
        atom.implicit_h = valence.saturating_sub(bond_valence + 1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_kekulized_benzene_is_aromatic() {
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        // All atoms should be aromatic.
        for node in mol.graph.node_indices() {
            assert!(
                mol.graph[node].is_aromatic,
                "atom {} should be aromatic",
                node.index()
            );
        }
        // All bonds should be aromatic.
        for edge in mol.graph.edge_indices() {
            assert_eq!(
                mol.graph[edge].order,
                BondOrder::Aromatic,
                "bond should be aromatic"
            );
        }
    }

    #[test]
    fn test_kekulized_benzene_implicit_h() {
        let mol = parse_smiles("C1=CC=CC=C1").unwrap();
        // Each aromatic C should have 1 implicit H.
        for node in mol.graph.node_indices() {
            assert_eq!(
                mol.graph[node].implicit_h, 1,
                "aromatic C should have 1 implicit H"
            );
        }
    }

    #[test]
    fn test_aromatic_smiles_still_works() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        for node in mol.graph.node_indices() {
            assert!(mol.graph[node].is_aromatic);
        }
    }

    #[test]
    fn test_cyclohexane_not_aromatic() {
        let mol = parse_smiles("C1CCCCC1").unwrap();
        for node in mol.graph.node_indices() {
            assert!(
                !mol.graph[node].is_aromatic,
                "cyclohexane atoms should not be aromatic"
            );
        }
    }

    #[test]
    fn test_kekulized_benzene_same_as_aromatic() {
        let kek = parse_smiles("C1=CC=CC=C1").unwrap();
        let aro = parse_smiles("c1ccccc1").unwrap();

        // Same number of atoms and bonds.
        assert_eq!(kek.atom_count(), aro.atom_count());
        assert_eq!(kek.bond_count(), aro.bond_count());

        // Same aromatic bond count.
        let kek_arom = kek
            .bonds()
            .filter(|(_, _, b)| b.order == BondOrder::Aromatic)
            .count();
        let aro_arom = aro
            .bonds()
            .filter(|(_, _, b)| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(kek_arom, aro_arom);

        // Same implicit H.
        let kek_h: u8 = kek.graph.node_weights().map(|a| a.implicit_h).sum();
        let aro_h: u8 = aro.graph.node_weights().map(|a| a.implicit_h).sum();
        assert_eq!(kek_h, aro_h);
    }

    #[test]
    fn test_pyrrole_aromatic() {
        // Kekulized pyrrole: C1=CC=CN1 (5-membered ring with NH)
        // Note: this is a simplified test; pyrrole in Kekulized form
        // is N1C=CC=C1 or similar.
        let mol = parse_smiles("C1=CC=C[NH]1").unwrap();
        // If the ring detection works correctly, atoms should be aromatic.
        let aromatic_count = mol.graph.node_weights().filter(|a| a.is_aromatic).count();
        assert_eq!(aromatic_count, 5, "pyrrole should have 5 aromatic atoms");
    }

    #[test]
    fn test_caffeine_pyrimidinedione_aromatic() {
        // Caffeine's 6-ring is a pyrimidinedione: 2 ring N's flanked by single
        // bonds (pyrrole-like, 2 π each), 2 exocyclic carbonyls (sp² via =O,
        // 0 π each), 2 fused C's with in-ring double (1 π each). Hückel π = 6.
        let mol = parse_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        let aromatic_atoms = mol.graph.node_weights().filter(|a| a.is_aromatic).count();
        // Caffeine has 9 aromatic atoms across the imidazole (5) and
        // pyrimidinedione (6) fused-system, sharing 2 atoms → 5 + 6 - 2 = 9.
        assert_eq!(
            aromatic_atoms, 9,
            "caffeine should have 9 aromatic atoms (imidazole + pyrimidinedione fused), got {aromatic_atoms}"
        );
    }

    #[test]
    fn test_ellagic_acid_fully_aromatic() {
        // Ellagic acid: 4 fused 6-rings (2 benzenes + 2 fused 6-membered
        // lactones). Aromaticity propagates through the fused system via the
        // fixpoint loop — lactone rings only register aromatic after their
        // adjacent benzene rings are marked.
        let mol = parse_smiles("C1=C2C3=C(C(=C1O)O)OC(=O)C4=CC(=C(C(=C43)OC2=O)O)O").unwrap();
        let aromatic_atoms = mol.graph.node_weights().filter(|a| a.is_aromatic).count();
        assert!(
            aromatic_atoms >= 14,
            "ellagic acid should have ≥14 aromatic atoms via fused-ring propagation, got {aromatic_atoms}"
        );
    }
}
