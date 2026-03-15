/// Aromatic perception for Kekulized SMILES input.
///
/// When a molecule is parsed from Kekulized SMILES (e.g. `C1=CC=CC=C1`),
/// aromaticity is not explicitly encoded. This module detects aromatic rings
/// from the bond topology and marks atoms/bonds accordingly, so that
/// Kekulized and aromatic SMILES produce identical descriptor values.
use petgraph::graph::NodeIndex;

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

    // Detect rings (we compute fresh since the lazy cache may not be populated,
    // and we're about to mutate the molecule which would invalidate it anyway).
    let ring_info = RingInfo::detect(&mol.graph);
    let rings = ring_info.rings().to_vec();

    // Collect which atoms and bonds to mark aromatic.
    let mut aromatic_atoms = std::collections::HashSet::new();
    let mut aromatic_bond_pairs: Vec<(NodeIndex, NodeIndex)> = Vec::new();

    for ring in &rings {
        let ring_size = ring.len();
        // Only consider 5- and 6-membered rings for simple Hückel aromaticity.
        if ring_size != 5 && ring_size != 6 {
            continue;
        }

        if is_aromatic_ring(mol, ring) {
            for &atom in ring {
                aromatic_atoms.insert(atom);
            }
            for i in 0..ring_size {
                let a = ring[i];
                let b = ring[(i + 1) % ring_size];
                aromatic_bond_pairs.push((a, b));
            }
        }
    }

    if aromatic_atoms.is_empty() {
        return;
    }

    // Mark atoms as aromatic.
    for &idx in &aromatic_atoms {
        mol.graph[idx].is_aromatic = true;
    }

    // Mark bonds as aromatic.
    for (a, b) in &aromatic_bond_pairs {
        if let Some(edge) = mol.graph.find_edge(*a, *b) {
            let bond = &mut mol.graph[edge];
            bond.order = BondOrder::Aromatic;
            bond.is_aromatic = true;
        }
    }

    // Recalculate implicit hydrogens for atoms whose aromaticity changed.
    recalculate_implicit_h(mol, &aromatic_atoms);
}

/// Check whether a ring has the bond pattern of an aromatic system.
///
/// Requirements:
/// - All atoms must be sp2-capable elements (C, N, O, S, Se).
/// - The ring must have alternating single/double bonds, OR all aromatic bonds.
/// - For 6-membered rings: classic Hückel (6 pi electrons).
/// - For 5-membered rings: need a heteroatom lone pair donor (N, O, S) or
///   the ring must have 2 double bonds + 1 heteroatom contributing 2 e-.
fn is_aromatic_ring(mol: &Molecule, ring: &[NodeIndex]) -> bool {
    let ring_size = ring.len();

    // Check all atoms are sp2-capable.
    for &idx in ring {
        let elem = mol.graph[idx].element;
        if !is_sp2_capable(elem) {
            return false;
        }
    }

    // Collect bond orders around the ring.
    let mut bond_orders: Vec<BondOrder> = Vec::with_capacity(ring_size);
    for i in 0..ring_size {
        let a = ring[i];
        let b = ring[(i + 1) % ring_size];
        let edge = mol.graph.find_edge(a, b);
        match edge {
            Some(e) => bond_orders.push(mol.graph[e].order),
            None => return false, // ring nodes not actually bonded
        }
    }

    // Check for alternating single/double bond pattern.
    let has_alternating = is_alternating_single_double(&bond_orders);

    if ring_size == 6 && has_alternating {
        // Classic benzene-like: alternating single/double in 6-membered ring
        // with all sp2-capable atoms -> 6 pi electrons -> aromatic.
        return true;
    }

    if ring_size == 5 && is_alternating_5_ring(ring, &bond_orders, mol) {
        return true;
    }

    false
}

/// Check for strict alternating single/double bonds.
fn is_alternating_single_double(orders: &[BondOrder]) -> bool {
    if orders.is_empty() {
        return false;
    }

    let n = orders.len();
    // Must have both single and double bonds.
    let n_single = orders.iter().filter(|&&o| o == BondOrder::Single).count();
    let n_double = orders.iter().filter(|&&o| o == BondOrder::Double).count();

    if n_single + n_double != n {
        return false; // has triple or other non-standard bond
    }

    // Check alternation: each bond must differ from its neighbors.
    for i in 0..n {
        let next = (i + 1) % n;
        if orders[i] == orders[next] {
            return false;
        }
    }
    true
}

/// Check whether a 5-membered ring is aromatic.
///
/// A 5-membered aromatic ring (like pyrrole, furan, thiophene) has:
/// - 2 double bonds and 1 heteroatom (N, O, S) that contributes a lone pair
///   to make 6 pi electrons (2*2 + 2 = 6), OR
/// - Alternating single/double with a heteroatom.
fn is_alternating_5_ring(
    ring: &[NodeIndex],
    bond_orders: &[BondOrder],
    mol: &Molecule,
) -> bool {
    let n_double = bond_orders
        .iter()
        .filter(|&&o| o == BondOrder::Double)
        .count();

    if n_double != 2 {
        return false;
    }

    // Need exactly one heteroatom (N, O, S) that donates a lone pair.
    // The heteroatom should be the one flanked by two single bonds.
    let mut heteroatom_count = 0;
    for (i, &idx) in ring.iter().enumerate() {
        let elem = mol.graph[idx].element;
        if matches!(elem, Element::N | Element::O | Element::S | Element::Se) {
            // Check that this heteroatom is flanked by single bonds
            // (it contributes the lone pair, not a pi bond).
            let prev_bond = bond_orders[(i + ring.len() - 1) % ring.len()];
            let next_bond = bond_orders[i];
            if prev_bond == BondOrder::Single && next_bond == BondOrder::Single {
                heteroatom_count += 1;
            }
        }
    }

    heteroatom_count == 1
}

/// Check if an element is capable of sp2 hybridization.
fn is_sp2_capable(elem: Element) -> bool {
    matches!(elem, Element::C | Element::N | Element::O | Element::S | Element::Se)
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
        let aromatic_count = mol
            .graph
            .node_weights()
            .filter(|a| a.is_aromatic)
            .count();
        assert_eq!(aromatic_count, 5, "pyrrole should have 5 aromatic atoms");
    }
}
