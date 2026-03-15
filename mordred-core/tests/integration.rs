/// Integration tests for mordred-core with known molecules.
///
/// Reference values computed from molecular structure (not Python mordred)
/// to verify our SMILES parser + descriptors work end-to-end.
use mordred_core::descriptor::bond_count::AromaticBondCount;
use mordred_core::descriptor::constitutional::{
    AtomCount, BondCount, HeavyAtomCount, MolecularWeight,
};
use mordred_core::descriptor::topological::{WienerIndex, ZagrebIndex1};
use mordred_core::{Descriptor, DescriptorSet, parse_smiles};

// ─── Benzene: c1ccccc1 ───

#[test]
fn benzene_heavy_atoms() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 6.0);
}

#[test]
fn benzene_bonds() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    // 6 aromatic bonds + 6 implicit H bonds = 12
    assert_eq!(BondCount.calculate(&mol).unwrap(), 12.0);
}

#[test]
fn benzene_mw() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    let mw = MolecularWeight.calculate(&mol).unwrap();
    // C6H6 = 6*12.0 + 6*1.00783 ≈ 78.047
    assert!((mw - 78.047).abs() < 0.01, "benzene MW={}", mw);
}

#[test]
fn benzene_wiener() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    let w = WienerIndex.calculate(&mol).unwrap();
    // Wiener index for 6-cycle:
    // distances: 6 pairs at dist 1, 6 pairs at dist 2, 3 pairs at dist 3
    // W = 6*1 + 6*2 + 3*3 = 6 + 12 + 9 = 27
    assert_eq!(w, 27.0);
}

#[test]
fn benzene_zagreb1() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    let z = ZagrebIndex1.calculate(&mol).unwrap();
    // Each atom has degree 2: 6 * 2^2 = 24
    assert_eq!(z, 24.0);
}

// ─── Ethanol: CCO ───

#[test]
fn ethanol_atoms() {
    let mol = parse_smiles("CCO").unwrap();
    assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 3.0);
    // C(3H) + C(2H) + O(1H) = 3 + 3+2+1 = 9 total atoms
    assert_eq!(AtomCount.calculate(&mol).unwrap(), 9.0);
}

#[test]
fn ethanol_mw() {
    let mol = parse_smiles("CCO").unwrap();
    let mw = MolecularWeight.calculate(&mol).unwrap();
    // C2H6O = 2*12.0 + 6*1.00783 + 15.995 ≈ 46.042
    assert!((mw - 46.042).abs() < 0.01, "ethanol MW={}", mw);
}

#[test]
fn ethanol_wiener() {
    let mol = parse_smiles("CCO").unwrap();
    let w = WienerIndex.calculate(&mol).unwrap();
    // 3 atoms in a chain: W = 1 + 2 + 1 = 4
    assert_eq!(w, 4.0);
}

// ─── Caffeine: Cn1cnc2c1c(=O)n(c(=O)n2C)C ───

#[test]
fn caffeine_heavy_atoms() {
    let mol = parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
    // C8H10N4O2 — 14 heavy atoms
    assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 14.0);
}

#[test]
fn caffeine_mw() {
    let mol = parse_smiles("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
    let mw = MolecularWeight.calculate(&mol).unwrap();
    // C8H10N4O2 = 8*12.0 + 10*1.00783 + 4*14.003 + 2*15.995 ≈ 194.080
    assert!((mw - 194.080).abs() < 0.1, "caffeine MW={}", mw);
}

// ─── DescriptorSet integration ───

#[test]
fn descriptor_set_calculates_all() {
    let set = DescriptorSet::all();
    let mol = parse_smiles("CCO").unwrap();
    let results = set.calculate(&mol);
    // Should have 37 descriptors (4 constitutional + 3 topological + 2 connectivity + 13 atom count + 4 bond count + 11 ring count)
    assert_eq!(results.len(), 37);
    // All should succeed
    for (name, result) in &results {
        assert!(result.is_ok(), "descriptor {} failed", name);
    }
}

#[test]
fn descriptor_set_names() {
    let set = DescriptorSet::all();
    let names = set.names();
    assert!(names.contains(&"MW"));
    assert!(names.contains(&"WPath"));
    assert!(names.contains(&"Chi0"));
}

// ─── Acetic acid: CC(=O)O ───

#[test]
fn acetic_acid_structure() {
    let mol = parse_smiles("CC(=O)O").unwrap();
    assert_eq!(mol.atom_count(), 4); // C, C, O, O
    assert_eq!(mol.bond_count(), 3); // C-C, C=O, C-O
    assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 4.0);
}

// ─── Water: O ───

#[test]
fn water() {
    let mol = parse_smiles("O").unwrap();
    assert_eq!(HeavyAtomCount.calculate(&mol).unwrap(), 1.0);
    assert_eq!(AtomCount.calculate(&mol).unwrap(), 3.0); // O + 2H
    let mw = MolecularWeight.calculate(&mol).unwrap();
    assert!((mw - 18.011).abs() < 0.01);
}

// ─── Kekulized benzene: C1=CC=CC=C1 (aromatic perception) ───

#[test]
fn kekulized_benzene_aromatic_bonds() {
    let mol = parse_smiles("C1=CC=CC=C1").unwrap();
    assert_eq!(
        AromaticBondCount.calculate(&mol).unwrap(),
        6.0,
        "Kekulized benzene should have 6 aromatic bonds"
    );
}

#[test]
fn kekulized_benzene_no_double_bonds() {
    use mordred_core::descriptor::bond_count::DoubleBondCount;
    let mol = parse_smiles("C1=CC=CC=C1").unwrap();
    assert_eq!(
        DoubleBondCount.calculate(&mol).unwrap(),
        0.0,
        "Kekulized benzene should have 0 double bonds after aromaticity perception"
    );
}

#[test]
fn kekulized_benzene_mw_matches_aromatic() {
    let kek = parse_smiles("C1=CC=CC=C1").unwrap();
    let aro = parse_smiles("c1ccccc1").unwrap();
    let mw_kek = MolecularWeight.calculate(&kek).unwrap();
    let mw_aro = MolecularWeight.calculate(&aro).unwrap();
    assert!(
        (mw_kek - mw_aro).abs() < 0.001,
        "MW should match: Kekulized={} aromatic={}",
        mw_kek,
        mw_aro
    );
}

#[test]
fn kekulized_benzene_atom_count_matches_aromatic() {
    let kek = parse_smiles("C1=CC=CC=C1").unwrap();
    let aro = parse_smiles("c1ccccc1").unwrap();
    assert_eq!(
        AtomCount.calculate(&kek).unwrap(),
        AtomCount.calculate(&aro).unwrap(),
        "Total atom count should match between Kekulized and aromatic benzene"
    );
}

// ─── Ibuprofen: CC(C)CC1=CC=C(C=C1)C(C)C(=O)O ───

#[test]
fn ibuprofen_phenyl_aromatic_bonds() {
    let mol = parse_smiles("CC(C)CC1=CC=C(C=C1)C(C)C(=O)O").unwrap();
    assert_eq!(
        AromaticBondCount.calculate(&mol).unwrap(),
        6.0,
        "Ibuprofen phenyl ring should have 6 aromatic bonds"
    );
}

// ─── Cyclohexane: C1CCCCC1 (must NOT be aromatic) ───

#[test]
fn cyclohexane_no_aromatic_bonds() {
    let mol = parse_smiles("C1CCCCC1").unwrap();
    assert_eq!(
        AromaticBondCount.calculate(&mol).unwrap(),
        0.0,
        "Cyclohexane should have 0 aromatic bonds"
    );
}
