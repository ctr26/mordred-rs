/// Integration tests for mordred-core with known molecules.
///
/// Reference values computed from molecular structure (not Python mordred)
/// to verify our SMILES parser + descriptors work end-to-end.
use mordred_core::descriptor::connectivity::{Chi0, Chi1};
use mordred_core::descriptor::constitutional::{
    AtomCount, BondCount, HeavyAtomCount, MolecularWeight,
};
use mordred_core::descriptor::topological::{WienerIndex, ZagrebIndex1, ZagrebIndex2};
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
    // 6 aromatic bonds in a ring
    assert_eq!(BondCount.calculate(&mol).unwrap(), 6.0);
}

#[test]
fn benzene_mw() {
    let mol = parse_smiles("c1ccccc1").unwrap();
    let mw = MolecularWeight.calculate(&mol).unwrap();
    // C6H6 = 6*12.011 + 6*1.008 = 78.114
    assert!((mw - 78.114).abs() < 0.1, "benzene MW={}", mw);
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
    // C2H6O = 2*12.011 + 6*1.008 + 15.999 = 46.069
    assert!((mw - 46.069).abs() < 0.1, "ethanol MW={}", mw);
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
    // Caffeine MW ≈ 194.19
    assert!((mw - 194.19).abs() < 1.0, "caffeine MW={}", mw);
}

// ─── DescriptorSet integration ───

#[test]
fn descriptor_set_calculates_all() {
    let set = DescriptorSet::all();
    let mol = parse_smiles("CCO").unwrap();
    let results = set.calculate(&mol);
    // Should have 31 descriptors (4 constitutional + 3 topological + 2 connectivity + 7 atom count + 4 bond count + 11 ring count)
    assert_eq!(results.len(), 31);
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
    assert!(names.contains(&"WienerIndex"));
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
    assert!((mw - 18.015).abs() < 0.01);
}
