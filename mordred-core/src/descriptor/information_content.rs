/// Information content descriptors (IC0, TIC0, SIC0, CIC0, BIC0).
///
/// Shannon entropy-based descriptors that quantify the diversity of
/// atom/bond equivalence classes in a molecule.
///
/// Reference: mordred Python `InformationContent.py`
/// Theory: Basak, S.C. et al. (1988) — graph-theoretic information content.
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Shannon entropy H = -Σ p_i * log2(p_i) from a frequency distribution.
fn shannon_entropy(counts: &[u32], total: u32) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let n = total as f64;
    counts
        .iter()
        .filter(|&&c| c > 0)
        .map(|&c| {
            let p = c as f64 / n;
            -p * p.log2()
        })
        .sum()
}

/// Compute atom equivalence class counts based on element type (order 0).
/// Includes implicit hydrogens in both the H class count and the total.
fn atom_equivalence_counts(mol: &Molecule) -> (Vec<u32>, u32) {
    use crate::molecule::Element;
    let props = mol.properties();
    let mut counts: Vec<u32> = Vec::new();
    for (i, &c) in props.element_counts.iter().enumerate() {
        if i == Element::H.discriminant_index() {
            // Use full hydrogen count (explicit + implicit)
            if props.hydrogen_count > 0 {
                counts.push(props.hydrogen_count);
            }
        } else if c > 0 {
            counts.push(c);
        }
    }
    let total = props.heavy_atom_count + props.hydrogen_count;
    (counts, total)
}

/// Sum of bond orders (aromatic = 1.5; each implicit H contributes one single).
///
/// Mirrors upstream's `sum(b.GetBondTypeAsDouble() for b in mol.GetBonds())`
/// under `explicit_hydrogens=True`. For pure 6-membered aromatic rings the
/// 1.5-weight sum equals the kekulized single+double sum; they may differ for
/// odd-membered aromatic rings (resolved by upstream via `kekulize=True` on
/// the RDKit side).
fn sum_bond_orders(mol: &Molecule) -> f64 {
    let explicit: f64 = mol.bonds().map(|(_, _, b)| b.order.as_f64()).sum();
    explicit + f64::from(mol.properties().implicit_h_sum)
}

/// IC0: Mean information content (Shannon entropy of atom element classes).
pub struct InformationContent;

impl Descriptor for InformationContent {
    fn name(&self) -> &str {
        "IC0"
    }
    fn description(&self) -> &str {
        "Mean information content (order 0)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let (counts, total) = atom_equivalence_counts(mol);
        Ok(shannon_entropy(&counts, total))
    }
}

/// TIC0: Total information content = N * IC0.
pub struct TotalInformationContent;

impl Descriptor for TotalInformationContent {
    fn name(&self) -> &str {
        "TIC0"
    }
    fn description(&self) -> &str {
        "Total information content (order 0)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let (counts, total) = atom_equivalence_counts(mol);
        let ic = shannon_entropy(&counts, total);
        Ok(total as f64 * ic)
    }
}

/// SIC0: Structural information content = IC0 / log2(N).
pub struct StructuralInformationContent;

impl Descriptor for StructuralInformationContent {
    fn name(&self) -> &str {
        "SIC0"
    }
    fn description(&self) -> &str {
        "Structural information content (order 0)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let (counts, total) = atom_equivalence_counts(mol);
        if total <= 1 {
            return Ok(0.0);
        }
        let ic = shannon_entropy(&counts, total);
        Ok(ic / (total as f64).log2())
    }
}

/// CIC0: Complementary information content = log2(N) - IC0.
pub struct ComplementaryInformationContent;

impl Descriptor for ComplementaryInformationContent {
    fn name(&self) -> &str {
        "CIC0"
    }
    fn description(&self) -> &str {
        "Complementary information content (order 0)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let (counts, total) = atom_equivalence_counts(mol);
        if total <= 1 {
            return Ok(0.0);
        }
        let ic = shannon_entropy(&counts, total);
        Ok((total as f64).log2() - ic)
    }
}

/// BIC0: Bonding information content = IC0 / log2(Σ bond_orders).
///
/// Upstream definition (mordred Python `BondingIC`): the atom-level Shannon
/// entropy IC0 normalized by the log of the sum of weighted bond orders, not a
/// separate bond-class entropy. With `explicit_hydrogens=True` upstream
/// includes implicit-H single bonds in the sum.
pub struct BondInformationContent;

impl Descriptor for BondInformationContent {
    fn name(&self) -> &str {
        "BIC0"
    }
    fn description(&self) -> &str {
        "Bond information content (order 0)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let (counts, total) = atom_equivalence_counts(mol);
        if total <= 1 {
            return Ok(0.0);
        }
        let b = sum_bond_orders(mol);
        if b <= 1.0 {
            return Ok(0.0);
        }
        let ic = shannon_entropy(&counts, total);
        Ok(ic / b.log2())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_ic0_methane() {
        // CH4: 1 C + 4 H = 5 atoms, 2 classes
        // p(C)=1/5, p(H)=4/5
        // IC = -(1/5)*log2(1/5) - (4/5)*log2(4/5)
        let mol = parse_smiles("C").unwrap();
        let ic = InformationContent.calculate(&mol).unwrap();
        let expected =
            -(1.0_f64 / 5.0) * (1.0_f64 / 5.0).log2() - (4.0_f64 / 5.0) * (4.0_f64 / 5.0).log2();
        assert!(
            (ic - expected).abs() < 1e-10,
            "IC0 methane: got {ic}, expected {expected}"
        );
    }

    #[test]
    fn test_ic0_pure_element() {
        // All same element: entropy should be 0 (only H atoms)
        // H2: [H][H] — 2 H atoms, 1 class → IC = 0
        let mol = parse_smiles("[H][H]").unwrap();
        let ic = InformationContent.calculate(&mol).unwrap();
        assert!((ic - 0.0).abs() < 1e-10, "IC0 H2 should be 0, got {ic}");
    }

    #[test]
    fn test_tic0_methane() {
        let mol = parse_smiles("C").unwrap();
        let tic = TotalInformationContent.calculate(&mol).unwrap();
        let expected_ic =
            -(1.0_f64 / 5.0) * (1.0_f64 / 5.0).log2() - (4.0_f64 / 5.0) * (4.0_f64 / 5.0).log2();
        let expected = 5.0 * expected_ic;
        assert!(
            (tic - expected).abs() < 1e-10,
            "TIC0 methane: got {tic}, expected {expected}"
        );
    }

    #[test]
    fn test_sic0_methane() {
        let mol = parse_smiles("C").unwrap();
        let sic = StructuralInformationContent.calculate(&mol).unwrap();
        let expected_ic =
            -(1.0_f64 / 5.0) * (1.0_f64 / 5.0).log2() - (4.0_f64 / 5.0) * (4.0_f64 / 5.0).log2();
        let expected = expected_ic / 5.0_f64.log2();
        assert!(
            (sic - expected).abs() < 1e-10,
            "SIC0 methane: got {sic}, expected {expected}"
        );
    }

    #[test]
    fn test_cic0_methane() {
        let mol = parse_smiles("C").unwrap();
        let cic = ComplementaryInformationContent.calculate(&mol).unwrap();
        let expected_ic =
            -(1.0_f64 / 5.0) * (1.0_f64 / 5.0).log2() - (4.0_f64 / 5.0) * (4.0_f64 / 5.0).log2();
        let expected = 5.0_f64.log2() - expected_ic;
        assert!(
            (cic - expected).abs() < 1e-10,
            "CIC0 methane: got {cic}, expected {expected}"
        );
    }

    #[test]
    fn test_sic_plus_cic_equals_log2n() {
        // SIC0 * log2(N) + CIC0 = log2(N)
        // i.e., IC0 + CIC0 = log2(N)
        let mol = parse_smiles("CCO").unwrap();
        let ic = InformationContent.calculate(&mol).unwrap();
        let cic = ComplementaryInformationContent.calculate(&mol).unwrap();
        let (_, total) = atom_equivalence_counts(&mol);
        let log2n = (total as f64).log2();
        assert!(
            (ic + cic - log2n).abs() < 1e-10,
            "IC0 + CIC0 should equal log2(N): {ic} + {cic} vs {log2n}"
        );
    }

    #[test]
    fn test_bic0_ethanol() {
        // Ethanol CCO: 9 atoms (2C + 1O + 6H), IC0 = entropy([2,1,6]/9) ≈ 1.22440.
        // Sum bond orders = 2 explicit singles + 6 implicit H = 8; log2(8) = 3.
        // BIC0 ≈ 1.22440 / 3 ≈ 0.40813.
        let mol = parse_smiles("CCO").unwrap();
        let bic = BondInformationContent.calculate(&mol).unwrap();
        let expected = 0.408_131_67;
        assert!(
            (bic - expected).abs() < 1e-6,
            "BIC0 ethanol: got {bic}, expected {expected}"
        );
    }

    #[test]
    fn test_bic0_formaldehyde() {
        // Formaldehyde C=O: 4 atoms (C+O+2H), IC0 = entropy([1,1,2]/4) = 1.5.
        // Sum bond orders = 1 double (2.0) + 2 implicit H singles = 4; log2(4) = 2.
        // BIC0 = 1.5 / 2 = 0.75.
        let mol = parse_smiles("C=O").unwrap();
        let bic = BondInformationContent.calculate(&mol).unwrap();
        let expected = 0.75;
        assert!(
            (bic - expected).abs() < 1e-6,
            "BIC0 formaldehyde: got {bic}, expected {expected}"
        );
    }

    #[test]
    fn test_ic0_ethanol() {
        // Ethanol CCO: 2C + 1O + 6H = 9 atoms
        // p(C) = 2/9, p(O) = 1/9, p(H) = 6/9
        let mol = parse_smiles("CCO").unwrap();
        let ic = InformationContent.calculate(&mol).unwrap();
        let expected = -(2.0 / 9.0) * (2.0_f64 / 9.0).log2()
            - (1.0 / 9.0) * (1.0_f64 / 9.0).log2()
            - (6.0 / 9.0) * (6.0_f64 / 9.0).log2();
        assert!(
            (ic - expected).abs() < 1e-10,
            "IC0 ethanol: got {ic}, expected {expected}"
        );
    }

    #[test]
    fn test_sic0_single_atom() {
        // Single atom: N ≤ 1 returns 0
        let mol = parse_smiles("[He]").unwrap();
        let sic = StructuralInformationContent.calculate(&mol).unwrap();
        assert!((sic - 0.0).abs() < 1e-10);
    }
}
