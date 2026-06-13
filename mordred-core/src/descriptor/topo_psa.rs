/// Topological polar surface area (TPSA) — Ertl, Rohde & Selzer (2000).
///
/// Sums per-atom PSA contributions from a fixed lookup keyed by
/// `(element, n_H, n_neighbors, n_single, n_double, n_triple, n_aromatic,
/// formal_charge, in_3-ring)`. Two flavours:
///
/// - `TopoPSA(NO)` — N and O atoms only (this is the classic Ertl TPSA).
/// - `TopoPSA` — N, O, plus S and P contributions.
///
/// Reference: Ertl, P.; Rohde, B.; Selzer, P. J. Med. Chem. 2000, 43, 3714–3717.
/// (DOI: 10.1021/jm000942e). The numeric table is the RDKit `CalcTPSA`
/// implementation in `Code/GraphMol/Descriptors/MolSurf.cpp`, which is the
/// table mordred-py delegates to.
use petgraph::visit::EdgeRef;

use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{BondOrder, Element, Molecule};

/// Per-atom bond-environment counts used to key the Ertl PSA table.
struct AtomEnv {
    n_h: i32,
    n_neighbors: i32,
    n_single: i32,
    n_double: i32,
    n_triple: i32,
    n_aromatic: i32,
}

fn atom_env(mol: &Molecule, idx: petgraph::graph::NodeIndex) -> AtomEnv {
    let atom = mol.atom(idx);
    let mut n_h = atom.implicit_h as i32;
    let mut n_neighbors = 0i32;
    let mut n_single = 0i32;
    let mut n_double = 0i32;
    let mut n_triple = 0i32;
    let mut n_aromatic = 0i32;

    for edge in mol.graph.edges(idx) {
        let other = if edge.source() == idx {
            edge.target()
        } else {
            edge.source()
        };
        let bond = edge.weight();
        // Explicit H atoms count toward nH, not toward degree (matches the
        // RDKit reference: H bonds decrement nNbrs and increment nHs).
        if mol.atom(other).element == Element::H {
            n_h += 1;
            continue;
        }
        n_neighbors += 1;
        if bond.is_aromatic {
            n_aromatic += 1;
        } else {
            match bond.order {
                BondOrder::Single => n_single += 1,
                BondOrder::Double => n_double += 1,
                BondOrder::Triple => n_triple += 1,
                BondOrder::Aromatic => n_aromatic += 1,
            }
        }
    }

    AtomEnv {
        n_h,
        n_neighbors,
        n_single,
        n_double,
        n_triple,
        n_aromatic,
    }
}

/// Ertl/RDKit PSA contribution for N and O atoms.
fn psa_contrib_no(mol: &Molecule, idx: petgraph::graph::NodeIndex) -> f64 {
    let atom = mol.atom(idx);
    let env = atom_env(mol, idx);
    let chg = atom.charge as i32;
    let in3 = mol
        .ring_info()
        .smallest_ring_size(idx)
        .is_some_and(|s| s == 3);

    let mut tmp: f64 = -1.0;
    match atom.element {
        Element::N => {
            match env.n_neighbors {
                1 => {
                    if env.n_h == 0 && chg == 0 && env.n_triple == 1 {
                        tmp = 23.79;
                    } else if env.n_h == 1 && chg == 0 && env.n_double == 1 {
                        tmp = 23.85;
                    } else if env.n_h == 2 && chg == 0 && env.n_single == 1 {
                        tmp = 26.02;
                    } else if env.n_h == 2 && chg == 1 && env.n_double == 1 {
                        tmp = 25.59;
                    } else if env.n_h == 3 && chg == 1 && env.n_single == 1 {
                        tmp = 27.64;
                    }
                }
                2 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 1 && env.n_double == 1 {
                        tmp = 12.36;
                    } else if env.n_h == 0 && chg == 0 && env.n_triple == 1 && env.n_double == 1 {
                        tmp = 13.60;
                    } else if env.n_h == 1 && chg == 0 && env.n_single == 2 && in3 {
                        tmp = 21.94;
                    } else if env.n_h == 1 && chg == 0 && env.n_single == 2 && !in3 {
                        tmp = 12.03;
                    } else if env.n_h == 0 && chg == 1 && env.n_triple == 1 && env.n_single == 1 {
                        tmp = 4.36;
                    } else if env.n_h == 1 && chg == 1 && env.n_double == 1 && env.n_single == 1 {
                        tmp = 13.97;
                    } else if env.n_h == 2 && chg == 1 && env.n_single == 2 {
                        tmp = 16.61;
                    } else if env.n_h == 0 && chg == 0 && env.n_aromatic == 2 {
                        tmp = 12.89;
                    } else if env.n_h == 1 && chg == 0 && env.n_aromatic == 2 {
                        tmp = 15.79;
                    } else if env.n_h == 1 && chg == 1 && env.n_aromatic == 2 {
                        tmp = 14.14;
                    }
                }
                3 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 3 && in3 {
                        tmp = 3.01;
                    } else if env.n_h == 0 && chg == 0 && env.n_single == 3 && !in3 {
                        tmp = 3.24;
                    } else if env.n_h == 0 && chg == 0 && env.n_single == 1 && env.n_double == 2 {
                        tmp = 11.68;
                    } else if env.n_h == 0 && chg == 1 && env.n_single == 2 && env.n_double == 1 {
                        tmp = 3.01;
                    } else if env.n_h == 1 && chg == 1 && env.n_single == 3 {
                        tmp = 4.44;
                    } else if env.n_h == 0 && chg == 0 && env.n_aromatic == 3 {
                        tmp = 4.41;
                    } else if env.n_h == 0 && chg == 0 && env.n_single == 1 && env.n_aromatic == 2 {
                        tmp = 4.93;
                    } else if env.n_h == 0 && chg == 0 && env.n_double == 1 && env.n_aromatic == 2 {
                        tmp = 8.39;
                    } else if env.n_h == 0 && chg == 1 && env.n_aromatic == 3 {
                        tmp = 4.10;
                    } else if env.n_h == 0 && chg == 1 && env.n_single == 1 && env.n_aromatic == 2 {
                        tmp = 3.88;
                    }
                }
                4 => {
                    if env.n_h == 0 && env.n_single == 4 && chg == 1 {
                        tmp = 0.0;
                    }
                }
                _ => {}
            }
            if tmp < 0.0 {
                tmp = 30.5 - env.n_neighbors as f64 * 8.2 + env.n_h as f64 * 1.5;
                if tmp < 0.0 {
                    tmp = 0.0;
                }
            }
        }
        Element::O => {
            match env.n_neighbors {
                1 => {
                    if env.n_h == 0 && chg == 0 && env.n_double == 1 {
                        tmp = 17.07;
                    } else if env.n_h == 1 && chg == 0 && env.n_single == 1 {
                        tmp = 20.23;
                    } else if env.n_h == 0 && chg == -1 && env.n_single == 1 {
                        tmp = 23.06;
                    }
                }
                2 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 2 && in3 {
                        tmp = 12.53;
                    } else if env.n_h == 0 && chg == 0 && env.n_single == 2 && !in3 {
                        tmp = 9.23;
                    } else if env.n_h == 0 && chg == 0 && env.n_aromatic == 2 {
                        tmp = 13.14;
                    }
                }
                _ => {}
            }
            if tmp < 0.0 {
                tmp = 28.5 - env.n_neighbors as f64 * 8.6 + env.n_h as f64 * 1.5;
                if tmp < 0.0 {
                    tmp = 0.0;
                }
            }
        }
        _ => return 0.0,
    }
    tmp
}

/// Additional PSA contribution for S and P atoms (only used when `no_only=false`).
fn psa_contrib_sp(mol: &Molecule, idx: petgraph::graph::NodeIndex) -> f64 {
    let atom = mol.atom(idx);
    let env = atom_env(mol, idx);
    let chg = atom.charge as i32;

    match atom.element {
        Element::P => {
            match env.n_neighbors {
                2 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 1 && env.n_double == 1 {
                        return 34.14;
                    }
                }
                3 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 3 {
                        return 13.59;
                    } else if env.n_h == 1 && chg == 0 && env.n_single == 2 && env.n_double == 1 {
                        return 23.47;
                    }
                }
                4 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 3 && env.n_double == 1 {
                        return 9.81;
                    }
                }
                _ => {}
            }
            0.0
        }
        Element::S => {
            match env.n_neighbors {
                1 => {
                    if env.n_h == 0 && chg == 0 && env.n_double == 1 {
                        return 32.09;
                    } else if env.n_h == 1 && chg == 0 && env.n_single == 1 {
                        return 38.80;
                    }
                }
                2 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 2 {
                        return 25.30;
                    } else if env.n_h == 0 && chg == 0 && env.n_aromatic == 2 {
                        return 28.24;
                    }
                }
                3 => {
                    if env.n_h == 0 && chg == 0 && env.n_aromatic == 2 && env.n_double == 1 {
                        return 21.70;
                    } else if env.n_h == 0 && chg == 0 && env.n_single == 2 && env.n_double == 1 {
                        return 19.21;
                    }
                }
                4 => {
                    if env.n_h == 0 && chg == 0 && env.n_single == 2 && env.n_double == 2 {
                        return 8.38;
                    }
                }
                _ => {}
            }
            0.0
        }
        _ => 0.0,
    }
}

/// Topological polar surface area.
///
/// `no_only = true` reproduces upstream `TopoPSA(NO)` (Ertl TPSA over N and O
/// only). `no_only = false` reproduces upstream `TopoPSA`, which adds the
/// Ertl S and P contributions.
pub struct TopoPSA {
    no_only: bool,
}

impl TopoPSA {
    pub fn new(no_only: bool) -> Self {
        Self { no_only }
    }
}

impl Descriptor for TopoPSA {
    fn name(&self) -> &str {
        if self.no_only {
            "TopoPSA(NO)"
        } else {
            "TopoPSA"
        }
    }

    fn description(&self) -> &str {
        if self.no_only {
            "Topological polar surface area (N, O only)"
        } else {
            "Topological polar surface area (N, O, S, P)"
        }
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let mut total = 0.0;
        for (idx, atom) in mol.atoms() {
            match atom.element {
                Element::N | Element::O => total += psa_contrib_no(mol, idx),
                Element::S | Element::P if !self.no_only => total += psa_contrib_sp(mol, idx),
                _ => {}
            }
        }
        Ok(total)
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
    fn topo_psa_hexane_zero() {
        let mol = parse_smiles("CCCCCC").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(approx(v, 0.0, 1e-6), "TopoPSA(NO) Hexane: got {v}");
        let v = TopoPSA::new(false).calculate(&mol).unwrap();
        assert!(approx(v, 0.0, 1e-6), "TopoPSA Hexane: got {v}");
    }

    #[test]
    fn topo_psa_caffeine() {
        // Upstream reference: 61.82 (both NO-only and full).
        let mol = parse_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(approx(v, 61.82, 1e-2), "TopoPSA(NO) Caffeine: got {v}");
        let v = TopoPSA::new(false).calculate(&mol).unwrap();
        assert!(approx(v, 61.82, 1e-2), "TopoPSA Caffeine: got {v}");
    }

    #[test]
    fn topo_psa_epicatechin() {
        // 5 OH + 1 ring O — upstream 110.38.
        let mol = parse_smiles("C1[C@H]([C@H](OC2=CC(=CC(=C21)O)O)C3=CC(=C(C=C3)O)O)O").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(approx(v, 110.38, 1e-2), "TopoPSA(NO) Epicatechin: got {v}");
    }

    #[test]
    fn topo_psa_dmso_sulfur() {
        // DMSO: S=O contributes 32.09 (S, 1 neighbour, double bond), plus the
        // sulfur's =O contributes 17.07 from the oxygen → 36.28 total when
        // no_only=false. Upstream NO-only=17.07.
        let mol = parse_smiles("CS(=O)C").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(approx(v, 17.07, 1e-2), "TopoPSA(NO) DMSO: got {v}");
        let v = TopoPSA::new(false).calculate(&mol).unwrap();
        // S has 3 neighbours (2 C single + 1 O double), nH=0, chg=0:
        // matches "n_single == 2 && n_double == 1" → 19.21; +17.07 from O = 36.28.
        assert!(approx(v, 36.28, 1e-2), "TopoPSA DMSO: got {v}");
    }

    #[test]
    fn topo_psa_thiophene_aromatic_sulfur() {
        // Thiophene: NO-only=0, full=28.24 (aromatic S, 2 aromatic neighbours).
        let mol = parse_smiles("c1cscc1").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(approx(v, 0.0, 1e-2), "TopoPSA(NO) Thiophene: got {v}");
        let v = TopoPSA::new(false).calculate(&mol).unwrap();
        assert!(approx(v, 28.24, 1e-2), "TopoPSA Thiophene: got {v}");
    }

    #[test]
    fn topo_psa_methylphosphonic_acid() {
        // CP(=O)(O)O — P with 3 single + 1 double, nH=0 → 9.81; 2 OH (20.23·2=40.46)
        // + 1 =O (17.07). Total = 9.81 + 40.46 + 17.07 = 67.34 (upstream).
        // NO-only = 40.46 + 17.07 = 57.53 (upstream).
        let mol = parse_smiles("CP(=O)(O)O").unwrap();
        let v = TopoPSA::new(true).calculate(&mol).unwrap();
        assert!(
            approx(v, 57.53, 1e-2),
            "TopoPSA(NO) MethylphosphonicAcid: got {v}"
        );
        let v = TopoPSA::new(false).calculate(&mol).unwrap();
        assert!(
            approx(v, 67.34, 1e-2),
            "TopoPSA MethylphosphonicAcid: got {v}"
        );
    }
}
