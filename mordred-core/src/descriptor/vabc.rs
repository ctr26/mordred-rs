/// van der Waals volume (ABC) — Zhao 2003.
///
/// `Vabc = Σ (4/3·π·r_i³) - 5.92·Nb - 14.7·NRa - 3.8·NRA`
///
/// where:
/// * `r_i` is the Bondi van der Waals radius (Å) of atom `i` — implicit
///   hydrogens are summed explicitly with `r_H = 1.20`;
/// * `Nb` is the total bond count, including bonds to implicit hydrogens
///   (matches upstream `BondCount(type="any")` with `explicit_hydrogens=True`);
/// * `NRa` is the number of SSSR rings whose atoms are all aromatic;
/// * `NRA` is the number of SSSR rings with at least one non-aromatic atom.
///
/// Per-atom contribution is computed at runtime from the exact Bondi-radius
/// table copied verbatim from upstream `mordred/VdwVolumeABC.py` to avoid
/// rounding drift (rounded constants caused a previous parity miss of ~0.03).
///
/// Reference: Zhao, Abraham & Zissimos, J. Org. Chem. 68 (2003) 7368-7373,
/// :doi:`10.1021/jo034808o`. Upstream source: `mordred/VdwVolumeABC.py`.
use std::f64::consts::PI;

use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{Element, Molecule};

/// Bondi van der Waals radii (Å), copied verbatim from upstream
/// `mordred/VdwVolumeABC.py::bondi_radii`. The per-atom volume is then
/// `4/3·π·r³` (no rounding) to reproduce upstream's `atom_contrib` dict.
fn bondi_radius(elem: Element) -> Option<f64> {
    Some(match elem {
        Element::H => 1.20,
        Element::C => 1.70,
        Element::N => 1.55,
        Element::O => 1.52,
        Element::F => 1.47,
        Element::Cl => 1.75,
        Element::Br => 1.85,
        Element::P => 1.80,
        Element::S => 1.80,
        // "As" => 1.85 — upstream supports it, but not in mordred-rs Element enum.
        Element::B => 2.13,
        Element::Si => 2.10,
        Element::Se => 1.90,
        _ => return None,
    })
}

#[inline]
fn atom_contrib(elem: Element) -> Option<f64> {
    bondi_radius(elem).map(|r| 4.0 / 3.0 * PI * r * r * r)
}

pub struct Vabc;

impl Descriptor for Vabc {
    fn name(&self) -> &str {
        "Vabc"
    }
    fn description(&self) -> &str {
        "ABC van der Waals volume"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        // Cache implicit-H contribution since it's by far the hottest term.
        let v_h = atom_contrib(Element::H).expect("H always has a Bondi radius");

        // Σ heavy-atom contributions + implicit-H contributions.
        let mut ac = 0.0;
        for (_, atom) in mol.atoms() {
            ac += atom_contrib(atom.element).ok_or_else(|| {
                MordredError::CalculationError(format!(
                    "no Bondi vdW radius for element {:?}",
                    atom.element
                ))
            })?;
            ac += f64::from(atom.implicit_h) * v_h;
        }

        // Total bond count includes implicit-H bonds (matches BondCount("any")).
        let n_bonds = mol.bond_count() as u32 + mol.properties().implicit_h_sum;

        // Split SSSR rings into aromatic (all atoms aromatic) and non-aromatic.
        let rings = mol.ring_info().rings();
        let mut n_ra = 0u32; // aromatic
        let mut n_rna = 0u32; // non-aromatic
        for ring in rings {
            if ring.iter().all(|&i| mol.atom(i).is_aromatic) {
                n_ra += 1;
            } else {
                n_rna += 1;
            }
        }

        Ok(ac - 5.92 * f64::from(n_bonds) - 14.7 * f64::from(n_ra) - 3.8 * f64::from(n_rna))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    // Upstream tolerance: digit=4 → |Δ| ≤ 1.5e-4.
    const TOL: f64 = 1.5e-4;

    fn vabc(smiles: &str) -> f64 {
        let mol = parse_smiles(smiles).unwrap();
        Vabc.calculate(&mol).unwrap()
    }

    #[test]
    fn vabc_hexane() {
        // CCCCCC: 6C + 14H, 19 total bonds, 0 rings → 112.3323665
        let v = vabc("CCCCCC");
        assert!((v - 112.3323665).abs() < TOL, "hexane: {v}");
    }

    #[test]
    fn vabc_benzene() {
        // c1ccccc1: 6C + 6H, 12 total bonds, 1 aromatic ring → 81.16653165
        let v = vabc("c1ccccc1");
        assert!((v - 81.16653165).abs() < TOL, "benzene: {v}");
    }

    #[test]
    fn vabc_caffeine() {
        // CN1C=NC2=C1C(=O)N(C(=O)N2C)C → 151.43307774
        let v = vabc("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
        assert!((v - 151.43307774).abs() < TOL, "caffeine: {v}");
    }
}
