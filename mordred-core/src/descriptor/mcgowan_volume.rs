/// McGowan characteristic volume.
///
/// `VMcGowan = Σ V_atomic_i - 6.56 · n_bonds_total`, where bonds_total counts
/// every bond including implicit-hydrogen bonds. Per-element atomic volumes
/// are the McGowan-Abraham characteristic-volume table (cm³/mol).
///
/// Reference: McGowan & Mellors (1986); upstream `mordred/McGowanVolume.py`.
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::{Element, Molecule};

fn atomic_volume(elem: Element) -> Option<f64> {
    Some(match elem {
        Element::H => 8.71,
        Element::B => 18.32,
        Element::C => 16.35,
        Element::N => 14.39,
        Element::O => 12.43,
        Element::F => 10.48,
        Element::Si => 26.83,
        Element::P => 24.87,
        Element::S => 22.91,
        Element::Cl => 20.95,
        Element::Br => 26.21,
        Element::I => 34.53,
        Element::Se => 28.73,
        // Less common elements not in McGowan's original table.
        _ => return None,
    })
}

pub struct McGowanVolume;

impl Descriptor for McGowanVolume {
    fn name(&self) -> &str {
        "VMcGowan"
    }
    fn description(&self) -> &str {
        "McGowan characteristic volume (cm³/mol)"
    }
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let props = mol.properties();
        let mut sum = 0.0;
        for (_, atom) in mol.atoms() {
            sum += atomic_volume(atom.element).ok_or_else(|| {
                MordredError::CalculationError(format!(
                    "no McGowan volume for element {:?}",
                    atom.element
                ))
            })?;
            // Implicit hydrogens each contribute H's volume.
            sum += f64::from(atom.implicit_h) * 8.71;
        }
        let total_bonds = mol.bond_count() as f64 + f64::from(props.implicit_h_sum);
        Ok(sum - 6.56 * total_bonds)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_hexane() {
        // Hexane C6H14: 6·16.35 + 14·8.71 - 6.56·19 = 95.4
        let mol = parse_smiles("CCCCCC").unwrap();
        let v = McGowanVolume.calculate(&mol).unwrap();
        assert!((v - 95.4).abs() < 1e-6, "got {v}");
    }

    #[test]
    fn test_benzene() {
        // C6H6: 6·16.35 + 6·8.71 - 6.56·12 = 71.64
        let mol = parse_smiles("c1ccccc1").unwrap();
        let v = McGowanVolume.calculate(&mol).unwrap();
        assert!((v - 71.64).abs() < 1e-6, "got {v}");
    }

    #[test]
    fn test_caffeine() {
        let mol = parse_smiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C").unwrap();
        let v = McGowanVolume.calculate(&mol).unwrap();
        assert!((v - 136.32).abs() < 1e-4, "got {v}");
    }
}
