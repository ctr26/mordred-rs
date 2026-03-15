/// Ring count descriptors.
///
/// Reference: mordred Python `RingCount.py`
use crate::descriptor::Descriptor;
use crate::error::MordredError;
use crate::molecule::Molecule;

/// Total number of rings in the SSSR.
pub struct RingCount;

impl Descriptor for RingCount {
    fn name(&self) -> &str {
        "nRing"
    }

    fn description(&self) -> &str {
        "Number of rings (SSSR)"
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        Ok(mol.num_rings() as f64)
    }
}

/// Number of rings of a specific size.
pub struct RingSizeCount {
    size: usize,
    name: String,
    description: String,
}

impl RingSizeCount {
    pub fn new(size: usize) -> Self {
        Self {
            size,
            name: format!("n{}Ring", size),
            description: format!("Number of {}-membered rings", size),
        }
    }
}

impl Descriptor for RingSizeCount {
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &str {
        &self.description
    }

    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError> {
        let count = mol.ring_info().rings_of_size(self.size).len();
        Ok(count as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::parse_smiles;

    #[test]
    fn test_benzene_ring_count() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(RingCount.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_naphthalene_ring_count() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        assert_eq!(RingCount.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_ethane_ring_count() {
        let mol = parse_smiles("CC").unwrap();
        assert_eq!(RingCount.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_benzene_ring_size_6() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let desc = RingSizeCount::new(6);
        assert_eq!(desc.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_benzene_ring_size_5() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        let desc = RingSizeCount::new(5);
        assert_eq!(desc.calculate(&mol).unwrap(), 0.0);
    }

    #[test]
    fn test_cyclopropane_ring_size_3() {
        let mol = parse_smiles("C1CC1").unwrap();
        let desc = RingSizeCount::new(3);
        assert_eq!(desc.calculate(&mol).unwrap(), 1.0);
    }

    #[test]
    fn test_naphthalene_ring_size_6() {
        let mol = parse_smiles("c1ccc2ccccc2c1").unwrap();
        let desc = RingSizeCount::new(6);
        assert_eq!(desc.calculate(&mol).unwrap(), 2.0);
    }

    #[test]
    fn test_descriptor_name() {
        assert_eq!(RingCount.name(), "nRing");
        let desc = RingSizeCount::new(5);
        assert_eq!(desc.name(), "n5Ring");
    }
}
