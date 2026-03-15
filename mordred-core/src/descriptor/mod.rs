pub mod atom_count;
pub mod bond_count;
pub mod connectivity;
pub mod constitutional;
pub mod ring_count;
pub mod topological;

use crate::error::MordredError;
use crate::molecule::Molecule;

/// A molecular descriptor that computes a single numeric value.
pub trait Descriptor: Send + Sync {
    /// Short name of the descriptor (e.g. "MW", "nAtom").
    fn name(&self) -> &str;

    /// Human-readable description.
    fn description(&self) -> &str;

    /// Calculate the descriptor value for a molecule.
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError>;
}

/// A collection of descriptors for batch calculation.
pub struct DescriptorSet {
    descriptors: Vec<Box<dyn Descriptor>>,
}

impl DescriptorSet {
    /// Create a new empty descriptor set.
    pub fn new() -> Self {
        Self {
            descriptors: Vec::new(),
        }
    }

    /// Create a descriptor set with all built-in descriptors.
    pub fn all() -> Self {
        let mut set = Self::new();
        set.add_all_constitutional();
        set.add_all_topological();
        set.add_all_connectivity();
        set.add_all_atom_counts();
        set.add_all_bond_counts();
        set.add_all_ring_counts();
        set
    }

    /// Add a descriptor.
    pub fn add(&mut self, descriptor: Box<dyn Descriptor>) {
        self.descriptors.push(descriptor);
    }

    /// Add all constitutional descriptors.
    pub fn add_all_constitutional(&mut self) {
        self.add(Box::new(constitutional::AtomCount));
        self.add(Box::new(constitutional::HeavyAtomCount));
        self.add(Box::new(constitutional::BondCount));
        self.add(Box::new(constitutional::MolecularWeight));
    }

    /// Add all topological descriptors.
    pub fn add_all_topological(&mut self) {
        self.add(Box::new(topological::WienerIndex));
        self.add(Box::new(topological::ZagrebIndex1));
        self.add(Box::new(topological::ZagrebIndex2));
    }

    /// Add all connectivity descriptors.
    pub fn add_all_connectivity(&mut self) {
        self.add(Box::new(connectivity::Chi0));
        self.add(Box::new(connectivity::Chi1));
    }

    /// Add all element-specific atom count descriptors.
    pub fn add_all_atom_counts(&mut self) {
        self.add(Box::new(atom_count::CarbonCount));
        self.add(Box::new(atom_count::NitrogenCount));
        self.add(Box::new(atom_count::OxygenCount));
        self.add(Box::new(atom_count::SulfurCount));
        self.add(Box::new(atom_count::PhosphorusCount));
        self.add(Box::new(atom_count::HalogenCount));
        self.add(Box::new(atom_count::HeteroatomCount));
    }

    /// Add all ring count descriptors.
    pub fn add_all_ring_counts(&mut self) {
        self.add(Box::new(ring_count::RingCount));
        for size in 3..=12 {
            self.add(Box::new(ring_count::RingSizeCount::new(size)));
        }
    }

    /// Add all bond type count descriptors.
    pub fn add_all_bond_counts(&mut self) {
        self.add(Box::new(bond_count::SingleBondCount));
        self.add(Box::new(bond_count::DoubleBondCount));
        self.add(Box::new(bond_count::TripleBondCount));
        self.add(Box::new(bond_count::AromaticBondCount));
    }

    /// Calculate all descriptors for a molecule.
    pub fn calculate(&self, mol: &Molecule) -> Vec<(&str, Result<f64, MordredError>)> {
        self.descriptors
            .iter()
            .map(|d| (d.name(), d.calculate(mol)))
            .collect()
    }

    /// List all descriptor names.
    pub fn names(&self) -> Vec<&str> {
        self.descriptors.iter().map(|d| d.name()).collect()
    }

    /// List all descriptors with descriptions.
    pub fn list(&self) -> Vec<(&str, &str)> {
        self.descriptors
            .iter()
            .map(|d| (d.name(), d.description()))
            .collect()
    }

    /// Number of descriptors in the set.
    pub fn len(&self) -> usize {
        self.descriptors.len()
    }

    /// Returns true if the set is empty.
    pub fn is_empty(&self) -> bool {
        self.descriptors.is_empty()
    }
}

impl Default for DescriptorSet {
    fn default() -> Self {
        Self::all()
    }
}
