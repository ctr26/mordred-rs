# Architecture

## Overview

mordred-rs is a Cargo workspace with three crates:

```
mordred-rs/
├── mordred-core/       # Core library
│   └── src/
│       ├── molecule/   # Molecular representation
│       │   ├── element.rs   # Periodic table data
│       │   ├── atom.rs      # Atom struct
│       │   ├── bond.rs      # Bond types
│       │   ├── mol.rs       # Molecule graph (petgraph)
│       │   └── smiles.rs    # SMILES parser
│       ├── descriptor/ # Descriptor system
│       │   ├── mod.rs            # Descriptor trait + DescriptorSet
│       │   ├── constitutional.rs # AtomCount, BondCount, MW, HeavyAtomCount
│       │   ├── topological.rs    # WienerIndex, Zagreb indices
│       │   └── connectivity.rs   # Chi0, Chi1 (Randic indices)
│       ├── error.rs    # Error types
│       └── lib.rs      # Public API
├── mordred-py/         # Python bindings
│   └── src/lib.rs      # PyO3 Calculator class
└── mordred-cli/        # CLI
    └── src/main.rs     # clap-based CLI
```

## Design Decisions

### Molecule as a Graph

The `Molecule` struct wraps `petgraph::UnGraph<Atom, Bond>`. This gives us:
- O(1) atom/bond access
- Built-in graph algorithms (Floyd-Warshall for distance matrix)
- Easy iteration over neighbours, edges

### Descriptor Trait

```rust
pub trait Descriptor: Send + Sync {
    fn name(&self) -> &str;
    fn description(&self) -> &str;
    fn calculate(&self, mol: &Molecule) -> Result<f64, MordredError>;
}
```

Each descriptor is a zero-sized struct implementing this trait. `DescriptorSet` collects them as `Box<dyn Descriptor>` for batch calculation.

The `Send + Sync` bounds allow future parallelisation with rayon.

### SMILES Parser

A hand-written recursive descent parser that handles:
- Organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I)
- Bracket atoms with charges, isotopes, H counts (`[Fe]`, `[NH4+]`)
- Aromatic atoms (c, n, o, s, p)
- Bond types (single, double, triple, aromatic)
- Branches with `(` and `)`
- Ring closures (single digit and `%nn`)
- Disconnected fragments (`.`)

Implicit hydrogens are filled after parsing based on default valence rules, with a correction for aromatic atoms (which have an extra pi bond).

### Error Handling

`MordredError` (via thiserror) covers SMILES parsing failures, calculation errors, and invalid molecules. Descriptors return `Result<f64, MordredError>` so callers can handle partial failures gracefully.

## Adding a New Descriptor

1. Create a struct in the appropriate category module
2. Implement `Descriptor` for it
3. Register it in `DescriptorSet::add_all_*` methods
4. Add unit tests + integration test entries
