# mordred-rs

A fast molecular descriptor calculator in Rust — port of [mordred](https://github.com/mordred-descriptor/mordred).

## Features

- **9 molecular descriptors** across 3 categories (constitutional, topological, connectivity)
- **SMILES parser** with support for organic subset, bracket atoms, branches, ring closures, aromaticity
- **Python bindings** via PyO3/maturin
- **CLI tool** with JSON and CSV output

## Quick Start

### CLI

```bash
# Calculate all descriptors for a molecule
cargo run -p mordred-cli -- "CCO"

# List available descriptors
cargo run -p mordred-cli -- --list

# Process a file of SMILES
cargo run -p mordred-cli -- --file molecules.smi --format csv --output results.csv
```

### Rust Library

```rust
use mordred_core::{DescriptorSet, parse_smiles};

let mol = parse_smiles("c1ccccc1").unwrap(); // benzene
let descriptors = DescriptorSet::all();

for (name, result) in descriptors.calculate(&mol) {
    match result {
        Ok(val) => println!("{}: {:.4}", name, val),
        Err(e) => println!("{}: error ({})", name, e),
    }
}
```

### Python

```bash
cd mordred-py && maturin develop
```

```python
from mordred_rs import Calculator

calc = Calculator()
result = calc.calculate("CCO")
print(result)  # {'nAtom': 9.0, 'MW': 46.069, ...}

# Batch calculation
results = calc.calculate_batch(["CCO", "c1ccccc1", "CC(=O)O"])
```

## Available Descriptors

| Name | Category | Description |
|------|----------|-------------|
| nAtom | Constitutional | Total atom count (including implicit H) |
| nHeavyAtom | Constitutional | Heavy atom count (non-hydrogen) |
| nBond | Constitutional | Bond count (between explicit atoms) |
| MW | Constitutional | Molecular weight (average isotope) |
| WienerIndex | Topological | Sum of shortest path distances |
| Zagreb1 | Topological | First Zagreb index (sum of squared degrees) |
| Zagreb2 | Topological | Second Zagreb index (sum of degree products over edges) |
| Chi0 | Connectivity | Zero-order Randic connectivity index |
| Chi1 | Connectivity | First-order Randic connectivity index |

## Project Structure

```
mordred-rs/
├── mordred-core/    # Core library: molecule, descriptors
├── mordred-py/      # Python bindings (PyO3 + maturin)
└── mordred-cli/     # Command-line interface
```

## Building

```bash
cargo build --workspace --exclude mordred-py
cargo test --workspace --exclude mordred-py
```

For Python bindings:

```bash
cd mordred-py && maturin develop --release
```

## License

MIT
