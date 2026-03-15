# Contributing to mordred-rs

## Development Setup

### Prerequisites

- Rust 1.87+ (install via [rustup](https://rustup.rs/))
- Python 3.9+ (for Python bindings)
- maturin (`pip install maturin`)

### Building

```bash
# Build all crates
make build

# Run tests
make test

# Run all checks (required before pushing)
make prepush
```

### Code Style

We follow Google style guides:

- **Rust**: Google C++ style principles adapted for Rust, enforced via `rustfmt` and `clippy`
- **Python**: Google Python style guide, enforced via `ruff`
- **C++**: Google C++ style guide for the FFI wrapper

### Commit Messages

We use [Conventional Commits](https://www.conventionalcommits.org/):

```
<type>(<scope>): <description>

[optional body]

[optional footer]
```

**Types**: `feat`, `fix`, `docs`, `style`, `refactor`, `perf`, `test`, `build`, `ci`, `chore`, `revert`

**Scopes**: `core`, `cli`, `py`, `ffi`

**Examples**:
```
feat(core): add SSSR ring detection algorithm
fix(py): handle RDKit Mol objects without GetSmiles method
docs: update README with Python API examples
test(core): add benzene ring count validation
refactor(ffi): simplify result access functions
```

### Branch Naming

Use conventional branch names:

```
feat/<description>
fix/<description>
refactor/<description>
docs/<description>
```

### Pull Requests

1. Create a feature branch: `git checkout -b feat/my-feature`
2. Make changes following the style guides
3. Run `make prepush` to verify
4. Commit with conventional commit messages
5. Push and open a PR

### Adding New Descriptors

1. Create a new file in `mordred-core/src/descriptor/`
2. Implement the `Descriptor` trait for each descriptor
3. Add a registration method in `descriptor/mod.rs`
4. Call it in `DescriptorSet::all()`
5. Add unit tests in the same file
6. Add integration tests in `tests/integration.rs`
7. Verify descriptor names match the original Python mordred

### Running Python Bindings Locally

```bash
make install-python
python -c "from mordred import Calculator, descriptors; calc = Calculator(descriptors); print(calc('CCO'))"
```
