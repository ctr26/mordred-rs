# mordred-rs Roadmap

A Rust reimplementation of [mordred](https://github.com/mordred-descriptor/mordred), addressing upstream milestone 2.0.0: **remove RDKit dependency**.

## Status

### Phase 1 ✅ Complete

| Component | Status |
|-----------|--------|
| Cargo workspace | ✅ core, cli, py |
| SMILES parser | ✅ organic subset, brackets, aromaticity, branches, rings, charges |
| Molecule graph | ✅ petgraph with distance matrix |
| Python bindings | ✅ PyO3 Calculator class |
| CLI | ✅ clap, JSON/CSV output |

**Descriptors implemented (31):**
- Constitutional (4): `nAtom`, `nHeavyAtom`, `nBond`, `MW`
- Topological (3): `WienerIndex`, `Zagreb1`, `Zagreb2`
- Connectivity (2): `Chi0`, `Chi1`
- AtomCount (7): `nC`, `nN`, `nO`, `nS`, `nP`, `nHalo`, `nHetero`
- BondCount (4): `nSingle`, `nDouble`, `nTriple`, `nAromatic`
- RingCount (11): `nRing`, `nRing3`–`nRing12`

---

## Phase 2: Core 2D Descriptors

Priority descriptors that don't require RDKit or 3D coordinates.

### From upstream issues

| Issue | Descriptor | Reference | Status |
|-------|------------|-----------|--------|
| [#30](https://github.com/mordred-descriptor/mordred/issues/30) | LogP (XLogP3, JPLogP) | [doi:10.1021/ci700257y](https://doi.org/10.1021/ci700257y) | TODO |
| [#26](https://github.com/mordred-descriptor/mordred/issues/26) | Wiener Polarity Index | [arxiv:1801.05963](https://arxiv.org/pdf/1801.05963.pdf) | TODO |
| [#15](https://github.com/mordred-descriptor/mordred/issues/15) | Gasteiger charge (rdkit-free) | — | TODO |

### Additional 2D descriptors from upstream mordred

| Module | Descriptors | Status |
|--------|-------------|--------|
| Aromatic | nAromatic, AromaticRatio | TODO |
| RingCount | nRing, nRing3-12, nHeteroRing | ✅ nRing + nRing3-12 done; nHeteroRing TODO |
| AtomCount | nC, nN, nO, nS, nP, nHalo, nHetero | ✅ |
| BondCount | nSingle, nDouble, nTriple, nAromatic | ✅ |
| RotatableBond | nRot | TODO |
| HydrogenBond | nHBDon, nHBAcc | TODO |
| Lipinski | Lipinski filter (MW, HBD, HBA, LogP) | TODO |
| CarbonTypes | nSpC, nSp2C, nSp3C | TODO |
| BalabanJ | J index | TODO |
| BertzCT | Bertz complexity | TODO |
| KappaShapeIndex | Kappa1, Kappa2, Kappa3 | TODO |
| PathCount | nPath1-10 | TODO |
| EccentricConnectivityIndex | ECI | TODO |
| TopoPSA | Topological polar surface area | TODO |
| McGowanVolume | McGowan characteristic volume | TODO |

### Infrastructure needed

- [x] SSSR algorithm (smallest set of smallest rings)
- [ ] Atom typing system (for LogP, EState)
- [x] Hybridization detection from bond orders (bond order data available)

---

## Phase 3: EState & Autocorrelation

| Module | Descriptors | Notes |
|--------|-------------|-------|
| EState | EState indices, SaaX counts | Requires intrinsic state calculation |
| SLogP | Crippen LogP, MR | Atom contribution tables |
| Autocorrelation | ATS, AATS, ATSC, MATS, GATS | 2D autocorrelation |
| Chi | Chi2-Chi10 | Higher-order connectivity |
| InformationContent | IC, SIC, CIC, BIC | Shannon entropy based |

Upstream issue: [#63](https://github.com/mordred-descriptor/mordred/issues/63) — AtomTypeEState indices should return zero not empty.

---

## Phase 4: Production

| Task | Notes | Status |
|------|-------|--------|
| GitHub repo | Push to ctr26/mordred-rs | ✅ |
| CI | GitHub Actions (cargo test, clippy, maturin) | ✅ |
| Branch protection | Squash-only merges, dev branch strategy | ✅ |
| Benchmarks | criterion.rs vs Python mordred | TODO |
| SDF parser | [#61](https://github.com/mordred-descriptor/mordred/issues/61) — PubChem SDF hardening | TODO |
| Rayon | Parallel batch calculation | TODO |
| PyPI | maturin publish | ✅ (configured) |
| crates.io | cargo publish | ✅ (configured) |
| Docs | [#8](https://github.com/mordred-descriptor/mordred/issues/8) — numerical examples | TODO |

---

## Phase 5: 3D Descriptors (Future)

Requires 3D coordinates (conformer generation or external input).

From upstream milestone 1.4.0:

| Issue | Descriptor |
|-------|------------|
| [#14](https://github.com/mordred-descriptor/mordred/issues/14) | WHIM |
| [#13](https://github.com/mordred-descriptor/mordred/issues/13) | RDF |
| [#12](https://github.com/mordred-descriptor/mordred/issues/12) | 3D autocorrelation |

Additional 3D modules from upstream:
- CPSA (charged partial surface area)
- MoRSE
- GravitationalIndex
- MomentOfInertia
- GeometricalIndex
- PBF

---

## Phase 6: Advanced (Future)

| Issue | Feature |
|-------|---------|
| [#70](https://github.com/mordred-descriptor/mordred/issues/70) | QM descriptors |
| [#49](https://github.com/mordred-descriptor/mordred/issues/49) | Molecular framework export |
| — | BCUT (Burden eigenvalues) |

---

## Upstream Mordred Issues Addressed

| Issue | Title | How mordred-rs addresses it |
|-------|-------|----------------------------|
| [#15](https://github.com/mordred-descriptor/mordred/issues/15) | Split backend (remove rdkit) | **Core goal** — pure Rust, no RDKit |
| [#113](https://github.com/mordred-descriptor/mordred/issues/113) | numpy.product deprecated | N/A — no numpy |
| [#107](https://github.com/mordred-descriptor/mordred/issues/107) | numpy.float deprecated | N/A — no numpy |
| [#84](https://github.com/mordred-descriptor/mordred/issues/84) | networkx compatibility | N/A — petgraph instead |
| [#82](https://github.com/mordred-descriptor/mordred/issues/82) | Multiprocessing overloads cores | Rayon with work-stealing |
| [#19](https://github.com/JacksonBurns/mordred-community/issues/19) | Speed up calculations | Rust = 10-100x faster |

---

## References

- Upstream: https://github.com/mordred-descriptor/mordred
- Community fork: https://github.com/JacksonBurns/mordred-community
- Paper: Moriwaki et al. (2018) J Cheminform 10:4 [doi:10.1186/s13321-018-0258-y](https://doi.org/10.1186/s13321-018-0258-y)
