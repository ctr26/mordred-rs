# NEXT — mordred-rs

> Resumable plan. Updated **2026-06-13** after audit + fresh `cargo test` run.
> Run `git status`, `git log --oneline -10`, and the parity gate first.

## What this is
Rust reimplementation of the `mordred` molecular-descriptor library, exposed as
a **Python API-compatible shim** (drop-in for the original). Scope: **full
upstream parity** including 3D, BCUT, QM. Pitch target: `mordred-community`.

## Parity gate (the source of truth)
Vendored 8293-case corpus runs the original mordred's pytest against this impl.
**Current state (2026-06-13): 624 pass / 1 fail / 294 skip / 7374 xfail.**

The single remaining failure:
- `InformationContent.yaml::EllagicAcid::BIC0` — expected 0.282, got 0.2759
  (tolerance 0.0015, off by 0.0061 absolute, ~2.1% relative).

## Diagnosis of the EllagicAcid BIC0 failure

**Root cause (most likely):** the Hückel-via-fixpoint aromaticity perception in
`mordred-core/src/molecule/aromaticity.rs::perceive_aromaticity` over-marks the
two **lactone 6-rings** of ellagic acid as aromatic. RDKit (which upstream
mordred wraps) treats only the two benzene rings as aromatic; the lactone
rings stay Kekulized.

### Triangulation
- EllagicAcid formula C₁₄H₆O₈ → 14 C + 8 O + 6 H = 28 atoms
- IC0 (entropy of atom classes) = -0.5·log₂(0.5) - 0.286·log₂(0.286) -
  0.214·log₂(0.214) ≈ 1.493 — robust to aromaticity choice
- Solving `BIC0 = IC0 / log₂(B)`:
  - Expected 0.282 ⇒ B ≈ 39.3
  - Current  0.276 ⇒ B ≈ 42.5
- With all 19 ring bonds Aromatic (1.5): explicit Σ ≈ 36.5 + implicit 6H = **42.5** ✅ matches current
- With only the 12 benzene ring bonds Aromatic and the 7 lactone ring bonds
  Kekulized (alternating singles/doubles in the lactone-only segments):
  Σ ≈ 39.x — matches expected

### Why the fixpoint over-marks the lactones
A lactone 6-ring has structure `-O-C(=O)-Cb-Cb-Cb-Cb-` (Cb shared with the
adjacent benzene). Under `ring_pi_contribution`:
- lactone `O` → Some(2) (pyrrole-like lone-pair donor)
- carbonyl `C` with exocyclic `=O` → Some(0)
- 4 `Cb` shared with the (already-aromatic) benzene → Some(1) each = 4
- Σπ = 6 → satisfies Hückel (4n+2) → ring marked aromatic

But RDKit/mordred apply additional Kekulé constraints (or aromaticity model =
`rdkit` vs `mdl`) that reject this ring.

### Candidate fixes (pick one, validate with full parity run)

1. **Reject rings that contain an exocyclic-double-bonded C as a ring atom.**
   Simplest. **Verified to regress caffeine pyrimidinedione** (both carbonyl
   C's would trigger rejection). Don't use.
2. **Use RDKit's aromaticity model**: skip rings where any in-ring sp² atom's
   exocyclic-π hybridization is incompatible with the ring's π system
   (i.e., a `C=O` carbonyl out-of-plane). Requires modeling RDKit's
   Kekulé enforcement; non-trivial.
3. **Configurable aromaticity model**: expose `AromaticityModel::Rdkit` (the
   default upstream mordred uses) vs `AromaticityModel::Mdl`/`Hueckel` and run
   parity in Rdkit mode. Most work, cleanest result. **Recommended path
   forward — needs its own session.**

Geometric distinguishers between ellagic acid lactone (over-marked) and
caffeine pyrimidinedione (correctly marked) all triggered on caffeine too in
back-of-envelope checks. The real difference is RDKit's broader Kekulé
constraint, not a single-atom property.

## ~~Untracked WIP descriptor files~~ — landed 2026-06-14

All four files now wired into `mordred-core/src/descriptor/mod.rs` with
their unit tests green (125/0 in `cargo test -p mordred-core --lib`):

| File | LoC | Descriptor | Tests | Status |
|---|---:|---|---:|---|
| `balaban_j.rs` | 202 | BalabanJ topological index | 6 | landed |
| `path_count.rs` | 373 | MPC/TMPC/PIPC + nPath\* counts | 8 | landed (1-line API fix: `node_bound` → `node_count`) |
| `topo_psa.rs` | 387 | Topological polar surface area (Ertl) | 6 | landed |
| `vabc.rs` | 129 | VABC topological volume | 3 | landed |

Parity gate (Python pytest in `mordred-py`) not exercised in the loop session
that landed these — `make test-python` needs `maturin develop` + pytest setup.
**Run that before the next PR-ready push** to confirm none of the 7374 xfails
regress to fail with the new descriptor outputs.

## Blockers / watch-outs
- **gh auth dead**: `gh auth status` reports "token in default is invalid".
  Re-mint a token with `repo,workflow,read:org` scopes before any PR ops.
- **Push blocked** until gh re-auth: `git push` over HTTPS prompts for creds
  and aborts non-interactively.
- **Tracking PR #34**: open per memory but state unverified — check after re-auth.
- **Don't auto-push**. User pushes.

## Next actions (ordered)

1. Re-mint gh token, `git push` outstanding commits (`test/upstream-parity-gate`
   has new BalabanJ/TopoPSA/VABC/PathCount commits beyond `origin/HEAD`).
2. Run the parity gate (`cd mordred-py && maturin develop && pytest tests/`)
   to confirm:
   - Existing 624/1 baseline preserved.
   - The 4 newly-landed descriptors haven't promoted any of the 7374 xfails
     to unexpected fails.
3. **BIC0 / EllagicAcid**: Option 3 (configurable RDKit-aromaticity model)
   is the recommended path; deserves its own session. Don't attempt Options 1
   or 2 — see analysis above.
4. Update / close issues #31–#33 based on parity gate state.
5. Push, advance PR #34 to ready-for-review.
