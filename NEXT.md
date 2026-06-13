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
   Simplest. May regress caffeine pyrimidinedione (test at line 291) — check.
2. **Use RDKit's aromaticity model**: skip rings where any in-ring sp² atom's
   exocyclic-π hybridization is incompatible with the ring's π system
   (i.e., a `C=O` carbonyl out-of-plane).
3. **Configurable aromaticity model**: expose `AromaticityModel::Rdkit` (the
   default upstream mordred uses) vs `AromaticityModel::Mdl`/`Hueckel` and run
   parity in Rdkit mode. Most work, cleanest result.

Option 1 is the smallest patch. Risk: regression on caffeine
(`test_caffeine_pyrimidinedione_aromatic` at aromaticity.rs:291). Validate by
running the full `cargo test --workspace` AND the parity gate after the change.

## Untracked WIP descriptor files (decide before next push)

Four files staged in `mordred-core/src/descriptor/` (1090 LoC total):

| File | LoC | Descriptor | Status |
|---|---:|---|---|
| `balaban_j.rs` | 202 | BalabanJ topological index | needs unit tests + parity check |
| `path_count.rs` | 373 | nPath1–nPath10 (n-bond walk counts) | needs parity gate |
| `topo_psa.rs` | 387 | Topological polar surface area | needs Ertl-table validation |
| `vabc.rs` | 129 | VABC (van der Waals adjusted) | needs unit tests |

**Action options:**
- (a) Land each as separate commit with tests + parity entry — recommended for
  rolling them into PR #34's scope.
- (b) `git stash` or move to `wip/` branch if not ready for review.
- (c) Delete if superseded — confirm none of these match a more-recent
  implementation already on the branch.

Do NOT commit them with the aromaticity fix in the same commit.

## Blockers / watch-outs
- **gh auth dead**: `gh auth status` reports "token in default is invalid".
  Re-mint a token with `repo,workflow,read:org` scopes before any PR ops.
- **Push blocked** until gh re-auth: `git push` over HTTPS prompts for creds
  and aborts non-interactively.
- **Tracking PR #34**: open per memory but state unverified — check after re-auth.
- **Don't auto-push**. User pushes.

## Next actions (ordered)

1. Re-mint gh token, `git push` outstanding commits (none beyond
   `test/upstream-parity-gate` tip — branch is up-to-date with origin).
2. Triage the 4 untracked WIP descriptor files: land or stash.
3. Apply Option 1 fix above; run `cargo test --workspace`; if green, run the
   parity gate and confirm 625 pass / 0 fail; commit.
4. Update / close issues #31–#33 based on parity gate state.
5. Push, advance PR #34 to ready-for-review.
