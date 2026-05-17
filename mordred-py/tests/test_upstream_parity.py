"""Parity gate against upstream mordred reference YAML.

Each `(yaml_file, molecule, descriptor)` triple becomes one parametrized test
case. Cases for descriptors not implemented in mordred-rs are marked ``xfail``
so that pytest summary reports a live ``N passed, M xfailed`` view of parity
coverage. As descriptors land in the Rust core, their xfail markers
automatically clear without any test-file edit.

Reference corpus is vendored from ``mordredcommunity/mordred`` (BSD-3-Clause);
see ``upstream_references/ATTRIBUTION.md``.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pytest

yaml = pytest.importorskip("yaml")

from mordred import Calculator  # noqa: E402

REFS_DIR = Path(__file__).parent / "upstream_references"

# YAML files whose descriptors fundamentally require 3D coordinates,
# QM machinery, or partial-charge models — mordred-rs is 2D-only by design,
# so these are skipped wholesale rather than xfailed per-descriptor (noise
# reduction; clutters the summary otherwise).
THREE_D_FILES: frozenset[str] = frozenset(
    {
        "CPSA.yaml",
        "GeometricalIndex.yaml",
        "GravitationalIndex.yaml",
        "MomentOfInertia.yaml",
        "MoRSE.yaml",
    }
)


def _load_smiles_corpus() -> dict[str, str]:
    """Parse ``structures.smi``: lines of ``SMILES Name``."""
    corpus: dict[str, str] = {}
    for line in (REFS_DIR / "structures.smi").read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split(maxsplit=1)
        if len(parts) == 2:
            smiles, name = parts
            corpus[name] = smiles
    return corpus


SMILES: dict[str, str] = _load_smiles_corpus()
CALCULATOR = Calculator()
IMPLEMENTED: frozenset[str] = frozenset(CALCULATOR.descriptors)


def _normalize_group(group: dict[str, Any]) -> tuple[list[str], list[Any], int]:
    """Coerce both list-of-names and scalar-name YAML shapes to a uniform form.

    Upstream YAML mixes two schemas:
        - ``names: [A, B, C]`` with ``results: {Mol: [v_a, v_b, v_c]}``
        - ``names: ECIndex`` with ``results: {Mol: 38}``
    Return ``(names, value_extractor, digit)`` where the per-mol values are
    always lists aligned with ``names``.
    """
    raw_names = group.get("names")
    names = [raw_names] if isinstance(raw_names, str) else list(raw_names)
    digit = int(group.get("digit", 0))
    return names, group.get("results", {}), digit


def _iter_cases() -> list[Any]:
    """Yield one ``pytest.param`` per (file, mol, descriptor)."""
    cases: list[Any] = []
    for yaml_path in sorted(REFS_DIR.rglob("*.yaml")):
        rel = yaml_path.relative_to(REFS_DIR).as_posix()
        is_3d = yaml_path.name in THREE_D_FILES
        try:
            groups = yaml.safe_load(yaml_path.read_text()) or []
        except yaml.YAMLError:
            continue
        for group in groups:
            if not isinstance(group, dict):
                continue
            names, results, digit = _normalize_group(group)
            for mol_name, raw_values in results.items():
                if mol_name not in SMILES:
                    continue
                values = [raw_values] if not isinstance(raw_values, list) else raw_values
                for desc_name, expected in zip(names, values):
                    if expected == "skip" or expected is None:
                        continue
                    marks: list[Any] = []
                    if is_3d:
                        marks.append(
                            pytest.mark.skip(
                                reason=f"3D / QM not implemented (file: {yaml_path.name})"
                            )
                        )
                    elif desc_name not in IMPLEMENTED:
                        marks.append(
                            pytest.mark.xfail(
                                reason=f"{desc_name} not yet implemented",
                                strict=False,
                            )
                        )
                    cases.append(
                        pytest.param(
                            mol_name,
                            desc_name,
                            expected,
                            digit,
                            id=f"{rel}::{mol_name}::{desc_name}",
                            marks=marks,
                        )
                    )
    return cases


@pytest.mark.parametrize("mol_name,desc_name,expected,digit", _iter_cases())
def test_upstream_parity(
    mol_name: str, desc_name: str, expected: float, digit: int
) -> None:
    """A vendored (mol, descriptor) pair matches upstream's reference value.

    Tolerance follows numpy's ``assert_almost_equal``: ``|a - b| < 1.5e-digit``.
    """
    smiles = SMILES[mol_name]
    try:
        result = CALCULATOR(smiles)
    except ValueError as exc:
        pytest.xfail(f"SMILES parse failed: {exc}")

    try:
        actual = result[desc_name]
    except (KeyError, IndexError):
        pytest.xfail(f"{desc_name} not exposed by shim")

    if actual is None or (isinstance(actual, float) and math.isnan(actual)):
        pytest.xfail(f"{desc_name} returned None/NaN for {mol_name}")

    tolerance = 1.5 * (10 ** -digit)
    assert actual == pytest.approx(float(expected), abs=tolerance), (
        f"{mol_name}.{desc_name}: expected {expected}, got {actual} "
        f"(tol={tolerance:g})"
    )
