"""Benchmark Python mordred (original) and mordred-community against mordred-rs."""

import json
import sys
import time

MOLECULES = [
    ("methane", "C"),
    ("ethanol", "CCO"),
    ("acetic_acid", "CC(=O)O"),
    ("water", "O"),
    ("benzene", "c1ccccc1"),
    ("naphthalene", "c1ccc2ccccc2c1"),
    ("phenol", "c1ccc(cc1)O"),
    ("aniline", "c1ccc(cc1)N"),
    ("toluene", "Cc1ccccc1"),
    ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
    ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
    ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
    ("glycine", "NCC(=O)O"),
    ("alanine", "CC(N)C(=O)O"),
    ("phenylalanine", "NC(Cc1ccccc1)C(=O)O"),
    ("tryptophan", "NC(Cc1c[nH]c2ccccc12)C(=O)O"),
    ("chloroform", "ClC(Cl)Cl"),
    ("mixed_halogen", "FC(Cl)(Br)I"),
    ("fluorobenzene", "Fc1ccccc1"),
    ("pyridine", "c1ccncc1"),
    ("furan", "c1ccoc1"),
    ("thiophene", "c1ccsc1"),
    ("imidazole", "c1cnc[nH]1"),
    ("citric_acid", "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
    ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
    ("cholesterol", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"),
    ("testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
    ("tripeptide_gly", "NCC(=O)NCC(=O)NCC(=O)O"),
    ("diethyl_ether", "CCOCC"),
    ("dmso", "CS(=O)C"),
    ("indole", "c1ccc2[nH]ccc2c1"),
    ("quinoline", "c1ccc2ncccc2c1"),
    ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
    ("cyclopropane", "C1CC1"),
    ("cubane", "C12C3C4C1C5C3C4C25"),
    ("adamantane", "C1C2CC3CC1CC(C2)C3"),
    ("metformin", "CN(C)C(=N)NC(=N)N"),
    ("sildenafil", "CCCc1nn(C)c2c1nc(nc2OCC)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"),
]

N_ITER = 10  # iterations over all molecules


def bench_mordred_original():
    """Benchmark original mordred (1.2.0)."""
    from rdkit import Chem
    from mordred import Calculator, descriptors

    calc = Calculator(descriptors, ignore_3D=True)
    n_desc = len(calc)

    # Pre-parse with RDKit
    mols = [(name, Chem.MolFromSmiles(smi)) for name, smi in MOLECULES]
    mols = [(name, mol) for name, mol in mols if mol is not None]

    # Warmup
    for _, mol in mols[:3]:
        calc(mol)

    # Per-molecule timing
    per_mol = {}
    for name, mol in mols:
        start = time.perf_counter()
        for _ in range(N_ITER):
            calc(mol)
        elapsed = time.perf_counter() - start
        per_mol[name] = elapsed / N_ITER * 1000  # ms

    # Aggregate
    start = time.perf_counter()
    for _ in range(N_ITER):
        for _, mol in mols:
            calc(mol)
    elapsed = time.perf_counter() - start
    total_mols = N_ITER * len(mols)
    avg_ms = elapsed / total_mols * 1000

    return {
        "name": f"mordred (Python) v1.2.0",
        "descriptors": n_desc,
        "molecules": len(mols),
        "per_molecule_ms": per_mol,
        "aggregate_avg_ms": avg_ms,
        "aggregate_total_s": elapsed,
        "total_molecules": total_mols,
        "throughput_mol_per_sec": total_mols / elapsed,
    }


def bench_mordred_community():
    """Benchmark mordred-community."""
    from rdkit import Chem

    # mordredcommunity uses same import path but different internals
    import mordredcommunity
    from mordredcommunity import Calculator, descriptors

    calc = Calculator(descriptors, ignore_3D=True)
    n_desc = len(calc)

    mols = [(name, Chem.MolFromSmiles(smi)) for name, smi in MOLECULES]
    mols = [(name, mol) for name, mol in mols if mol is not None]

    # Warmup
    for _, mol in mols[:3]:
        calc(mol)

    per_mol = {}
    for name, mol in mols:
        start = time.perf_counter()
        for _ in range(N_ITER):
            calc(mol)
        elapsed = time.perf_counter() - start
        per_mol[name] = elapsed / N_ITER * 1000

    start = time.perf_counter()
    for _ in range(N_ITER):
        for _, mol in mols:
            calc(mol)
    elapsed = time.perf_counter() - start
    total_mols = N_ITER * len(mols)
    avg_ms = elapsed / total_mols * 1000

    return {
        "name": f"mordred-community v{mordredcommunity.__version__}",
        "descriptors": n_desc,
        "molecules": len(mols),
        "per_molecule_ms": per_mol,
        "aggregate_avg_ms": avg_ms,
        "aggregate_total_s": elapsed,
        "total_molecules": total_mols,
        "throughput_mol_per_sec": total_mols / elapsed,
    }


def bench_mordred_community_2d_only():
    """Benchmark mordred-community with only 2D descriptors (fairer comparison)."""
    from rdkit import Chem

    import mordredcommunity
    from mordredcommunity import Calculator, descriptors

    calc = Calculator(descriptors, ignore_3D=True)
    n_desc = len(calc)

    mols = [(name, Chem.MolFromSmiles(smi)) for name, smi in MOLECULES]
    mols = [(name, mol) for name, mol in mols if mol is not None]

    # Warmup
    for _, mol in mols[:3]:
        calc(mol)

    per_mol = {}
    for name, mol in mols:
        start = time.perf_counter()
        for _ in range(N_ITER):
            calc(mol)
        elapsed = time.perf_counter() - start
        per_mol[name] = elapsed / N_ITER * 1000

    start = time.perf_counter()
    for _ in range(N_ITER):
        for _, mol in mols:
            calc(mol)
    elapsed = time.perf_counter() - start
    total_mols = N_ITER * len(mols)
    avg_ms = elapsed / total_mols * 1000

    return {
        "name": f"mordred-community v{mordredcommunity.__version__} (ignore_3D=True)",
        "descriptors": n_desc,
        "molecules": len(mols),
        "per_molecule_ms": per_mol,
        "aggregate_avg_ms": avg_ms,
        "aggregate_total_s": elapsed,
        "total_molecules": total_mols,
        "throughput_mol_per_sec": total_mols / elapsed,
    }


if __name__ == "__main__":
    results = []

    print("=" * 70)
    print("BENCHMARKING PYTHON MORDRED VARIANTS")
    print(f"Molecules: {len(MOLECULES)}, Iterations: {N_ITER}")
    print("=" * 70)

    # Original mordred
    print("\n>>> Benchmarking mordred (original) v1.2.0 ...")
    try:
        r = bench_mordred_original()
        results.append(r)
        print(f"    {r['descriptors']} descriptors, {r['aggregate_avg_ms']:.2f} ms/mol, "
              f"{r['throughput_mol_per_sec']:.0f} mol/s")
    except Exception as e:
        print(f"    FAILED: {e}")

    # mordred-community
    print("\n>>> Benchmarking mordred-community ...")
    try:
        r = bench_mordred_community()
        results.append(r)
        print(f"    {r['descriptors']} descriptors, {r['aggregate_avg_ms']:.2f} ms/mol, "
              f"{r['throughput_mol_per_sec']:.0f} mol/s")
    except Exception as e:
        print(f"    FAILED: {e}")

    # Print per-molecule table
    print("\n" + "=" * 70)
    print("PER-MOLECULE TIMES (ms)")
    print("=" * 70)

    header = f"{'Molecule':20}"
    for r in results:
        header += f" | {r['name']:>30}"
    print(header)
    print("-" * len(header))

    for name, _ in MOLECULES:
        row = f"{name:20}"
        for r in results:
            t = r["per_molecule_ms"].get(name, float("nan"))
            row += f" | {t:>30.2f}"
        print(row)

    # Save JSON
    with open("bench_python_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to bench_python_results.json")
