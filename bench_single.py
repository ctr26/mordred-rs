"""Benchmark a single Python mordred installation."""

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

N_ITER = 10


def run_benchmark():
    from rdkit import Chem
    from mordred import Calculator, descriptors

    # Detect version
    try:
        from mordred import __version__ as mv
        version = f"mordred v{mv}"
    except ImportError:
        version = "mordred (unknown version)"
    # Detect if it's community fork
    try:
        from importlib.metadata import version as pkg_version
        cv = pkg_version("mordredcommunity")
        version = f"mordred-community v{cv}"
    except Exception:
        try:
            cv = pkg_version("mordred")
            version = f"mordred (original) v{cv}"
        except Exception:
            pass

    calc = Calculator(descriptors, ignore_3D=True)
    n_desc = len(calc)

    mols = []
    for name, smi in MOLECULES:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mols.append((name, mol))

    # Warmup
    for _, mol in mols[:3]:
        calc(mol)

    # Per-molecule
    per_mol = {}
    for name, mol in mols:
        start = time.perf_counter()
        for _ in range(N_ITER):
            calc(mol)
        elapsed = time.perf_counter() - start
        per_mol[name] = elapsed / N_ITER * 1_000_000  # us

    # Aggregate
    start = time.perf_counter()
    for _ in range(N_ITER):
        for _, mol in mols:
            calc(mol)
    elapsed = time.perf_counter() - start
    total_mols = N_ITER * len(mols)
    avg_us = elapsed / total_mols * 1_000_000

    result = {
        "name": version,
        "descriptors": n_desc,
        "molecules": len(mols),
        "per_molecule_us": per_mol,
        "aggregate_avg_us": avg_us,
        "aggregate_total_s": elapsed,
        "total_molecules": total_mols,
        "throughput_mol_per_sec": total_mols / elapsed,
    }

    print(json.dumps(result))


if __name__ == "__main__":
    run_benchmark()
