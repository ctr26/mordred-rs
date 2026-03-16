"""Compare mordred-rs descriptor values against Python mordred for shared descriptors."""

import json
import subprocess
import sys
import os

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

# Mapping from mordred-rs names to Python mordred descriptor names
# Python mordred uses different class names but same output names for most
RS_TO_PY = {
    # Constitutional
    "nAtom": "nAtom",
    "nHeavyAtom": "nHeavyAtom",
    "nBonds": "nBonds",
    "MW": "MW",
    # Topological
    "WPath": "WPath",
    "Zagreb1": "Zagreb1",
    "Zagreb2": "Zagreb2",
    # Atom counts
    "nC": "nC",
    "nN": "nN",
    "nO": "nO",
    "nS": "nS",
    "nP": "nP",
    "nX": "nX",
    "nH": "nH",
    "nB": "nB",
    "nF": "nF",
    "nCl": "nCl",
    "nBr": "nBr",
    "nI": "nI",
    "nHetero": "nHet",
    # Bond counts
    "nBondsS": "nBondsS",
    "nBondsD": "nBondsD",
    "nBondsT": "nBondsT",
    "nBondsA": "nBondsA",
    # Ring counts
    "nRing": "nRing",
    "n3Ring": "n3Ring",
    "n4Ring": "n4Ring",
    "n5Ring": "n5Ring",
    "n6Ring": "n6Ring",
    "n7Ring": "n7Ring",
    "n8Ring": "n8Ring",
    "n9Ring": "n9Ring",
    "n10Ring": "n10Ring",
    "n11Ring": "n11Ring",
    "n12Ring": "n12Ring",
    # Connectivity
    "Chi0": "Chi0",
    "Chi1": "Chi1",
}


def get_rust_values():
    """Get descriptor values from mordred-rs CLI."""
    results = {}
    for name, smi in MOLECULES:
        out = subprocess.run(
            ["target/release/mordred", smi],
            capture_output=True, text=True
        )
        if out.returncode == 0:
            data = json.loads(out.stdout)
            data.pop("SMILES", None)
            results[name] = data
    return results


def get_python_values():
    """Get descriptor values from Python mordred."""
    from rdkit import Chem
    from mordred import Calculator, descriptors

    calc = Calculator(descriptors, ignore_3D=True)
    desc_names = [str(d) for d in calc.descriptors]

    results = {}
    for name, smi in MOLECULES:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        r = calc(mol)
        d = {}
        for i, desc_name in enumerate(desc_names):
            val = r[i]
            if val is not None and not isinstance(val, Exception):
                try:
                    d[desc_name] = float(val)
                except (TypeError, ValueError):
                    pass
        results[name] = d
    return results


def compare():
    rust = get_rust_values()
    python = get_python_values()

    print("=" * 110)
    print("CORRECTNESS COMPARISON: mordred-rs vs Python mordred (like-for-like descriptors)")
    print("=" * 110)

    total_checks = 0
    matches = 0
    mismatches = []

    for mol_name, _ in MOLECULES:
        rs_vals = rust.get(mol_name, {})
        py_vals = python.get(mol_name, {})

        for rs_name, py_name in RS_TO_PY.items():
            rs_val = rs_vals.get(rs_name)
            py_val = py_vals.get(py_name)

            if rs_val is None or py_val is None:
                continue

            total_checks += 1
            # Use relative tolerance for MW, absolute for counts
            if rs_name == "MW":
                ok = abs(rs_val - py_val) < 0.01
            elif rs_name in ("Chi0", "Chi1"):
                ok = abs(rs_val - py_val) < 0.001
            else:
                ok = abs(rs_val - py_val) < 0.01

            if ok:
                matches += 1
            else:
                mismatches.append((mol_name, rs_name, py_name, rs_val, py_val))

    print(f"\nTotal comparisons: {total_checks}")
    print(f"Matches:           {matches}")
    print(f"Mismatches:        {len(mismatches)}")
    print(f"Match rate:        {matches/total_checks*100:.1f}%")

    if mismatches:
        print(f"\n{'Molecule':20} | {'rs_name':12} | {'py_name':12} | {'rs_value':>12} | {'py_value':>12} | {'delta':>10}")
        print("-" * 90)
        for mol_name, rs_name, py_name, rs_val, py_val in mismatches:
            delta = rs_val - py_val
            print(f"{mol_name:20} | {rs_name:12} | {py_name:12} | {rs_val:>12.4f} | {py_val:>12.4f} | {delta:>+10.4f}")

    # Print sample values for a few molecules
    print("\n" + "=" * 110)
    print("SAMPLE VALUES (benzene, caffeine, cholesterol)")
    print("=" * 110)
    for mol_name in ["benzene", "caffeine", "cholesterol"]:
        rs_vals = rust.get(mol_name, {})
        py_vals = python.get(mol_name, {})
        print(f"\n--- {mol_name} ---")
        print(f"{'Descriptor':12} | {'mordred-rs':>14} | {'Python mordred':>14} | {'Match':>6}")
        print("-" * 55)
        for rs_name, py_name in RS_TO_PY.items():
            rs_val = rs_vals.get(rs_name)
            py_val = py_vals.get(py_name)
            if rs_val is None and py_val is None:
                continue
            rs_str = f"{rs_val:.4f}" if rs_val is not None else "N/A"
            py_str = f"{py_val:.4f}" if py_val is not None else "N/A"
            if rs_val is not None and py_val is not None:
                ok = "OK" if abs(rs_val - py_val) < 0.01 else "DIFF"
            else:
                ok = "MISS"
            print(f"{rs_name:12} | {rs_str:>14} | {py_str:>14} | {ok:>6}")


if __name__ == "__main__":
    compare()
