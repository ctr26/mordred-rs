use std::time::Instant;

use mordred_core::{DescriptorSet, parse_smiles};

fn main() {
    let molecules: Vec<(&str, &str)> = vec![
        // Simple molecules
        ("methane", "C"),
        ("ethanol", "CCO"),
        ("acetic_acid", "CC(=O)O"),
        ("water", "O"),
        // Aromatic
        ("benzene", "c1ccccc1"),
        ("naphthalene", "c1ccc2ccccc2c1"),
        ("phenol", "c1ccc(cc1)O"),
        ("aniline", "c1ccc(cc1)N"),
        ("toluene", "Cc1ccccc1"),
        // Drug-like
        ("caffeine", "Cn1cnc2c1c(=O)n(c(=O)n2C)C"),
        ("aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
        ("ibuprofen", "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"),
        ("paracetamol", "CC(=O)Nc1ccc(O)cc1"),
        ("nicotine", "CN1CCC[C@H]1c1cccnc1"),
        // Amino acids
        ("glycine", "NCC(=O)O"),
        ("alanine", "CC(N)C(=O)O"),
        ("phenylalanine", "NC(Cc1ccccc1)C(=O)O"),
        ("tryptophan", "NC(Cc1c[nH]c2ccccc12)C(=O)O"),
        // Halogenated
        ("chloroform", "ClC(Cl)Cl"),
        ("mixed_halogen", "FC(Cl)(Br)I"),
        ("fluorobenzene", "Fc1ccccc1"),
        // Heterocyclic
        ("pyridine", "c1ccncc1"),
        ("furan", "c1ccoc1"),
        ("thiophene", "c1ccsc1"),
        ("imidazole", "c1cnc[nH]1"),
        // Larger molecules
        ("citric_acid", "OC(=O)CC(O)(CC(=O)O)C(=O)O"),
        ("glucose", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O"),
        (
            "cholesterol",
            "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
        ),
        ("testosterone", "CC12CCC3C(C1CCC2O)CCC4=CC(=O)CCC34C"),
        // Polymeric fragments
        ("tripeptide_gly", "NCC(=O)NCC(=O)NCC(=O)O"),
        ("diethyl_ether", "CCOCC"),
        ("dmso", "CS(=O)C"),
        // Fused rings
        ("indole", "c1ccc2[nH]ccc2c1"),
        ("quinoline", "c1ccc2ncccc2c1"),
        ("anthracene", "c1ccc2cc3ccccc3cc2c1"),
        // Strained / unusual
        ("cyclopropane", "C1CC1"),
        ("cubane", "C12C3C4C1C5C3C4C25"),
        ("adamantane", "C1C2CC3CC1CC(C2)C3"),
        // Larger drug
        ("metformin", "CN(C)C(=N)NC(=N)N"),
        (
            "sildenafil",
            "CCCc1nn(C)c2c1nc(nc2OCC)c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1",
        ),
    ];

    let set = DescriptorSet::all();

    // Warmup
    for (_, smi) in &molecules {
        if let Ok(mol) = parse_smiles(smi) {
            set.calculate(&mol);
        }
    }

    // Pre-parse all molecules
    let parsed: Vec<_> = molecules
        .iter()
        .filter_map(|(name, smi)| parse_smiles(smi).ok().map(|mol| (*name, *smi, mol)))
        .collect();

    println!(
        "mordred-rs benchmark ({} descriptors, {} molecules)",
        set.len(),
        parsed.len()
    );
    println!("=====================================================");

    // Per-molecule timings (calc only)
    println!("\n--- Per-molecule (calc only, 10000 iterations) ---");
    for (name, _smi, mol) in &parsed {
        let n = 10_000;
        let start = Instant::now();
        for _ in 0..n {
            let _ = set.calculate(mol);
        }
        let elapsed = start.elapsed();
        let us = elapsed.as_nanos() as f64 / n as f64 / 1000.0;
        println!("{:20}: {:8.2} us/mol", name, us);
    }

    // Aggregate: calc only
    let n = 5_000;
    let start = Instant::now();
    for _ in 0..n {
        for (_, _, mol) in &parsed {
            let _ = set.calculate(mol);
        }
    }
    let elapsed = start.elapsed();
    let total_mols = n * parsed.len();
    println!(
        "\n--- Aggregate (calc only) ---\nTotal: {} molecules in {:.2} ms ({:.2} us/mol)",
        total_mols,
        elapsed.as_millis(),
        elapsed.as_nanos() as f64 / total_mols as f64 / 1000.0
    );

    // Aggregate: parse + calc
    let smiles_list: Vec<&str> = molecules.iter().map(|(_, s)| *s).collect();
    let n = 2_000;
    let start = Instant::now();
    for _ in 0..n {
        for s in &smiles_list {
            if let Ok(mol) = parse_smiles(s) {
                let _ = set.calculate(&mol);
            }
        }
    }
    let elapsed = start.elapsed();
    let total_mols = n * smiles_list.len();
    println!(
        "\n--- Aggregate (parse + calc) ---\nTotal: {} molecules in {:.2} ms ({:.2} us/mol)",
        total_mols,
        elapsed.as_millis(),
        elapsed.as_nanos() as f64 / total_mols as f64 / 1000.0
    );

    // Throughput
    let duration_secs = elapsed.as_secs_f64();
    println!(
        "Throughput: {:.0} molecules/sec",
        total_mols as f64 / duration_secs
    );
}
