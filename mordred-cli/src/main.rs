use clap::Parser;
use mordred_core::{DescriptorSet, parse_smiles};
use rayon::prelude::*;
use std::io::{self, BufRead, BufWriter, Write};

#[derive(Parser)]
#[command(name = "mordred", about = "Molecular descriptor calculator", version)]
struct Cli {
    /// SMILES string to calculate descriptors for
    smiles: Option<String>,

    /// Input file with one SMILES per line
    #[arg(short, long)]
    file: Option<String>,

    /// Output file (CSV format)
    #[arg(short, long)]
    output: Option<String>,

    /// List available descriptors
    #[arg(short, long)]
    list: bool,

    /// Output format: json or csv
    #[arg(long, default_value = "json")]
    format: String,

    /// Number of SMILES to process per chunk in batch mode
    #[arg(long, default_value = "1024")]
    chunk_size: usize,
}

fn main() {
    let cli = Cli::parse();
    let descriptor_set = DescriptorSet::all();

    if cli.list {
        println!("Available descriptors:");
        for (name, desc) in descriptor_set.list() {
            println!("  {:<20} {}", name, desc);
        }
        return;
    }

    // Single SMILES mode
    if let Some(ref smi) = cli.smiles {
        output_single(&cli, &descriptor_set, smi);
        return;
    }

    // File batch mode — chunked streaming
    if let Some(ref path) = cli.file {
        let file = std::fs::File::open(path).unwrap_or_else(|e| {
            eprintln!("Error opening file '{}': {}", path, e);
            std::process::exit(1);
        });
        let reader = io::BufReader::new(file);
        process_stream(&cli, &descriptor_set, reader);
        return;
    }

    eprintln!("Error: provide a SMILES string or --file. Use --help for usage.");
    std::process::exit(1);
}

fn output_single(cli: &Cli, descriptor_set: &DescriptorSet, smi: &str) {
    match parse_smiles(smi) {
        Ok(mol) => {
            let results = descriptor_set.calculate(&mol);
            match cli.format.as_str() {
                "csv" => {
                    let mut writer: BufWriter<Box<dyn Write>> = make_writer(cli);
                    let names = descriptor_set.names();
                    write!(writer, "SMILES").unwrap();
                    for name in &names {
                        write!(writer, ",{}", name).unwrap();
                    }
                    writeln!(writer).unwrap();
                    write!(writer, "{}", smi).unwrap();
                    for (_, result) in &results {
                        match result {
                            Ok(val) => write!(writer, ",{:.6}", val).unwrap(),
                            Err(_) => write!(writer, ",").unwrap(),
                        }
                    }
                    writeln!(writer).unwrap();
                }
                _ => {
                    let mut map = serde_json::Map::new();
                    map.insert(
                        "SMILES".to_string(),
                        serde_json::Value::String(smi.to_string()),
                    );
                    for (name, result) in &results {
                        match result {
                            Ok(val) => {
                                map.insert(
                                    name.to_string(),
                                    serde_json::Value::Number(
                                        serde_json::Number::from_f64(*val)
                                            .unwrap_or(serde_json::Number::from(0)),
                                    ),
                                );
                            }
                            Err(_) => {
                                map.insert(name.to_string(), serde_json::Value::Null);
                            }
                        }
                    }
                    let json = serde_json::Value::Object(map);
                    if let Some(ref path) = cli.output {
                        let file = std::fs::File::create(path).unwrap_or_else(|e| {
                            eprintln!("Error creating output file: {}", e);
                            std::process::exit(1);
                        });
                        serde_json::to_writer_pretty(file, &json).unwrap();
                    } else {
                        println!("{}", serde_json::to_string_pretty(&json).unwrap());
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("Error parsing '{}': {}", smi, e);
        }
    }
}

fn process_stream(cli: &Cli, descriptor_set: &DescriptorSet, reader: io::BufReader<std::fs::File>) {
    let names = descriptor_set.names();
    let chunk_size = cli.chunk_size.max(1);
    let mut lines = reader.lines();
    let mut chunk: Vec<String> = Vec::with_capacity(chunk_size);

    match cli.format.as_str() {
        "csv" => {
            let mut writer: BufWriter<Box<dyn Write>> = make_writer(cli);

            // Header
            write!(writer, "SMILES").unwrap();
            for name in &names {
                write!(writer, ",{}", name).unwrap();
            }
            writeln!(writer).unwrap();

            loop {
                chunk.clear();
                for _ in 0..chunk_size {
                    match lines.next() {
                        Some(Ok(line)) => {
                            let trimmed = line.trim().to_string();
                            if let Some(smi) = trimmed.split_whitespace().next() {
                                if !smi.is_empty() {
                                    chunk.push(smi.to_string());
                                }
                            }
                        }
                        Some(Err(_)) => continue,
                        None => break,
                    }
                }

                if chunk.is_empty() {
                    break;
                }

                let results: Vec<(String, Option<Vec<Option<f64>>>)> = chunk
                    .par_iter()
                    .map(|smi| match parse_smiles(smi) {
                        Ok(mol) => {
                            let vals: Vec<Option<f64>> = descriptor_set
                                .calculate(&mol)
                                .into_iter()
                                .map(|(_, r)| r.ok())
                                .collect();
                            (smi.clone(), Some(vals))
                        }
                        Err(e) => {
                            eprintln!("Warning: failed to parse '{}': {}", smi, e);
                            (smi.clone(), None)
                        }
                    })
                    .collect();

                for (smi, vals) in &results {
                    write!(writer, "{}", smi).unwrap();
                    match vals {
                        Some(vals) => {
                            for v in vals {
                                match v {
                                    Some(val) => write!(writer, ",{:.6}", val).unwrap(),
                                    None => write!(writer, ",").unwrap(),
                                }
                            }
                        }
                        None => {
                            for _ in &names {
                                write!(writer, ",").unwrap();
                            }
                        }
                    }
                    writeln!(writer).unwrap();
                }
            }
        }
        _ => loop {
            chunk.clear();
            for _ in 0..chunk_size {
                match lines.next() {
                    Some(Ok(line)) => {
                        let trimmed = line.trim().to_string();
                        if let Some(smi) = trimmed.split_whitespace().next() {
                            if !smi.is_empty() {
                                chunk.push(smi.to_string());
                            }
                        }
                    }
                    Some(Err(_)) => continue,
                    None => break,
                }
            }

            if chunk.is_empty() {
                break;
            }

            let results: Vec<_> = chunk
                .par_iter()
                .filter_map(|smi| {
                    parse_smiles(smi).ok().map(|mol| {
                        let results = descriptor_set.calculate(&mol);
                        let mut map = serde_json::Map::new();
                        map.insert("SMILES".to_string(), serde_json::Value::String(smi.clone()));
                        for (name, result) in &results {
                            match result {
                                Ok(val) => {
                                    map.insert(
                                        name.to_string(),
                                        serde_json::Value::Number(
                                            serde_json::Number::from_f64(*val)
                                                .unwrap_or(serde_json::Number::from(0)),
                                        ),
                                    );
                                }
                                Err(_) => {
                                    map.insert(name.to_string(), serde_json::Value::Null);
                                }
                            }
                        }
                        serde_json::Value::Object(map)
                    })
                })
                .collect();

            for json in &results {
                println!("{}", serde_json::to_string_pretty(json).unwrap());
            }
        },
    }
}

fn make_writer(cli: &Cli) -> BufWriter<Box<dyn Write>> {
    BufWriter::new(if let Some(ref path) = cli.output {
        Box::new(std::fs::File::create(path).unwrap_or_else(|e| {
            eprintln!("Error creating output file: {}", e);
            std::process::exit(1);
        }))
    } else {
        Box::new(io::stdout())
    })
}
