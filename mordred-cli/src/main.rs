use clap::Parser;
use mordred_core::{DescriptorSet, parse_smiles};
use std::io::{self, BufRead};

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

    // Collect SMILES to process
    let smiles_list: Vec<String> = if let Some(ref smi) = cli.smiles {
        vec![smi.clone()]
    } else if let Some(ref path) = cli.file {
        let file = std::fs::File::open(path).unwrap_or_else(|e| {
            eprintln!("Error opening file '{}': {}", path, e);
            std::process::exit(1);
        });
        io::BufReader::new(file)
            .lines()
            .map_while(Result::ok)
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split_whitespace().next().unwrap_or("").to_string())
            .filter(|s| !s.is_empty())
            .collect()
    } else {
        eprintln!("Error: provide a SMILES string or --file. Use --help for usage.");
        std::process::exit(1);
    };

    let names: Vec<&str> = descriptor_set.names();

    match cli.format.as_str() {
        "csv" => {
            let mut writer: Box<dyn io::Write> = if let Some(ref path) = cli.output {
                Box::new(std::fs::File::create(path).unwrap_or_else(|e| {
                    eprintln!("Error creating output file: {}", e);
                    std::process::exit(1);
                }))
            } else {
                Box::new(io::stdout())
            };

            // Header
            write!(writer, "SMILES").unwrap();
            for name in &names {
                write!(writer, ",{}", name).unwrap();
            }
            writeln!(writer).unwrap();

            // Data
            for smi in &smiles_list {
                write!(writer, "{}", smi).unwrap();
                match parse_smiles(smi) {
                    Ok(mol) => {
                        let results = descriptor_set.calculate(&mol);
                        for (_, result) in &results {
                            match result {
                                Ok(val) => write!(writer, ",{:.6}", val).unwrap(),
                                Err(_) => write!(writer, ",").unwrap(),
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!("Warning: failed to parse '{}': {}", smi, e);
                        for _ in &names {
                            write!(writer, ",").unwrap();
                        }
                    }
                }
                writeln!(writer).unwrap();
            }
        }
        _ => {
            for smi in &smiles_list {
                match parse_smiles(smi) {
                    Ok(mol) => {
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
                    Err(e) => {
                        eprintln!("Error parsing '{}': {}", smi, e);
                    }
                }
            }
        }
    }
}
