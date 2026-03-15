use thiserror::Error;

/// Errors that can occur during descriptor calculation.
#[derive(Debug, Error)]
pub enum MordredError {
    #[error("SMILES parse error: {0}")]
    SmilesParseError(String),

    #[error("descriptor calculation error: {0}")]
    CalculationError(String),

    #[error("invalid molecule: {0}")]
    InvalidMolecule(String),
}
