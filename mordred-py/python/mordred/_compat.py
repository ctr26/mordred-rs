"""Compatibility utilities for RDKit integration."""


def mol_to_smiles(mol):
    """Convert an RDKit Mol object to a SMILES string.

    Args:
        mol: An RDKit ``Chem.Mol`` instance.

    Returns:
        A canonical SMILES string.

    Raises:
        ImportError: If RDKit is not installed.
    """
    try:
        from rdkit import Chem

        return Chem.MolToSmiles(mol)
    except ImportError:
        raise ImportError(
            "RDKit is required to convert Mol objects. "
            "Install it or pass SMILES strings instead."
        )
