"""Mordred molecular descriptor calculator -- Rust-powered drop-in replacement."""

from mordred._mordred_core import RustCalculator, RustResult
from mordred import descriptors


class Calculator:
    """Molecular descriptor calculator.

    Args:
        descs: Descriptor specification. Pass the ``descriptors`` module for
            all descriptors, or ``None`` for an empty calculator.
        ignore_3D: If ``True``, skip 3D descriptors. Currently always ``True``
            as only 2D descriptors are implemented.
    """

    def __init__(self, descs=None, ignore_3D=False):
        self._inner = RustCalculator()
        self._ignore_3D = ignore_3D
        # *descs* is accepted for API compat; currently all descriptors are
        # always loaded since we only have a fixed set.

    def __call__(self, mol):
        """Calculate descriptors for a molecule.

        Args:
            mol: SMILES string or RDKit Mol object.

        Returns:
            Result object with descriptor values.
        """
        smiles = _to_smiles(mol)
        return self._inner.calculate(smiles)

    def __len__(self):
        return len(self._inner)

    @property
    def descriptors(self):
        """Tuple of descriptor names."""
        return tuple(self._inner.descriptor_names())

    def pandas(self, mols, quiet=False, ipynb=False, nproc=1, id=-1):
        """Calculate descriptors for multiple molecules.

        Args:
            mols: Iterable of SMILES strings or RDKit Mol objects.
            quiet: Suppress progress output.
            ipynb: Use Jupyter-friendly progress bars.
            nproc: Number of processes (currently ignored).
            id: Identifier for progress bars.

        Returns:
            ``pandas.DataFrame`` with descriptor values.
        """
        import pandas as pd

        smiles_list = [_to_smiles(m) for m in mols]
        results = self._inner.calculate_batch(smiles_list)
        data = [r.to_dict() for r in results]
        return pd.DataFrame(data, columns=list(self._inner.descriptor_names()))

    def map(self, mols, nproc=1, nmols=None, quiet=False, ipynb=False, id=-1):
        """Calculate descriptors over molecules, yielding results.

        Args:
            mols: Iterable of SMILES strings or RDKit Mol objects.
            nproc: Number of processes (currently ignored).
            nmols: Expected number of molecules (for progress bars).
            quiet: Suppress progress output.
            ipynb: Use Jupyter-friendly progress bars.
            id: Identifier for progress bars.

        Yields:
            Result objects for each molecule.
        """
        for mol in mols:
            yield self(mol)


def _to_smiles(mol):
    """Convert a molecule input to a SMILES string.

    Args:
        mol: SMILES string or RDKit Mol object.

    Returns:
        A SMILES string.

    Raises:
        TypeError: If *mol* is not a string or RDKit Mol.
        ImportError: If an RDKit Mol is passed but RDKit is not installed.
    """
    if isinstance(mol, str):
        return mol
    if hasattr(mol, "GetSmiles"):
        try:
            from rdkit import Chem

            return Chem.MolToSmiles(mol)
        except ImportError:
            raise ImportError(
                "RDKit is required to pass Mol objects. "
                "Install it or pass SMILES strings instead."
            )
    raise TypeError(
        f"Expected SMILES string or RDKit Mol, got {type(mol).__name__}"
    )
