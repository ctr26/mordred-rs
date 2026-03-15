"""Tests for the mordred Calculator API."""

import pytest


def test_import():
    """Calculator and descriptors are importable from mordred."""
    from mordred import Calculator, descriptors
    assert Calculator is not None
    assert descriptors is not None


def test_calculator_creation():
    """Calculator can be created with descriptors module."""
    from mordred import Calculator, descriptors
    calc = Calculator(descriptors, ignore_3D=True)
    assert len(calc) > 0


def test_calculator_default():
    """Calculator with no args creates with all descriptors."""
    from mordred import Calculator
    calc = Calculator()
    assert len(calc) > 0


def test_calculate_smiles():
    """Calculate descriptors from a SMILES string."""
    from mordred import Calculator, descriptors
    calc = Calculator(descriptors)
    result = calc("CCO")
    assert len(result) > 0


def test_result_index_access():
    """Result supports integer indexing."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    val = result[0]
    assert isinstance(val, (float, type(None)))


def test_result_name_access():
    """Result supports string key access."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    val = result["nAtom"]
    assert isinstance(val, float)


def test_result_negative_index():
    """Result supports negative indexing."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    val = result[-1]
    assert isinstance(val, (float, type(None)))


def test_result_len():
    """Result has correct length."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    assert len(result) == len(calc)


def test_result_iter():
    """Result is iterable."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    values = list(result)
    assert len(values) == len(calc)


def test_result_to_dict():
    """Result can be converted to dict."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    d = result.to_dict()
    assert isinstance(d, dict)
    assert "nAtom" in d


def test_descriptors_property():
    """Calculator.descriptors returns tuple of names."""
    from mordred import Calculator
    calc = Calculator()
    descs = calc.descriptors
    assert isinstance(descs, tuple)
    assert "MW" in descs
    assert "nAtom" in descs


def test_known_values_ethanol():
    """Ethanol descriptors match known values."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("CCO")
    # CH3CH2OH: 3 heavy atoms, 9 total atoms, MW ≈ 46.07
    assert result["nHeavyAtom"] == 3.0
    assert result["nAtom"] == 9.0
    assert abs(result["MW"] - 46.069) < 0.1


def test_known_values_benzene():
    """Benzene descriptors match known values."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("c1ccccc1")
    assert result["nHeavyAtom"] == 6.0
    assert result["nC"] == 6.0
    assert abs(result["MW"] - 78.114) < 0.1
    assert result["nRing"] == 1.0


def test_invalid_smiles():
    """Invalid SMILES raises ValueError."""
    from mordred import Calculator
    calc = Calculator()
    with pytest.raises(ValueError):
        calc("not_a_smiles_XYZ123")


def test_pandas():
    """Calculator.pandas returns DataFrame."""
    pandas = pytest.importorskip("pandas")
    from mordred import Calculator
    calc = Calculator()
    df = calc.pandas(["CCO", "c1ccccc1"])
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == 2
    assert "MW" in df.columns


def test_map():
    """Calculator.map yields results."""
    from mordred import Calculator
    calc = Calculator()
    results = list(calc.map(["CCO", "C"]))
    assert len(results) == 2


def test_str_input():
    """Calculator accepts plain string SMILES."""
    from mordred import Calculator
    calc = Calculator()
    result = calc("C")
    assert result["nAtom"] == 5.0  # CH4
