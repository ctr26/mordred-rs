use pyo3::prelude::*;
use std::collections::HashMap;

use mordred_core::{DescriptorSet, parse_smiles};

/// Internal result type wrapping descriptor calculation results.
#[pyclass(name = "RustResult")]
struct RustResult {
    names: Vec<String>,
    values: Vec<Option<f64>>,
    name_to_index: HashMap<String, usize>,
}

#[pymethods]
impl RustResult {
    fn __len__(&self) -> usize {
        self.names.len()
    }

    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
        if let Ok(mut idx) = key.extract::<isize>() {
            let len = self.names.len() as isize;
            if idx < 0 {
                idx += len;
            }
            if idx < 0 || idx >= len {
                return Err(pyo3::exceptions::PyIndexError::new_err(
                    "index out of range",
                ));
            }
            Ok(self.values[idx as usize])
        } else if let Ok(name) = key.extract::<String>() {
            match self.name_to_index.get(&name) {
                Some(&idx) => Ok(self.values[idx]),
                None => Err(pyo3::exceptions::PyKeyError::new_err(name)),
            }
        } else {
            Err(pyo3::exceptions::PyTypeError::new_err(
                "index must be int or str",
            ))
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyResult<Py<RustResultIter>> {
        let iter = RustResultIter {
            values: slf.values.clone(),
            index: 0,
        };
        Py::new(slf.py(), iter)
    }

    fn __repr__(&self) -> String {
        let items: Vec<String> = self
            .names
            .iter()
            .zip(self.values.iter())
            .map(|(name, val)| match val {
                Some(v) => format!("{name}: {v}"),
                None => format!("{name}: None"),
            })
            .collect();
        format!("Result({{{}}})", items.join(", "))
    }

    /// Return results as a dict mapping descriptor name to value.
    fn to_dict(&self) -> HashMap<String, Option<f64>> {
        self.names
            .iter()
            .zip(self.values.iter())
            .map(|(name, val)| (name.clone(), *val))
            .collect()
    }

    /// Return descriptor names.
    fn keys(&self) -> Vec<String> {
        self.names.clone()
    }

    /// Return descriptor values in order.
    fn values(&self) -> Vec<Option<f64>> {
        self.values.clone()
    }
}

/// Iterator over RustResult values.
#[pyclass]
struct RustResultIter {
    values: Vec<Option<f64>>,
    index: usize,
}

#[pymethods]
impl RustResultIter {
    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }

    fn __next__(&mut self) -> Option<Option<f64>> {
        if self.index < self.values.len() {
            let val = self.values[self.index];
            self.index += 1;
            Some(val)
        } else {
            None
        }
    }
}

/// Build a RustResult from the core calculation output.
fn build_result(calc_output: Vec<(&str, Result<f64, mordred_core::MordredError>)>) -> RustResult {
    let mut names = Vec::with_capacity(calc_output.len());
    let mut values = Vec::with_capacity(calc_output.len());
    let mut name_to_index = HashMap::with_capacity(calc_output.len());

    for (i, (name, result)) in calc_output.into_iter().enumerate() {
        let name_string = name.to_string();
        name_to_index.insert(name_string.clone(), i);
        names.push(name_string);
        values.push(result.ok());
    }

    RustResult {
        names,
        values,
        name_to_index,
    }
}

/// Build a RustResult with all None values for a failed parse.
fn build_error_result(descriptor_set: &DescriptorSet) -> RustResult {
    let names: Vec<String> = descriptor_set
        .names()
        .into_iter()
        .map(|s| s.to_string())
        .collect();
    let values = vec![None; names.len()];
    let name_to_index: HashMap<String, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.clone(), i))
        .collect();

    RustResult {
        names,
        values,
        name_to_index,
    }
}

/// Internal calculator backed by Rust DescriptorSet.
#[pyclass(name = "RustCalculator")]
struct RustCalculator {
    descriptor_set: DescriptorSet,
}

#[pymethods]
impl RustCalculator {
    #[new]
    fn new() -> Self {
        Self {
            descriptor_set: DescriptorSet::all(),
        }
    }

    /// Calculate descriptors for a SMILES string, returns RustResult.
    fn calculate(&self, smiles: &str) -> PyResult<RustResult> {
        let mol = parse_smiles(smiles)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        let results = self.descriptor_set.calculate(&mol);
        Ok(build_result(results))
    }

    /// Calculate for a batch of SMILES, returns list of RustResult.
    fn calculate_batch(&self, smiles_list: Vec<String>) -> PyResult<Vec<RustResult>> {
        let results: Vec<RustResult> = smiles_list
            .iter()
            .map(|smi| match parse_smiles(smi) {
                Ok(mol) => build_result(self.descriptor_set.calculate(&mol)),
                Err(_) => build_error_result(&self.descriptor_set),
            })
            .collect();
        Ok(results)
    }

    /// List descriptor names.
    fn descriptor_names(&self) -> Vec<String> {
        self.descriptor_set
            .names()
            .into_iter()
            .map(|s| s.to_string())
            .collect()
    }

    /// Number of descriptors.
    fn __len__(&self) -> usize {
        self.descriptor_set.len()
    }
}

/// mordred._mordred_core native extension module.
#[pymodule]
fn _mordred_core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<RustCalculator>()?;
    m.add_class::<RustResult>()?;
    Ok(())
}
