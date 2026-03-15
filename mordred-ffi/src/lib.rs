use std::ffi::{CStr, CString};
use std::os::raw::c_char;
use std::ptr;

use mordred_core::{DescriptorSet, Molecule, parse_smiles};

// ---------------------------------------------------------------------------
// Opaque types exposed through the C API
// ---------------------------------------------------------------------------

/// Opaque handle to a descriptor calculator.
pub struct MordredCalculator {
    inner: DescriptorSet,
    /// Cached descriptor names as CStrings so that pointers returned by
    /// `mordred_calculator_descriptor_name` remain valid for the lifetime of
    /// the calculator.
    cached_names: Vec<CString>,
}

/// Opaque handle to a parsed molecule.
pub struct MordredMolecule {
    inner: Molecule,
}

/// Opaque handle to a calculation result set.
pub struct MordredResult {
    entries: Vec<MordredResultEntry>,
}

struct MordredResultEntry {
    name: CString,
    value: Option<f64>,
}

// ---------------------------------------------------------------------------
// Calculator lifecycle
// ---------------------------------------------------------------------------

/// Create a new calculator containing every built-in descriptor.
///
/// The returned pointer must be freed with `mordred_calculator_free`.
#[unsafe(no_mangle)]
pub extern "C" fn mordred_calculator_new() -> *mut MordredCalculator {
    let inner = DescriptorSet::all();
    let cached_names: Vec<CString> = inner
        .names()
        .into_iter()
        .map(|n| CString::new(n).expect("descriptor name contains interior NUL"))
        .collect();

    let calc = MordredCalculator {
        inner,
        cached_names,
    };
    Box::into_raw(Box::new(calc))
}

/// Free a calculator previously created by `mordred_calculator_new`.
///
/// Passing a null pointer is a no-op.
///
/// # Safety
///
/// `calc` must be a pointer returned by `mordred_calculator_new` and must not
/// have been freed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_calculator_free(calc: *mut MordredCalculator) {
    if !calc.is_null() {
        unsafe {
            drop(Box::from_raw(calc));
        }
    }
}

// ---------------------------------------------------------------------------
// Molecule lifecycle
// ---------------------------------------------------------------------------

/// Parse a SMILES string into a molecule handle.
///
/// On success a non-null pointer is returned and `*error_out` is left
/// unchanged.  On failure null is returned and, if `error_out` is non-null,
/// `*error_out` is set to a heap-allocated error string that must be freed
/// with `mordred_string_free`.
///
/// # Safety
///
/// `smiles` must be a valid null-terminated C string.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_molecule_from_smiles(
    smiles: *const c_char,
    error_out: *mut *mut c_char,
) -> *mut MordredMolecule {
    if smiles.is_null() {
        if !error_out.is_null() {
            let msg = CString::new("smiles pointer is null").unwrap();
            unsafe {
                *error_out = msg.into_raw();
            }
        }
        return ptr::null_mut();
    }

    let c_str = unsafe { CStr::from_ptr(smiles) };
    let smiles_str = match c_str.to_str() {
        Ok(s) => s,
        Err(e) => {
            if !error_out.is_null() {
                let msg = CString::new(format!("invalid UTF-8 in SMILES: {e}")).unwrap();
                unsafe {
                    *error_out = msg.into_raw();
                }
            }
            return ptr::null_mut();
        }
    };

    match parse_smiles(smiles_str) {
        Ok(mol) => Box::into_raw(Box::new(MordredMolecule { inner: mol })),
        Err(e) => {
            if !error_out.is_null() {
                let msg = CString::new(format!("{e}")).unwrap();
                unsafe {
                    *error_out = msg.into_raw();
                }
            }
            ptr::null_mut()
        }
    }
}

/// Free a molecule previously created by `mordred_molecule_from_smiles`.
///
/// Passing a null pointer is a no-op.
///
/// # Safety
///
/// `mol` must be a pointer returned by `mordred_molecule_from_smiles` and must
/// not have been freed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_molecule_free(mol: *mut MordredMolecule) {
    if !mol.is_null() {
        unsafe {
            drop(Box::from_raw(mol));
        }
    }
}

// ---------------------------------------------------------------------------
// Calculation
// ---------------------------------------------------------------------------

/// Calculate all descriptors for a molecule.
///
/// Returns a result handle that must be freed with `mordred_result_free`.
/// Returns null if either `calc` or `mol` is null.
///
/// # Safety
///
/// Both `calc` and `mol` must be valid, non-freed pointers obtained from the
/// corresponding constructor functions.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_calculate(
    calc: *const MordredCalculator,
    mol: *const MordredMolecule,
) -> *mut MordredResult {
    if calc.is_null() || mol.is_null() {
        return ptr::null_mut();
    }

    let calc = unsafe { &*calc };
    let mol = unsafe { &*mol };

    let raw_results = calc.inner.calculate(&mol.inner);

    let entries: Vec<MordredResultEntry> = raw_results
        .into_iter()
        .map(|(name, result)| {
            let cname = CString::new(name).expect("descriptor name contains interior NUL");
            MordredResultEntry {
                name: cname,
                value: result.ok(),
            }
        })
        .collect();

    Box::into_raw(Box::new(MordredResult { entries }))
}

/// Free a result set previously returned by `mordred_calculate`.
///
/// Passing a null pointer is a no-op.
///
/// # Safety
///
/// `result` must be a pointer returned by `mordred_calculate` and must not
/// have been freed already.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_result_free(result: *mut MordredResult) {
    if !result.is_null() {
        unsafe {
            drop(Box::from_raw(result));
        }
    }
}

// ---------------------------------------------------------------------------
// Result accessors
// ---------------------------------------------------------------------------

/// Return the number of entries in the result set.
///
/// Returns 0 if `result` is null.
///
/// # Safety
///
/// `result` must be a valid, non-freed pointer returned by `mordred_calculate`,
/// or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_result_len(result: *const MordredResult) -> usize {
    if result.is_null() {
        return 0;
    }
    unsafe { &*result }.entries.len()
}

/// Return the descriptor name for the entry at `index`.
///
/// The returned pointer is valid until `mordred_result_free` is called.
/// Returns null if `result` is null or `index` is out of bounds.
///
/// # Safety
///
/// `result` must be a valid, non-freed pointer returned by `mordred_calculate`,
/// or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_result_name(
    result: *const MordredResult,
    index: usize,
) -> *const c_char {
    if result.is_null() {
        return ptr::null();
    }
    let r = unsafe { &*result };
    match r.entries.get(index) {
        Some(entry) => entry.name.as_ptr(),
        None => ptr::null(),
    }
}

/// Return whether the calculation for the entry at `index` succeeded.
///
/// Returns `false` if `result` is null or `index` is out of bounds.
///
/// # Safety
///
/// `result` must be a valid, non-freed pointer returned by `mordred_calculate`,
/// or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_result_is_ok(result: *const MordredResult, index: usize) -> bool {
    if result.is_null() {
        return false;
    }
    let r = unsafe { &*result };
    match r.entries.get(index) {
        Some(entry) => entry.value.is_some(),
        None => false,
    }
}

/// Return the descriptor value for the entry at `index`.
///
/// Returns `NaN` if the calculation failed, `result` is null, or `index` is
/// out of bounds.
///
/// # Safety
///
/// `result` must be a valid, non-freed pointer returned by `mordred_calculate`,
/// or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_result_value(result: *const MordredResult, index: usize) -> f64 {
    if result.is_null() {
        return f64::NAN;
    }
    let r = unsafe { &*result };
    match r.entries.get(index) {
        Some(entry) => entry.value.unwrap_or(f64::NAN),
        None => f64::NAN,
    }
}

// ---------------------------------------------------------------------------
// Descriptor listing on the calculator
// ---------------------------------------------------------------------------

/// Return the number of descriptors in the calculator.
///
/// Returns 0 if `calc` is null.
///
/// # Safety
///
/// `calc` must be a valid, non-freed pointer returned by
/// `mordred_calculator_new`, or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_calculator_descriptor_count(
    calc: *const MordredCalculator,
) -> usize {
    if calc.is_null() {
        return 0;
    }
    unsafe { &*calc }.cached_names.len()
}

/// Return the name of the descriptor at `index`.
///
/// The returned pointer is valid until `mordred_calculator_free` is called.
/// Returns null if `calc` is null or `index` is out of bounds.
///
/// # Safety
///
/// `calc` must be a valid, non-freed pointer returned by
/// `mordred_calculator_new`, or null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_calculator_descriptor_name(
    calc: *const MordredCalculator,
    index: usize,
) -> *const c_char {
    if calc.is_null() {
        return ptr::null();
    }
    let c = unsafe { &*calc };
    match c.cached_names.get(index) {
        Some(name) => name.as_ptr(),
        None => ptr::null(),
    }
}

// ---------------------------------------------------------------------------
// String cleanup
// ---------------------------------------------------------------------------

/// Free a string previously allocated by mordred (e.g. an error message set by
/// `mordred_molecule_from_smiles`).
///
/// Passing a null pointer is a no-op.
///
/// # Safety
///
/// `s` must be a pointer previously returned via an error-out parameter, or
/// null.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn mordred_string_free(s: *mut c_char) {
    if !s.is_null() {
        unsafe {
            drop(CString::from_raw(s));
        }
    }
}
