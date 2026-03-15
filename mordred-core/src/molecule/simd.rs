/// SIMD-friendly molecular property computation.
///
/// Provides functions that extract SOA (Structure of Arrays) data from
/// a molecule's graph to enable LLVM auto-vectorization of reductions.
///
/// On stable Rust, element discriminants are extracted to a contiguous
/// `Vec<u8>` and processed with tight loops that LLVM can auto-vectorize.
///
/// On nightly Rust with `--features simd`, explicit `std::simd` intrinsics
/// are used for wider vectorization (not yet implemented — reserved for
/// future nightly builds).
use super::element::Element;

/// Build element histogram from a contiguous array of discriminant indices.
///
/// This is structured as a tight loop over a contiguous u8 slice,
/// which LLVM can optimize better than scattered graph accesses.
#[inline]
pub fn element_histogram(discriminants: &[u8]) -> [u32; 26] {
    let mut counts = [0u32; 26];
    for &d in discriminants {
        counts[d as usize] += 1;
    }
    counts
}

/// Sum implicit hydrogen values from a contiguous array.
///
/// Tight u8→u32 reduction loop for auto-vectorization.
#[inline]
pub fn sum_implicit_h(implicit_h_values: &[u8]) -> u32 {
    implicit_h_values.iter().map(|&h| h as u32).sum()
}

/// Accumulate molecular weight from contiguous mass values.
///
/// f64 accumulation in a tight loop enables LLVM to use SIMD FP adds.
#[inline]
pub fn sum_molecular_weight(masses: &[f64]) -> f64 {
    masses.iter().sum()
}

/// Count heavy atoms (non-hydrogen) from discriminant array.
#[inline]
pub fn count_heavy(discriminants: &[u8]) -> u32 {
    let h_idx = Element::H.discriminant_index() as u8;
    discriminants.iter().filter(|&&d| d != h_idx).count() as u32
}

/// Count halogens from discriminant array.
#[inline]
pub fn count_halogens(discriminants: &[u8]) -> u32 {
    let f_idx = Element::F.discriminant_index() as u8;
    let cl_idx = Element::Cl.discriminant_index() as u8;
    let br_idx = Element::Br.discriminant_index() as u8;
    let i_idx = Element::I.discriminant_index() as u8;
    discriminants
        .iter()
        .filter(|&&d| d == f_idx || d == cl_idx || d == br_idx || d == i_idx)
        .count() as u32
}

/// Count heteroatoms (non-C, non-H) from discriminant array.
#[inline]
pub fn count_heteroatoms(discriminants: &[u8]) -> u32 {
    let c_idx = Element::C.discriminant_index() as u8;
    let h_idx = Element::H.discriminant_index() as u8;
    discriminants
        .iter()
        .filter(|&&d| d != c_idx && d != h_idx)
        .count() as u32
}
