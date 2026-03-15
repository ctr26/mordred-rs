/// SMILES parser — supports a practical subset for descriptor calculation.
///
/// Supports: organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I),
/// bracket atoms [Fe], bonds (-, =, #, :), branches, ring closures,
/// aromatic atoms (c, n, o, s, p), charges, and hydrogen counts.
///
/// Optimized for throughput: operates on byte slices (SMILES is ASCII),
/// uses stack-allocated ring closure array, and pre-sizes the graph.
use petgraph::graph::NodeIndex;

use super::aromaticity::perceive_aromaticity;
use super::atom::Atom;
use super::bond::{Bond, BondOrder};
use super::element::Element;
use super::mol::Molecule;
use crate::error::MordredError;

/// Parse a SMILES string into a Molecule.
pub fn parse_smiles(smiles: &str) -> Result<Molecule, MordredError> {
    let mut parser = SmilesParser::new(smiles.as_bytes());
    parser.parse()?;
    parser.fill_implicit_hydrogens();
    perceive_aromaticity(&mut parser.mol);
    Ok(parser.mol)
}

struct SmilesParser<'a> {
    bytes: &'a [u8],
    pos: usize,
    mol: Molecule,
    stack: Vec<NodeIndex>,
    ring_opens: [Option<(NodeIndex, Option<BondOrder>)>; 100],
    pending_bond: Option<BondOrder>,
}

impl<'a> SmilesParser<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        // Pre-size graph: rough estimate of atoms ~ len/2, bonds ~ len/2
        let est = bytes.len() / 2 + 1;
        let mut mol = Molecule::new();
        mol.graph.reserve_nodes(est);
        mol.graph.reserve_edges(est);

        Self {
            bytes,
            pos: 0,
            mol,
            stack: Vec::with_capacity(est),
            ring_opens: [None; 100],
            pending_bond: None,
        }
    }

    #[inline(always)]
    fn peek(&self) -> Option<u8> {
        self.bytes.get(self.pos).copied()
    }

    #[inline(always)]
    fn advance(&mut self) -> Option<u8> {
        let b = self.bytes.get(self.pos).copied();
        if b.is_some() {
            self.pos += 1;
        }
        b
    }

    fn parse(&mut self) -> Result<(), MordredError> {
        while self.pos < self.bytes.len() {
            match self.bytes[self.pos] {
                b'(' => {
                    self.pos += 1;
                    if let Some(&top) = self.stack.last() {
                        self.stack.push(top);
                    }
                }
                b')' => {
                    self.pos += 1;
                    self.stack.pop();
                }
                b'-' => {
                    self.pos += 1;
                    self.pending_bond = Some(BondOrder::Single);
                }
                b'=' => {
                    self.pos += 1;
                    self.pending_bond = Some(BondOrder::Double);
                }
                b'#' => {
                    self.pos += 1;
                    self.pending_bond = Some(BondOrder::Triple);
                }
                b':' => {
                    self.pos += 1;
                    self.pending_bond = Some(BondOrder::Aromatic);
                }
                b'[' => {
                    let atom = self.parse_bracket_atom()?;
                    self.add_atom_and_bond(atom);
                }
                b'0'..=b'9' => {
                    let ring_num = self.bytes[self.pos] - b'0';
                    self.pos += 1;
                    self.handle_ring_closure(ring_num)?;
                }
                b'%' => {
                    self.pos += 1;
                    let d1 = self.advance().ok_or(MordredError::SmilesParseError(
                        "unexpected end after %".into(),
                    ))?;
                    let d2 = self.advance().ok_or(MordredError::SmilesParseError(
                        "unexpected end after %digit".into(),
                    ))?;
                    let ring_num = (d1 - b'0') * 10 + (d2 - b'0');
                    self.handle_ring_closure(ring_num)?;
                }
                b'/' | b'\\' => {
                    self.pos += 1;
                    // Cis/trans bond stereo — skip (no-op for 2D descriptors)
                }
                b'.' => {
                    self.pos += 1;
                    self.stack.clear();
                }
                c => {
                    let atom = self.parse_organic_atom(c)?;
                    self.add_atom_and_bond(atom);
                }
            }
        }
        Ok(())
    }

    fn parse_organic_atom(&mut self, c: u8) -> Result<Atom, MordredError> {
        self.pos += 1;
        let is_aromatic = c.is_ascii_lowercase();

        let element = if is_aromatic {
            match c {
                b'c' => Element::C,
                b'n' => Element::N,
                b'o' => Element::O,
                b's' => Element::S,
                b'p' => Element::P,
                _ => {
                    return Err(MordredError::SmilesParseError(format!(
                        "unknown aromatic atom: {}",
                        c as char
                    )));
                }
            }
        } else {
            // Check for two-letter elements first
            match (c, self.peek()) {
                (b'B', Some(b'r')) => {
                    self.pos += 1;
                    Element::Br
                }
                (b'C', Some(b'l')) => {
                    self.pos += 1;
                    Element::Cl
                }
                (b'S', Some(b'i')) => {
                    self.pos += 1;
                    Element::Si
                }
                (b'S', Some(b'e')) => {
                    self.pos += 1;
                    Element::Se
                }
                (b'N', Some(b'a')) => {
                    self.pos += 1;
                    Element::Na
                }
                (b'A', Some(b'l')) => {
                    self.pos += 1;
                    Element::Al
                }
                // Single letter elements
                (b'B', _) => Element::B,
                (b'C', _) => Element::C,
                (b'N', _) => Element::N,
                (b'O', _) => Element::O,
                (b'P', _) => Element::P,
                (b'S', _) => Element::S,
                (b'F', _) => Element::F,
                (b'I', _) => Element::I,
                _ => {
                    return Err(MordredError::SmilesParseError(format!(
                        "unknown element: {}",
                        c as char
                    )));
                }
            }
        };

        let mut atom = Atom::new(element);
        atom.is_aromatic = is_aromatic;
        Ok(atom)
    }

    fn parse_bracket_atom(&mut self) -> Result<Atom, MordredError> {
        self.pos += 1; // consume '['

        let mut isotope = None;
        // Parse optional isotope
        while self.peek().is_some_and(|b| b.is_ascii_digit()) {
            let num = isotope.unwrap_or(0u16);
            isotope = Some(num * 10 + (self.advance().unwrap() - b'0') as u16);
        }

        // Parse element symbol
        let first = self.advance().ok_or(MordredError::SmilesParseError(
            "unexpected end in bracket atom".into(),
        ))?;

        let is_aromatic = first.is_ascii_lowercase();
        let upper = if is_aromatic {
            first.to_ascii_uppercase()
        } else {
            first
        };

        // Check for two-letter element
        let element = if let Some(next) = self.peek() {
            if next.is_ascii_lowercase() {
                let sym = [upper, next];
                if let Some(elem) = element_from_two_bytes(sym) {
                    self.pos += 1;
                    elem
                } else {
                    element_from_byte(upper)?
                }
            } else {
                element_from_byte(upper)?
            }
        } else {
            element_from_byte(upper)?
        };

        let mut atom = Atom::new(element);
        atom.is_aromatic = is_aromatic;
        atom.isotope = isotope;

        // Parse optional H count
        if self.peek() == Some(b'H') {
            self.pos += 1;
            if self.peek().is_some_and(|b| b.is_ascii_digit()) {
                atom.implicit_h = self.advance().unwrap() - b'0';
            } else {
                atom.implicit_h = 1;
            }
        }

        // Parse optional charge
        if self.peek() == Some(b'+') {
            self.pos += 1;
            if self.peek().is_some_and(|b| b.is_ascii_digit()) {
                atom.charge = (self.advance().unwrap() - b'0') as i8;
            } else if self.peek() == Some(b'+') {
                self.pos += 1;
                atom.charge = 2;
            } else {
                atom.charge = 1;
            }
        } else if self.peek() == Some(b'-') {
            self.pos += 1;
            if self.peek().is_some_and(|b| b.is_ascii_digit()) {
                atom.charge = -((self.advance().unwrap() - b'0') as i8);
            } else if self.peek() == Some(b'-') {
                self.pos += 1;
                atom.charge = -2;
            } else {
                atom.charge = -1;
            }
        }

        // Skip to closing bracket (ignore stereo, etc.)
        while self.peek() != Some(b']') && self.peek().is_some() {
            self.pos += 1;
        }
        self.pos += 1; // consume ']'

        Ok(atom)
    }

    fn add_atom_and_bond(&mut self, atom: Atom) {
        let is_aromatic = atom.is_aromatic;
        let idx = self.mol.add_atom(atom);

        if let Some(&prev) = self.stack.last() {
            let order = self.pending_bond.take().unwrap_or_else(|| {
                if is_aromatic && self.mol.graph[prev].is_aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            });
            self.mol.add_bond(prev, idx, Bond::new(order));
        }
        self.pending_bond = None;

        if self.stack.is_empty() {
            self.stack.push(idx);
        } else {
            *self.stack.last_mut().unwrap() = idx;
        }
    }

    fn handle_ring_closure(&mut self, ring_num: u8) -> Result<(), MordredError> {
        let current = *self.stack.last().ok_or(MordredError::SmilesParseError(
            "ring closure with no current atom".into(),
        ))?;

        let idx = ring_num as usize;
        if let Some((open_idx, open_bond)) = self.ring_opens[idx].take() {
            let order = self.pending_bond.take().or(open_bond).unwrap_or_else(|| {
                if self.mol.graph[current].is_aromatic && self.mol.graph[open_idx].is_aromatic {
                    BondOrder::Aromatic
                } else {
                    BondOrder::Single
                }
            });
            let mut bond = Bond::new(order);
            bond.is_ring = true;
            self.mol.add_bond(open_idx, current, bond);
        } else {
            self.ring_opens[idx] = Some((current, self.pending_bond.take()));
        }
        Ok(())
    }

    fn fill_implicit_hydrogens(&mut self) {
        let indices: Vec<NodeIndex> = self.mol.graph.node_indices().collect();
        for idx in indices {
            let atom = &self.mol.graph[idx];
            if atom.isotope.is_some() {
                continue;
            }

            let valence = atom.element.default_valence();
            if valence == 0 {
                continue;
            }

            let bond_valence: u8 = self
                .mol
                .graph
                .edges(idx)
                .map(|e| e.weight().order.valence_contribution())
                .sum();

            let atom = &mut self.mol.graph[idx];
            if atom.implicit_h == 0 && !atom.is_aromatic {
                atom.implicit_h = valence.saturating_sub(bond_valence);
            } else if atom.is_aromatic && atom.implicit_h == 0 {
                atom.implicit_h = valence.saturating_sub(bond_valence + 1);
            }
        }
    }
}

#[inline]
fn element_from_byte(b: u8) -> Result<Element, MordredError> {
    match b {
        b'H' => Ok(Element::H),
        b'B' => Ok(Element::B),
        b'C' => Ok(Element::C),
        b'N' => Ok(Element::N),
        b'O' => Ok(Element::O),
        b'F' => Ok(Element::F),
        b'P' => Ok(Element::P),
        b'S' => Ok(Element::S),
        b'I' => Ok(Element::I),
        b'K' => Ok(Element::K),
        _ => Err(MordredError::SmilesParseError(format!(
            "unknown element: {}",
            b as char
        ))),
    }
}

#[inline]
fn element_from_two_bytes(sym: [u8; 2]) -> Option<Element> {
    match &sym {
        b"He" => Some(Element::He),
        b"Li" => Some(Element::Li),
        b"Be" => Some(Element::Be),
        b"Ne" => Some(Element::Ne),
        b"Na" => Some(Element::Na),
        b"Mg" => Some(Element::Mg),
        b"Al" => Some(Element::Al),
        b"Si" => Some(Element::Si),
        b"Cl" => Some(Element::Cl),
        b"Ar" => Some(Element::Ar),
        b"Ca" => Some(Element::Ca),
        b"Br" => Some(Element::Br),
        b"Fe" => Some(Element::Fe),
        b"Cu" => Some(Element::Cu),
        b"Zn" => Some(Element::Zn),
        b"Se" => Some(Element::Se),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_methane() {
        let mol = parse_smiles("C").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.graph[NodeIndex::new(0)].implicit_h, 4);
    }

    #[test]
    fn test_parse_ethanol() {
        let mol = parse_smiles("CCO").unwrap();
        assert_eq!(mol.atom_count(), 3);
        assert_eq!(mol.bond_count(), 2);
    }

    #[test]
    fn test_parse_benzene() {
        let mol = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(mol.atom_count(), 6);
        assert_eq!(mol.bond_count(), 6);
    }

    #[test]
    fn test_parse_double_bond() {
        let mol = parse_smiles("C=O").unwrap();
        assert_eq!(mol.atom_count(), 2);
        assert_eq!(mol.bond_count(), 1);
    }

    #[test]
    fn test_parse_branch() {
        let mol = parse_smiles("CC(=O)O").unwrap(); // acetic acid
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
    }

    #[test]
    fn test_parse_bracket_atom() {
        let mol = parse_smiles("[Fe]").unwrap();
        assert_eq!(mol.atom_count(), 1);
        assert_eq!(mol.graph[NodeIndex::new(0)].element, Element::Fe);
    }

    #[test]
    fn test_parse_charged() {
        let mol = parse_smiles("[NH4+]").unwrap();
        assert_eq!(mol.graph[NodeIndex::new(0)].charge, 1);
        assert_eq!(mol.graph[NodeIndex::new(0)].implicit_h, 4);
    }

    #[test]
    fn test_parse_stereo_cis_trans() {
        // cis-2-butene: C/C=C\C
        let mol = parse_smiles("C/C=C\\C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
    }

    #[test]
    fn test_parse_stereo_slash_forward() {
        // trans-2-butene: C/C=C/C
        let mol = parse_smiles("C/C=C/C").unwrap();
        assert_eq!(mol.atom_count(), 4);
        assert_eq!(mol.bond_count(), 3);
    }

    #[test]
    fn test_parse_chiral() {
        // L-alanine
        let mol = parse_smiles("[C@@H](N)(C)C(=O)O").unwrap();
        assert!(mol.atom_count() > 0);
    }

    #[test]
    fn test_parse_chiral_single_at() {
        // Chiral center with single @
        let mol = parse_smiles("[C@H](F)(Cl)Br").unwrap();
        assert_eq!(mol.atom_count(), 4);
    }
}
