/// SMILES parser — supports a practical subset for descriptor calculation.
///
/// Supports: organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I),
/// bracket atoms [Fe], bonds (-, =, #, :), branches, ring closures,
/// aromatic atoms (c, n, o, s, p), charges, and hydrogen counts.
use petgraph::graph::NodeIndex;

use super::aromaticity::perceive_aromaticity;
use super::atom::Atom;
use super::bond::{Bond, BondOrder};
use super::element::Element;
use super::mol::Molecule;
use crate::error::MordredError;

/// Parse a SMILES string into a Molecule.
pub fn parse_smiles(smiles: &str) -> Result<Molecule, MordredError> {
    let mut parser = SmilesParser::new(smiles);
    parser.parse()?;
    parser.fill_implicit_hydrogens();
    perceive_aromaticity(&mut parser.mol);
    Ok(parser.mol)
}

struct SmilesParser {
    chars: Vec<char>,
    pos: usize,
    mol: Molecule,
    stack: Vec<NodeIndex>,
    ring_opens: std::collections::HashMap<u8, (NodeIndex, Option<BondOrder>)>,
    pending_bond: Option<BondOrder>,
}

impl SmilesParser {
    fn new(smiles: &str) -> Self {
        Self {
            chars: smiles.chars().collect(),
            pos: 0,
            mol: Molecule::new(),
            stack: Vec::new(),
            ring_opens: std::collections::HashMap::new(),
            pending_bond: None,
        }
    }

    fn peek(&self) -> Option<char> {
        self.chars.get(self.pos).copied()
    }

    fn advance(&mut self) -> Option<char> {
        let c = self.chars.get(self.pos).copied();
        if c.is_some() {
            self.pos += 1;
        }
        c
    }

    fn parse(&mut self) -> Result<(), MordredError> {
        while self.pos < self.chars.len() {
            match self.peek() {
                Some('(') => {
                    self.advance();
                    // Push current atom onto branch stack
                    if let Some(&top) = self.stack.last() {
                        self.stack.push(top);
                    }
                }
                Some(')') => {
                    self.advance();
                    self.stack.pop();
                }
                Some('-') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Single);
                }
                Some('=') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Double);
                }
                Some('#') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Triple);
                }
                Some(':') => {
                    self.advance();
                    self.pending_bond = Some(BondOrder::Aromatic);
                }
                Some('[') => {
                    let atom = self.parse_bracket_atom()?;
                    self.add_atom_and_bond(atom);
                }
                Some(c) if c.is_ascii_digit() => {
                    self.advance();
                    let ring_num = c as u8 - b'0';
                    self.handle_ring_closure(ring_num)?;
                }
                Some('%') => {
                    // Two-digit ring closure
                    self.advance();
                    let d1 = self.advance().ok_or(MordredError::SmilesParseError(
                        "unexpected end after %".into(),
                    ))?;
                    let d2 = self.advance().ok_or(MordredError::SmilesParseError(
                        "unexpected end after %digit".into(),
                    ))?;
                    let ring_num = (d1 as u8 - b'0') * 10 + (d2 as u8 - b'0');
                    self.handle_ring_closure(ring_num)?;
                }
                Some('/') | Some('\\') => {
                    self.advance();
                    // Cis/trans bond stereo — skip (no-op for 2D descriptors)
                }
                Some('.') => {
                    self.advance();
                    // Disconnected fragment — reset stack connection
                    self.stack.clear();
                }
                Some(c) => {
                    let atom = self.parse_organic_atom(c)?;
                    self.add_atom_and_bond(atom);
                }
                None => break,
            }
        }
        Ok(())
    }

    fn parse_organic_atom(&mut self, c: char) -> Result<Atom, MordredError> {
        self.advance();
        let is_aromatic = c.is_ascii_lowercase();
        let symbol = if is_aromatic {
            c.to_ascii_uppercase().to_string()
        } else {
            // Check for two-letter elements
            let mut s = c.to_string();
            if let Some(next) = self.peek() {
                if next.is_ascii_lowercase() && !is_aromatic_organic(next) {
                    let two = format!("{}{}", c, next);
                    if Element::from_symbol(&two).is_some() {
                        self.advance();
                        s = two;
                    }
                }
            }
            s
        };

        let element = Element::from_symbol(&symbol).ok_or_else(|| {
            MordredError::SmilesParseError(format!("unknown element: {}", symbol))
        })?;

        let mut atom = Atom::new(element);
        atom.is_aromatic = is_aromatic;
        Ok(atom)
    }

    fn parse_bracket_atom(&mut self) -> Result<Atom, MordredError> {
        self.advance(); // consume '['

        let mut isotope = None;
        // Parse optional isotope
        if self.peek().is_some_and(|c| c.is_ascii_digit()) {
            let mut num = 0u16;
            while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                num = num * 10 + (self.advance().unwrap() as u16 - '0' as u16);
            }
            isotope = Some(num);
        }

        // Parse element symbol
        let first = self.advance().ok_or(MordredError::SmilesParseError(
            "unexpected end in bracket atom".into(),
        ))?;

        let is_aromatic = first.is_ascii_lowercase();
        let mut symbol = if is_aromatic {
            first.to_ascii_uppercase().to_string()
        } else {
            first.to_string()
        };

        // Check for second lowercase letter (two-letter element)
        if self.peek().is_some_and(|c| c.is_ascii_lowercase()) {
            let two = format!("{}{}", symbol.chars().next().unwrap(), self.peek().unwrap());
            if Element::from_symbol(&two).is_some() {
                self.advance();
                symbol = two;
            }
        }

        let element = Element::from_symbol(&symbol).ok_or_else(|| {
            MordredError::SmilesParseError(format!("unknown element: {}", symbol))
        })?;

        let mut atom = Atom::new(element);
        atom.is_aromatic = is_aromatic;
        atom.isotope = isotope;

        // Parse optional H count
        if self.peek() == Some('H') {
            self.advance();
            if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                atom.implicit_h = self.advance().unwrap() as u8 - b'0';
            } else {
                atom.implicit_h = 1;
            }
        }

        // Parse optional charge
        if self.peek() == Some('+') {
            self.advance();
            if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                atom.charge = (self.advance().unwrap() as u8 - b'0') as i8;
            } else if self.peek() == Some('+') {
                self.advance();
                atom.charge = 2;
            } else {
                atom.charge = 1;
            }
        } else if self.peek() == Some('-') {
            self.advance();
            if self.peek().is_some_and(|c| c.is_ascii_digit()) {
                atom.charge = -((self.advance().unwrap() as u8 - b'0') as i8);
            } else if self.peek() == Some('-') {
                self.advance();
                atom.charge = -2;
            } else {
                atom.charge = -1;
            }
        }

        // Skip to closing bracket (ignore stereo, etc.)
        while self.peek() != Some(']') && self.peek().is_some() {
            self.advance();
        }
        self.advance(); // consume ']'

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

        // Replace top of stack with current atom
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

        if let Some((open_idx, open_bond)) = self.ring_opens.remove(&ring_num) {
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
            self.ring_opens
                .insert(ring_num, (current, self.pending_bond.take()));
        }
        Ok(())
    }

    /// Fill implicit hydrogen counts for organic-subset atoms that don't have them set.
    fn fill_implicit_hydrogens(&mut self) {
        let indices: Vec<NodeIndex> = self.mol.graph.node_indices().collect();
        for idx in indices {
            let atom = &self.mol.graph[idx];
            // Only fill for atoms that came from organic subset (not bracket atoms)
            // Bracket atoms already have explicit H counts.
            // We detect bracket atoms by checking if implicit_h was already set or isotope is present.
            if atom.isotope.is_some() {
                continue; // bracket atom, already handled
            }

            let valence = atom.element.default_valence();
            if valence == 0 {
                continue; // can't compute implicit H
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
                // Aromatic atoms have an additional bond from the pi system
                // that isn't represented as an explicit edge.
                // e.g. aromatic C has valence 4, 2 aromatic bonds + 1 pi = 3, so 1 implicit H.
                atom.implicit_h = valence.saturating_sub(bond_valence + 1);
            }
        }
    }
}

fn is_aromatic_organic(c: char) -> bool {
    matches!(c, 'c' | 'n' | 'o' | 's' | 'p')
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
