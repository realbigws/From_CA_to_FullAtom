"""Classes to restrain secondary structure"""

import _modeller
from modeller import coordinates

def _parse_residue_list(residues, mdl):
    if not isinstance(residues, coordinates.ResidueList):
        raise TypeError("expecting a residue list, e.g. from " +
                        "model.residue_range()")
    if residues.mdl != mdl:
        raise ValueError("residue list is for a different " +
                         "model %s, %s" % (str(residues.mdl), str(mdl)))
    start = residues.offset + 1
    end = start + len(residues) - 1
    return (start, end)

class alpha(object):
    """Alpha (helix) secondary structure restraint"""
    def __init__(self, residues):
        self.__residues = residues

    def _add_restraint(self, rsr, mdl):
        (start, end) = _parse_residue_list(self.__residues, mdl)
        _modeller.mod_restraints_make_alpha(mdl.modpt, (start, end),
                                            mdl.env.libs.modpt)

class strand(object):
    """Strand secondary structure restraint"""
    def __init__(self, residues):
        self.__residues = residues

    def _add_restraint(self, rsr, mdl):
        (start, end) = _parse_residue_list(self.__residues, mdl)
        _modeller.mod_restraints_make_strand(mdl.modpt, (start, end),
                                             mdl.env.libs.modpt)

class sheet(object):
    """Beta-sheet secondary structure restraint"""
    def __init__(self, atom1, atom2, sheet_h_bonds):
        self.__atoms = (atom1, atom2)
        self.__sheet_h_bonds = sheet_h_bonds

    def _add_restraint(self, rsr, mdl):
        atoms = self.__atoms
        for (n, at) in enumerate(atoms):
            if not isinstance(at, coordinates.Atom):
                raise TypeError("atom%d must be an Atom object" % (n+1))
            if at.mdl != mdl:
                raise ValueError("atom%d is from a different model" % (n+1))
        _modeller.mod_restraints_make_sheet(mdl.modpt, [a.index for a in atoms],
                                            self.__sheet_h_bonds)
