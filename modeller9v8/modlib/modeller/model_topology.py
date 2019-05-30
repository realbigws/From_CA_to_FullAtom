import _modeller
from modeller.information import info
from modeller.util import modlist

class BondList(modlist.FixList):
    def __init__(self, mdl, getdimfunc, getfunc, natm):
        self.__mdl = mdl
        self.__mtp = _modeller.mod_model_mtp_get(mdl.modpt)
        self.__getdimfunc = getdimfunc
        self.__getfunc = getfunc
        self.__natm = natm
        modlist.FixList.__init__(self)

    def __len__(self):
        return self.__getdimfunc(self.__mtp)

    def _getfunc(self, indx):
        indarray = self.__getfunc(self.__mtp)
        inds = [_modeller.mod_int2_get(indarray, i, indx) \
                for i in range(self.__natm)]
        ats = self.__mdl.atoms
        return [ats[i-1] for i in inds]


def parse_atoms(atom_names):
    """Parse atom names for + or - prefixes"""
    parse_ats = []
    resoffs = []
    for a in atom_names:
        resoff = 0
        if a.startswith('+'):
            resoffs.append(1)
            parse_ats.append(a[1:])
        elif a.startswith('-'):
            resoffs.append(-1)
            parse_ats.append(a[1:])
        else:
            resoffs.append(0)
            parse_ats.append(a)
    return parse_ats, resoffs


class FindAtoms(object):
    def __init__(self, mdl, residue_type, atom_names, check_func=None, *args):
        self.__mdl = mdl
        self.residue_type = residue_type
        (self.atom_names, self.residue_offsets) = parse_atoms(atom_names)
        self.__index = 1
        self.check_func = check_func
        self.check_args = args
    def __iter__(self):
        return self
    def next(self):
        while True:
            atoms = self.__int_next()
            if self.check_func is None \
               or self.check_func(atoms, *self.check_args):
                return atoms
    def __int_next(self):
        (residue_index, atom_indices) = \
            _modeller.mod_model_find_atoms(self.__mdl.modpt, self.residue_type,
                                           self.atom_names,
                                           self.residue_offsets,
                                           self.__index,
                                           self.__mdl.env.libs.modpt)
        if residue_index == 0:
            raise StopIteration
        else:
            self.__index = residue_index + 1
            return [self.__mdl.atoms[i-1] for i in atom_indices]


class FindDihedrals(object):
    def __init__(self, mdl, residue_type, dihedral_type, check_func=None,
                 *args):
        self.__mdl = mdl
        self.residue_type = residue_type
        self.dihedral_type = dihedral_type
        self.__index = 1
        self.check_func = check_func
        self.check_args = args
    def __iter__(self):
        return self
    def next(self):
        while True:
            atoms = self.__int_next()
            if self.check_func is None \
               or self.check_func(atoms, *self.check_args):
                return atoms
    def __int_next(self):
        (residue_index, atom_indices) = \
            _modeller.mod_model_find_dihedrals(self.__mdl.modpt,
                                               self.residue_type,
                                               self.dihedral_type,
                                               self.__index,
                                               self.__mdl.env.libs.modpt)
        if residue_index == 0:
            raise StopIteration
        else:
            self.__index = residue_index + 1
            return [self.__mdl.atoms[i-1] for i in atom_indices]

def _print_atomlist(fh, mdl, xplor):
    """Utility function for write_psf() - write out atoms in PSF format"""
    print >> fh, "%8d !NATOM" % len(mdl.atoms)
    if xplor:
        for atom in mdl.atoms:
            res = atom.residue
            print >> fh, "%8d %-4s %-4d %-4s %-4s %-4s %10.6f %13.4f %11d" % \
                  (atom.index, res.chain.name, res.index, res.name, atom.name,
                   atom.type.name, atom.charge, atom.type.mass, 0)
    else:
        for atom in mdl.atoms:
            res = atom.residue
            print >> fh, "%8d %-4s %-4d %-4s %-4s %4d %10.6f %13.4f %11d" % \
                  (atom.index, res.chain.name, res.index, res.name, atom.name,
                   atom.type.index, atom.charge, atom.type.mass, 0)

def _print_bondlist(fh, mdl, bondlist, title, num_per_line):
    """Utility function for write_psf() - write out bonds in PSF format"""
    print >> fh, "\n%8d %s" % (len(bondlist), title)
    for start in range(0, len(bondlist), num_per_line):
        line = ""
        for bond in bondlist[start:start+num_per_line]:
            for atom in bond:
                line += "%8d" % atom.index
        print >> fh, line


def write_psf(filename, mdl, xplor=True):
    """Utility function for model.write_psf()"""
    if hasattr(filename, 'write'):
        fh = filename
    else:
        fh = open(filename, 'w')
    print >> fh, "PSF\n"
    print >> fh, "%8d !NTITLE" % 1
    if xplor:
        psftyp = "X-PLOR"
    else:
        psftyp = "CHARMM"
    print >> fh, "* %s type PSF generated by %s\n" % (psftyp, info.version)
    _print_atomlist(fh, mdl, xplor)
    _print_bondlist(fh, mdl, mdl.bonds, "!NBOND: bonds", 4)
    _print_bondlist(fh, mdl, mdl.angles, "!NTHETA: angles", 3)
    _print_bondlist(fh, mdl, mdl.dihedrals, "!NPHI: dihedrals", 2)
    _print_bondlist(fh, mdl, mdl.impropers, "!NIMPHI: impropers", 2)
    print >> fh, "\n%8d !NDON: donors" % 0
    print >> fh, "\n%8d !NACC: acceptors" % 0
