"""Classes for handling protein structure coordinates"""

import _modeller
import modeller.selection as selection
from modeller import sequence, modfile
import modeller.util.modutil as modutil
from modeller.util.logger import log

__docformat__ = "epytext en"

class Point(object):
    """An arbitary point in the Cartesian space of a model"""
    def __init__(self, mdl, x, y, z):
        self.mdl = mdl
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "<Point (%.2f, %.2f, %.2f)>" % (self.x, self.y, self.z)

    def select_sphere(self, radius):
        """Returns a selection of all atoms within the given distance"""
        inds = _modeller.mod_selection_sphere(self.mdl.modpt, self.x, self.y,
                                              self.z, radius)
        return selection.selection(AtomIndices(inds, self.mdl))


class Atom(Point):
    """A single atom in a protein structure"""

    def __init__(self, mdl, num):
        self.mdl = mdl
        self._num = num

    def get_atom_indices(self):
        return [self._num + 1], self.mdl

    def __repr__(self):
        chainid = self.residue.get_chain_suffix()
        return "<Atom %s:%s%s>" % (self.name, self.residue.num, chainid)

    def __get_x(self):
        x = _modeller.mod_coordinates_x_get(self.mdl.cdpt)
        return _modeller.mod_float1_get(x, self._num)
    def __set_x(self, val):
        x = _modeller.mod_coordinates_x_get(self.mdl.cdpt)
        _modeller.mod_float1_set(x, self._num, val)
    def __get_y(self):
        y = _modeller.mod_coordinates_y_get(self.mdl.cdpt)
        return _modeller.mod_float1_get(y, self._num)
    def __set_y(self, val):
        y = _modeller.mod_coordinates_y_get(self.mdl.cdpt)
        _modeller.mod_float1_set(y, self._num, val)
    def __get_z(self):
        z = _modeller.mod_coordinates_z_get(self.mdl.cdpt)
        return _modeller.mod_float1_get(z, self._num)
    def __set_z(self, val):
        z = _modeller.mod_coordinates_z_get(self.mdl.cdpt)
        _modeller.mod_float1_set(z, self._num, val)
    def __get_biso(self):
        biso = _modeller.mod_coordinates_biso_get(self.mdl.cdpt)
        return _modeller.mod_float1_get(biso, self._num)
    def __set_biso(self, val):
        biso = _modeller.mod_coordinates_biso_get(self.mdl.cdpt)
        _modeller.mod_float1_set(biso, self._num, val)
    def __get_occ(self):
        occ = _modeller.mod_coordinates_occ_get(self.mdl.cdpt)
        return _modeller.mod_float1_get(occ, self._num)
    def __set_occ(self, val):
        occ = _modeller.mod_coordinates_occ_get(self.mdl.cdpt)
        _modeller.mod_float1_set(occ, self._num, val)
    def __get_name(self):
        return _modeller.mod_coordinates_atmnam_get(self.mdl.cdpt, self._num)
    def __get_element(self):
        return _modeller.mod_coordinates_element_get(self.mdl.cdpt, self._num)
    def __set_element(self, val):
        return _modeller.mod_coordinates_element_set(self.mdl.cdpt, self._num,
                                                     val)
    def __get_residue(self):
        iresatm = _modeller.mod_coordinates_iresatm_get(self.mdl.cdpt)
        resind = _modeller.mod_int1_get(iresatm, self._num) - 1
        return self.mdl.residues[resind]
    def __get_index(self):
        return self._num + 1

    x = property(__get_x, __set_x, doc="x coordinate")
    y = property(__get_y, __set_y, doc="y coordinate")
    z = property(__get_z, __set_z, doc="z coordinate")
    biso = property(__get_biso, __set_biso, doc="Isotropic temperature factor")
    occ = property(__get_occ, __set_occ, doc="Occupancy")
    name = property(__get_name, doc="PDB name")
    element = property(__get_element, __set_element, doc="Element symbol")
    residue = property(__get_residue, doc="Residue object containing this atom")
    index = property(__get_index, doc="Internal atom index")


class Chain(sequence.Chain):
    """A single chain/segment in a protein structure."""

    def write(self, file, atom_file, align_code, comment='', format='PIR',
              chop_nonstd_termini=True):
        """Write this chain out to an alignment file"""
        f = _modeller.mod_chain_write
        if isinstance(file, str):
            file = modfile.File(file, 'w')
        return f(self.seq.seqpt, self.seq.cdpt, self._num, file.file_pointer,
                 atom_file, align_code, comment, format, chop_nonstd_termini,
                 self.seq.env.libs.modpt)

    def __get_atoms(self):
        (startres, endres) = self._get_resind()
        (startatm, endatm) = get_residue_atom_indices(self.seq, startres,
                                                      endres)
        suffix = ":%s" % self.name
        return AtomList(self.seq, startatm, endatm - startatm, suffix)

    def get_atom_indices(self):
        (startres, endres) = self._get_resind()
        (startatm, endatm) = get_residue_atom_indices(self.seq, startres,
                                                      endres)
        return (range(startatm+1, endatm+1), self.seq)

    atoms = property(__get_atoms, doc="List of all atoms in this chain")


class Residue(sequence.Residue):
    """A single residue in a protein structure"""

    def get_atom_indices(self):
        (startind, endind) = get_residue_atom_indices(self.mdl, self._num,
                                                      self._num + 1)
        return range(startind + 1, endind + 1), self.mdl

    def __repr__(self):
        chainid = self.get_chain_suffix()
        return "Residue %s%s (type %s)" % (self.num, chainid, self.name)

    def __str__(self):
        return "<%s>" % repr(self)

    def get_chain_suffix(self):
        """Returns a suffix - e.g. ':A' - to identify the chain this residue
           is in. If the chain has no ID, returns a blank string."""
        chainid = self.chain.name
        if chainid == ' ':
            return ''
        else:
            return ':' + chainid

    def __get_num(self):
        return _modeller.mod_coordinates_resnum_get(self.mdl.cdpt,
                                                    self._num).strip()
    def __set_num(self, val):
        _modeller.mod_coordinates_resnum_set(self.mdl.cdpt, self._num, val)
    def __get_atoms(self):
        (startind, endind) = get_residue_atom_indices(self.mdl, self._num,
                                                      self._num + 1)
        suffix = ":%s%s" % (self.num, self.get_chain_suffix())
        return AtomList(self.mdl, startind, endind - startind, suffix)

    def __get_alpha(self):
        return get_dihedral("alpha", 1, self._num, self.mdl)
    def __get_phi(self):
        return get_dihedral("phi", 2, self._num, self.mdl)
    def __get_psi(self):
        return get_dihedral("psi", 3, self._num, self.mdl)
    def __get_omega(self):
        return get_dihedral("omega", 4, self._num, self.mdl)
    def __get_chi1(self):
        return get_dihedral("chi1", 5, self._num, self.mdl)
    def __get_chi2(self):
        return get_dihedral("chi2", 6, self._num, self.mdl)
    def __get_chi3(self):
        return get_dihedral("chi3", 7, self._num, self.mdl)
    def __get_chi4(self):
        return get_dihedral("chi4", 8, self._num, self.mdl)
    def __get_chi5(self):
        return get_dihedral("chi5", 9, self._num, self.mdl)

    num = property(__get_num, __set_num, doc="PDB-style residue number")
    atoms = property(__get_atoms, doc="All atoms in this residue")
    alpha = property(__get_alpha, doc="Alpha dihedral angle")
    phi = property(__get_phi, doc="Phi dihedral angle")
    psi = property(__get_psi, doc="Psi dihedral angle")
    omega = property(__get_omega, doc="Omega dihedral angle")
    chi1 = property(__get_chi1, doc="Chi1 dihedral angle")
    chi2 = property(__get_chi2, doc="Chi2 dihedral angle")
    chi3 = property(__get_chi3, doc="Chi3 dihedral angle")
    chi4 = property(__get_chi4, doc="Chi4 dihedral angle")
    chi5 = property(__get_chi5, doc="Chi5 dihedral angle")


class ResidueList(object):
    """A list of L{Residue} objects."""

    def __init__(self, mdl, offset=0, length=None, suffix=""):
        self.mdl = mdl
        self.offset = offset
        self.length = length
        self.suffix = suffix

    def __repr__(self):
        ln = len(self)
        s = "%d residue" % ln
        if ln != 1:
            s += 's'
        return s

    def __str__(self):
        return "<List of " + repr(self) + ">"

    def __len__(self):
        if self.length is not None:
            return self.length
        else:
            return self.mdl.nres

    def get_atom_indices(self):
        (startind, endind) = get_residue_atom_indices(self.mdl, self.offset,
                                                      self.offset + len(self))
        return range(startind + 1, endind + 1), self.mdl

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx, self.mdl._indxres,
                                      (self.offset, self.length, self.suffix))
        if isinstance(ret, int):
            return self.mdl._residue_class(self.mdl, ret + self.offset)
        else:
            return [self[ind] for ind in ret]


class Coordinates(sequence.Sequence):
    """Holds protein coordinates"""

    #: class for atoms in this structure
    _atom_class = Atom
    _residue_class = Residue
    _residue_list_class = ResidueList
    _chain_class = Chain

    def __init__(self):
        sequence.Sequence.__init__(self)

    def get_list_atom_indices(self, objlist, num):
        inds = []
        nat = len(self.atoms)
        # If an object has a get_atom_indices method *and* is iterable, use
        # the method rather than the iterator, since it may be more efficient
        # (and, in the case of selections, may preserve order unlike
        # the iterator)
        if hasattr(objlist, 'get_atom_indices') \
           or not hasattr(objlist, '__iter__'):
            objlist = (objlist,)
        for obj in objlist:
            if hasattr(obj, '__iter__') \
               and not hasattr(obj, 'get_atom_indices'):
                inds.extend(self.get_list_atom_indices(obj, None))
            elif isinstance(obj, int):
                if obj < 1 or obj > nat:
                    raise IndexError(("Invalid atom index %d: should be " +
                                      "between %d and %d") % (obj, 1, nat))
                inds.append(obj)
            else:
                (objinds, mdl) = obj.get_atom_indices()
                if mdl != self:
                    raise TypeError("Incorrect model %s: expected %s"
                                    % (str(mdl), str(self)))
                inds.extend(objinds)
        if num is not None and len(inds) != num:
            raise ValueError("Expecting %d atom indices - got %d"
                             % (num, len(inds)))
        return inds

    def _indxatm(self, offset, length, suffix, indx):
        if isinstance(indx, str):
            func = _modeller.mod_coordinates_find_atom
            newindx = func(self.cdpt, self.seqpt, indx+suffix) - 1 - offset
            if newindx < 0 or (length is not None and newindx >= length):
                raise KeyError("No such atom: %s" % indx)
            return newindx
        raise TypeError("Atom IDs must be numbers or strings")

    def _indxres(self, offset, length, suffix, indx):
        if isinstance(indx, str):
            func = _modeller.mod_coordinates_find_residue
            newindx = func(self.cdpt, self.seqpt, indx+suffix) - 1 - offset
            if newindx < 0 or (length is not None and newindx >= length):
                raise KeyError("No such residue: %s" % indx)
            return newindx
        raise TypeError("Residue IDs must be numbers or strings")

    def residue_range(self, start, end):
        """Return a list of residues, running from start to end inclusive"""
        start = self.residues[start]._num
        end = self.residues[end]._num
        if end < start:
            raise ValueError("End residue is before start residue")
        return self._residue_list_class(self, start, end - start + 1)

    def atom_range(self, start, end):
        """Return a list of atoms, running from start to end inclusive"""
        start = self.atoms[start]._num
        end = self.atoms[end]._num
        if end < start:
            raise ValueError("End atom is before start atom")
        return AtomList(self, start, end - start + 1)

    def point(self, x, y, z):
        """Return a point in the Cartesian space of this model"""
        return Point(self, x, y, z)

    def make_chains(self, file, structure_types='structure',
                    minimal_resolution=99.0, minimal_chain_length=30,
                    max_nonstdres=10, chop_nonstd_termini=True,
                    minimal_stdres=30, alignment_format='PIR'):
        """Write out matching chains to separate files"""
        for chn in self.chains:
            if chn.filter(structure_types, minimal_resolution,
                          minimal_chain_length, max_nonstdres,
                          chop_nonstd_termini, minimal_stdres):
                (atom_file, code) = chn.atom_file_and_code(file)
                log.message('make_chains', "Wrote chain %s.chn" % code)
                chn.write(code + '.chn', atom_file, code,
                          'C; Produced by MODELLER', alignment_format,
                          chop_nonstd_termini)

    def __get_atoms(self):
        return AtomList(self)
    def __get_natm(self):
        return _modeller.mod_coordinates_natm_get(self.cdpt)

    natm = property(__get_natm, doc="Number of atoms in this structure")
    atoms = property(__get_atoms, doc="List of all atoms in this structure")


class AtomList(object):
    """A list of L{Atom} objects."""

    def __init__(self, mdl, offset=0, length=None, suffix=""):
        self.mdl = mdl
        self.offset = offset
        self.length = length
        self.suffix = suffix

    def __repr__(self):
        ln = len(self)
        s = "%d atom" % ln
        if ln != 1:
            s += 's'
        return s

    def __str__(self):
        return "<List of " + repr(self) + ">"

    def __len__(self):
        if self.length is not None:
            return self.length
        else:
            return self.mdl.natm

    def get_atom_indices(self):
        return range(self.offset + 1, self.offset + len(self) + 1), self.mdl

    def __getitem__(self, indx):
        ret = modutil.handle_seq_indx(self, indx, self.mdl._indxatm,
                                      (self.offset, self.length, self.suffix))
        if isinstance(ret, int):
            return self.mdl._atom_class(self.mdl, ret + self.offset)
        else:
            return [self[ind] for ind in ret]

    def __contains__(self, indx):
        try:
            ret = self[indx]
            return True
        except KeyError:
            return False


class AtomIndices:
    """Placeholder class to pass 'raw' atom indices directly to selections"""
    def __init__(self, inds, mdl):
        self.inds = inds
        self.mdl = mdl

    def get_atom_indices(self):
        return self.inds, self.mdl


def get_residue_atom_indices(mdl, start, end):
    iatmr1 = _modeller.mod_coordinates_iatmr1_get(mdl.cdpt)
    startind = _modeller.mod_int1_get(iatmr1, start) - 1
    if end < mdl.nres:
        endind = _modeller.mod_int1_get(iatmr1, end) - 1
    else:
        endind = mdl.natm
    return (startind, endind)


def get_dihedral(type, idihtyp, num, mdl):
    """Get a residue dihedral angle, or None if not defined for this residue"""
    if _modeller.mod_coordinates_has_dihedral(mdl.cdpt, mdl.seqpt,
                                              mdl.env.libs.modpt, idihtyp, num):
        return Dihedral(type, idihtyp, num, mdl)
    else:
        return None

class Dihedral(object):
    """A residue dihedral angle (e.g. alpha, omega, phi, psi, chi1)"""

    def __init__(self, type, idihtyp, num, mdl):
        self.type = type
        self.idihtyp = idihtyp
        self._num = num
        self.mdl = mdl

    def __repr__(self):
        return "%s dihedral" % self.type
    def __str__(self):
        return "<%s>" % repr(self)

    def __get_value(self):
        return _modeller.mod_coordinates_dihedral_get(self.mdl.cdpt,
                                                      self.mdl.seqpt,
                                                      self.mdl.env.libs.modpt,
                                                      self.idihtyp, self._num)
    def __get_dihclass(self):
        dih = self.value
        return _modeller.mod_sequence_dihclass_get(self.mdl.seqpt,
                                                   self.mdl.env.libs.modpt,
                                                   self.idihtyp, dih, self._num)
    def __get_atoms(self):
        ats = _modeller.mod_coordinates_dihatoms_get(self.mdl.cdpt,
                                                     self.mdl.seqpt,
                                                     self.mdl.env.libs.modpt,
                                                     self.idihtyp, self._num)
        return [self.mdl.atoms[i-1] for i in ats]

    value = property(__get_value, doc="Current value, in degrees")
    dihclass = property(__get_dihclass, doc="Current dihedral class")
    atoms = property(__get_atoms, doc="Atoms defining this dihedral angle")
