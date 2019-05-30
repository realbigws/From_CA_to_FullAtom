import _modeller
import util.modutil as modutil
from modeller.excluded_pair import excluded_pair
from modeller.rigid_body import rigid_body
from modeller.pseudo_atom import pseudo_atom
from modeller.symmetry import symmetry
from modeller.pseudo_atom_list import PseudoAtomList
from modeller.rigid_body_list import RigidBodyList
from modeller.excluded_pair_list import ExcludedPairList
from modeller.nonbonded_pair_list import NonbondedPairList
from modeller.symmetry_list import SymmetryList
import modeller.physical as physical
import modeller.alignment as alignment

class Restraints(object):
    """Holds all restraints which can act on a model"""

    def __init__(self, mdl):
        self.__mdl = mdl

    def __get_modpt(self):
        return _modeller.mod_model_rsr_get(self.__mdl.modpt)

    def __len__(self):
        return _modeller.mod_restraints_ncsr_get(self.modpt)

    def __str__(self):
        return "<List of %d restraints>" % len(self)

    def add(self, *args):
        """Add one or more specified restraints"""
        for rsr in args:
            if isinstance(rsr, excluded_pair):
                raise TypeError("use restraints.excluded_pairs.append() " +
                                "instead")
            elif isinstance(rsr, rigid_body):
                raise TypeError("use restraints.rigid_bodies.append() instead")
            elif isinstance(rsr, pseudo_atom):
                raise TypeError("use restraints.pseudo_atoms.append() instead")
            elif isinstance(rsr, symmetry):
                raise TypeError("use restraints.symmetry.append() instead")
            elif hasattr(rsr, "_add_restraint"):
                return rsr._add_restraint(self, self.__mdl)
            else:
                raise TypeError("must provide restraint objects here")

    def condense(self):
        """Delete unselected or redundant restraints"""
        self.unpick_redundant()
        self.remove_unpicked()

    def remove_unpicked(self):
        """Remove unselected restraints"""
        return _modeller.mod_restraints_remove_unpicked(self.modpt)

    def unpick(self, *atom_ids):
        """Unselect restraints acting on specified atoms"""
        inds = self.__mdl.get_list_atom_indices(atom_ids, None)
        return _modeller.mod_restraints_unpick(self.__mdl.modpt, inds)

    def unpick_redundant(self):
        """Unselect redundant cosine dihedral restraints"""
        return _modeller.mod_restraints_unpick_redundant(self.modpt)

    def unpick_all(self):
        """Unselect all restraints"""
        return _modeller.mod_restraints_unpick_all(self.modpt)

    def clear(self):
        """Delete all restraints"""
        return _modeller.mod_restraints_clear(self.__mdl.modpt)

    def append(self, file):
        """Append restraints from a file, and select them"""
        return _modeller.mod_restraints_read(self.__mdl.modpt, file)

    def read(self, file):
        """Read restraints from a file, and select them"""
        self.clear()
        return self.append(file)

    def write(self, file):
        """Write currently selected restraints to a file"""
        return _modeller.mod_restraints_write(self.__mdl.modpt, file)

    def pick(self, atmsel, residue_span_range=(0, 99999), restraint_sel_atoms=1,
             restraints_filter=physical.values(default=-999)):
        """Select specified restraints"""
        (inds, mdl) = atmsel.get_atom_indices()
        if mdl is not self.__mdl:
            raise ValueError("selection refers to a different model")
        return _modeller.mod_restraints_pick(self.__mdl.modpt, inds,
                                             residue_span_range,
                                             restraint_sel_atoms,
                                             restraints_filter)

    def spline(self, form, feature, group, spline_dx=0.5, spline_range=4.0,
               spline_min_points=5, output='', edat=None):
        """Approximate selected restraints by splines"""
        if edat is None:
            edat = self.__mdl.env.edat
        spline_select = (form.get_type(), feature.get_type(), group.get_type())
        return _modeller.mod_restraints_spline(self.__mdl.modpt, edat.modpt,
                                               self.__mdl.env.libs.modpt,
                                               spline_dx, spline_range,
                                               spline_select,
                                               spline_min_points, output)

    def reindex(self, mdl):
        """Renumber restraints for a new model"""
        return _modeller.mod_restraints_reindex(self.__mdl.modpt, mdl.modpt)

    def make(self, atmsel, restraint_type, spline_on_site,
             residue_span_range=(0, 99999),
             residue_span_sign=True, restraint_sel_atoms=1,
             basis_pdf_weight='LOCAL', basis_relative_weight=0.05,
             intersegment=True, dih_lib_only=False, spline_dx=0.5,
             spline_min_points=5, spline_range=4.0, mnch_lib=1,
             accessibility_type=8, surftyp=1, distngh=6.0, aln=None,
             edat=None, io=None):
        """Calculates and selects new restraints of a specified type"""
        if not hasattr(atmsel, "get_atom_indices"):
            raise TypeError("First argument needs to be an atom selection")
        if not isinstance(restraint_type, str):
            raise TypeError("restraint_type must be a string")
        restyp = restraint_type.upper()

        (inds, mdl) = atmsel.get_atom_indices()
        if mdl is None:
            raise ValueError("selection is empty")
        if mdl is not self.__mdl:
            raise ValueError("selection refers to a different model")

        if edat is None:
            edat = self.__mdl.env.edat
        if io is None:
            io = self.__mdl.env.io
        if aln is None:
            if restyp in \
                ('CHI1_DIHEDRAL', 'CHI2_DIHEDRAL', 'CHI3_DIHEDRAL',
                 'CHI4_DIHEDRAL', 'PHI_DIHEDRAL', 'PSI_DIHEDRAL',
                 'OMEGA_DIHEDRAL', 'PHI-PSI_BINORMAL'):
                raise ValueError("You must provide an alignment for this " +
                                 "restraint type")
            else:
                aln = alignment.alignment(self.__mdl.env)
        func = _modeller.mod_restraints_make
        return func(self.__mdl.modpt, edat.modpt, aln.modpt, io.modpt,
                    self.__mdl.env.libs.modpt, inds, (), (), residue_span_range,
                    restraint_type, restraint_sel_atoms, 0,
                    basis_pdf_weight, 0, 0., basis_relative_weight,
                    spline_on_site, residue_span_sign, intersegment,
                    dih_lib_only, spline_dx, spline_min_points, spline_range,
                    mnch_lib, accessibility_type, (0,0), (0,0,0), surftyp,
                    distngh, 0.0)

    def make_distance(self, atmsel1, atmsel2, aln, spline_on_site,
                      restraint_group, maximal_distance,
                      residue_span_range=(0, 99999), residue_span_sign=True,
                      distance_rsr_model=1, basis_pdf_weight='LOCAL',
                      basis_relative_weight=0.05, spline_dx=0.5,
                      spline_min_points=5, spline_range=4.0,
                      accessibility_type=8, restraint_stdev=(0.1, 1.0),
                      restraint_stdev2=(0., 0., 0.), surftyp=1, distngh=6.0,
                      edat=None, io=None, exclude_distance=0.0):
        """Calculates and selects new distance restraints. Any pair with
           distance less than C{exclude_distance} is also excluded from the
           nonbonded list."""
        if not isinstance(restraint_group, physical.physical_type):
            raise TypeError("restraint_group should be a physical_type object")

        (inds1, mdl) = atmsel1.get_atom_indices()
        if mdl is None:
            raise ValueError("selection is empty")
        elif mdl is not self.__mdl:
            raise ValueError("selection refers to a different model")
        (inds2, mdl) = atmsel2.get_atom_indices()
        if mdl is None:
            raise ValueError("selection is empty")
        elif mdl is not self.__mdl:
            raise ValueError("selection refers to a different model")

        if edat is None:
            edat = self.__mdl.env.edat
        if io is None:
            io = self.__mdl.env.io
        func = _modeller.mod_restraints_make
        return func(self.__mdl.modpt, edat.modpt, aln.modpt, io.modpt,
                    self.__mdl.env.libs.modpt, (), inds1, inds2,
                    residue_span_range, 'DISTANCE', 0,
                    restraint_group.get_type(), basis_pdf_weight,
                    distance_rsr_model, maximal_distance, basis_relative_weight,
                    spline_on_site, residue_span_sign, False, False, spline_dx,
                    spline_min_points, spline_range, False, accessibility_type,
                    restraint_stdev, restraint_stdev2, surftyp, distngh,
                    exclude_distance)

    def __get_pseudo_atoms(self):
        return PseudoAtomList(self.__mdl)
    def __get_rigid_bodies(self):
        return RigidBodyList(self.__mdl)
    def __set_rigid_bodies(self, obj):
        modlist.set_varlist(self.rigid_bodies, obj)
    def __del_rigid_bodies(self):
        modlist.del_varlist(self.rigid_bodies)
    def __get_excluded_pairs(self):
        return ExcludedPairList(self.__mdl)
    def __set_excluded_pairs(self, obj):
        modutil.set_varlist(self.excluded_pairs, obj)
    def __del_excluded_pairs(self):
        modutil.del_varlist(self.excluded_pairs)
    def __get_symmetry(self):
        return SymmetryList(self.__mdl)

    modpt = property(__get_modpt)
    pseudo_atoms = property(__get_pseudo_atoms, doc="Pseudo atoms")
    rigid_bodies = property(__get_rigid_bodies, __set_rigid_bodies,
                            __del_rigid_bodies, doc="Rigid bodies")
    excluded_pairs = property(__get_excluded_pairs, __set_excluded_pairs,
                              __del_excluded_pairs, doc="Excluded pairs")
    nonbonded_pairs = property(lambda x: NonbondedPairList(x.__mdl),
                               doc="Nonbonded pairs")
    symmetry = property(__get_symmetry, doc="Symmetry restraints")
