import coordinates
import _modeller

class pseudo_atom(coordinates.Point):
    """Base class for all pseudo and virtual atoms"""

    _builtin_index = None
    numatoms = 0

    def __init__(self, *atoms):
        self._atoms = atoms
        self.__mdl = None

    def __repr__(self):
        return "<Pseudo atom of %s>" % str(self._atoms)

    def _check_atoms(self, atoms):
        if len(atoms) != self.numatoms:
            raise ValueError(("This pseudo atom is defined on %d atoms " +
                              "- got %d") % (self.numatoms, len(atoms)))

    def _get_base_atoms(self, mdl):
        """Return the atom indices this pseudo atom is defined on"""
        atoms = mdl.get_list_atom_indices(self._atoms, None)
        self._check_atoms(atoms)
        return atoms

    def _set_atom_index(self, ind, mdl):
        self.__mdl = mdl
        self.__index = ind

    def get_atom_indices(self):
        self.__check_model()
        return ([self.__index], self.__mdl)

    def __check_model(self):
        if self.__mdl is None:
            raise ValueError("You must first add this pseudo atom to a model "
                             + "with the model.restraints.pseudo_atoms.append()"
                             + " method")

    def get_type(self):
        return self._builtin_index

    def __update_position(self):
        self.__check_model()
        _modeller.mod_pseudo_atom_update(self.__mdl.modpt, self.__index)

    def __get_x(self):
        self.__update_position()
        x = _modeller.mod_coordinates_x_get(self.__mdl.cdpt)
        return _modeller.mod_float1_get(x, self.__index-1)
    def __get_y(self):
        self.__update_position()
        y = _modeller.mod_coordinates_y_get(self.__mdl.cdpt)
        return _modeller.mod_float1_get(y, self.__index-1)
    def __get_z(self):
        self.__update_position()
        z = _modeller.mod_coordinates_z_get(self.__mdl.cdpt)
        return _modeller.mod_float1_get(z, self.__index-1)
    def __get_atoms(self):
        return self._atoms
    x = property(__get_x, doc="x coordinate")
    y = property(__get_y, doc="y coordinate")
    z = property(__get_z, doc="z coordinate")
    atoms = property(__get_atoms, doc="Real atoms defining this atom")


class gravity_center(pseudo_atom):
    """Gravity center of all atoms"""
    _builtin_index = 1

    def _check_atoms(self, atoms):
        if len(atoms) < 1:
            raise ValueError("Must specify at least 1 atom for a " +
                             "gravity center")

    def __repr__(self):
        return "<Gravity center of %s>" % str(self._atoms)



class ch2(pseudo_atom):
    """Pseudo aliphatic proton on a tetrahedral carbon (>CH2)
       not assigned stereospecifically; its position is
       between the two real protons; defined by the central
       C and the other two substituents"""
    _builtin_index = 4
    numatoms = 3

    def __repr__(self):
        return "<ch2 of %s>" % str(self._atoms)


class ch31(pseudo_atom):
    """Pseudo aliphatic proton on a tetrahedral carbon (-CH3),
       defined by the central C and the heavy atom X in X-CH3;
       its position is the average of the three real protons"""
    _builtin_index = 6
    numatoms = 2

    def __repr__(self):
        return "<ch31 of %s>" % str(self._atoms)


class ch32(pseudo_atom):
    """Pseudo aliphatic proton between two unassigned -CH3
       groups; defined by X in CH3 - X - CH3 and the two
       C atoms from the two CH3 groups (Val, Leu!);
       its position is the average of the six real protons"""
    _builtin_index = 7
    numatoms = 3

    def __repr__(self):
        return "<ch32 of %s>" % str(self._atoms)
