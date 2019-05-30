"""All physical restraint types, used for creating restraints, and in the
   schedule."""

class physical_type(object):
    """Physical restraint type"""
    _num = 0
    _types = []
    _shortname = ""
    def __init__(self, doc, shortname):
        physical_type._num += 1
        self._num = physical_type._num
        self.__doc__ = doc
        self._shortname = shortname
        physical_type._types.append(self)

    def get_type(self):
        """Get the numeric type of this physical restraint type"""
        return self._num

    def __str__(self):
        return '<physical restraint type: %s>' % self.__doc__

    def __repr__(self):
        return '<%s>' % self.__doc__
    modpt = property(get_type)


class values(object):
    """Scaling or cutoff values for all physical restraint types"""

    def __init__(self, default=1.0, **keys):
        self._default = default
        self._dict = {}
        for (term,val) in keys.items():
            term = eval("%s" % term)
            self[term] = val

    def set_from_list(self, lval):
        """Set the state from a Modeller internal list of values"""
        for (num, val) in enumerate(lval):
            typ = physical_type._types[num]
            self[typ] = val

    def __len__(self):
        return len(physical_type._types)

    def __repr__(self):
        if self._default is not None:
            terms = ["default=%f" % self._default]
        else:
            terms = []
        for term in physical_type._types:
            num = term.get_type()
            if num in self._dict:
                terms.append("%s=%f" % (term._shortname, self._dict[num]))
        return "physical.values(" + ", ".join(terms) + ")"

    def keys(self):
        keys = []
        for term in physical_type._types:
            if term.get_type() in self._dict:
                keys.append(term)
        return keys

    def __mul__(self, other):
        obj = values()
        if not isinstance(other, values):
            raise TypeError("Must use physical.values objects here")
        obj._default = self._default * other._default
        for x in self._dict:
            if x in other._dict:
                oval = other._dict[x]
            else:
                oval = other._default
            obj._dict[x] = self._dict[x] * oval
        for x in other._dict:
            if x not in self._dict:
                obj._dict[x] = self._default * other._dict[x]
        return obj

    def _typecheck(self, item):
        if not isinstance(item, physical_type):
            raise TypeError("keys must be physical_type objects or 'default'")

    def __contains__(self, item):
        if isinstance(item, str) and item == 'default':
            return True
        self._typecheck(item)
        return item in physical_type._types

    def __getitem__(self, key):
        # Allow read-only access as a list, for the convenience of old code
        # which expects a vector of numbers
        if isinstance(key, int):
            return self[physical_type._types[key]]
        if isinstance(key, str) and key == 'default':
            return self._default
        self._typecheck(key)
        try:
            return self._dict[key.get_type()]
        except KeyError:
            return self._default

    def __setitem__(self, key, value):
        if isinstance(key, str) and key == 'default':
            self._default = value
        self._typecheck(key)
        self._dict[key.get_type()] = value

def from_list(lval):
    """Initialize an object from a Modeller internal list of values"""
    obj = values(default=None)
    obj.set_from_list(lval)
    return obj

bond = physical_type("Bond length potential", "bond")
angle = physical_type("Bond angle potential", "angle")
dihedral = physical_type("Stereochemical cosine torsion potential", "dihedral")
improper = physical_type("Stereochemical improper torsion potential",
                         "improper")
soft_sphere = physical_type("Soft-sphere overlap restraints", "soft_sphere")
lennard_jones = physical_type("Lennard-Jones 6-12 potential", "lennard_jones")
coulomb = physical_type("Coulomb point-point electrostatic potential",
                        "coulomb")
h_bond = physical_type("H-bonding potential", "h_bond")
ca_distance = physical_type("Distance restraints 1 (CA-CA)", "ca_distance")
n_o_distance = physical_type("Distance restraints 2 (N-O)", "n_o_distance")
phi_dihedral = physical_type("Mainchain Phi dihedral restraints",
                             "phi_dihedral")
psi_dihedral = physical_type("Mainchain Psi dihedral restraints",
                             "psi_dihedral")
omega_dihedral = physical_type("Mainchain Omega dihedral restraints",
                               "omega_dihedral")
chi1_dihedral = physical_type("Sidechain Chi_1 dihedral restraints",
                              "chi1_dihedral")
chi2_dihedral = physical_type("Sidechain Chi_2 dihedral restraints",
                              "chi2_dihedral")
chi3_dihedral = physical_type("Sidechain Chi_3 dihedral restraints",
                              "chi3_dihedral")
chi4_dihedral = physical_type("Sidechain Chi_4 dihedral restraints",
                              "chi4_dihedral")
disulfide_distance = physical_type("Disulfide distance restraints",
                                   "disulfide_distance")
disulfide_angle = physical_type("Disulfide angle restraints", "disulfide_angle")
disulfide_dihedral = physical_type("Disulfide dihedral angle restraints",
                                   "disulfide_dihedral")
lower_distance = physical_type("Lower bound distance restraints",
                               "lower_distance")
upper_distance = physical_type("Upper bound distance restraints",
                               "upper_distance")
sd_mn_distance = physical_type("Distance restraints 3 (SDCH-MNCH)",
                               "sd_mn_distance")
chi5_dihedral = physical_type("Sidechain Chi_5 dihedral restraints",
                              "chi5_dihedral")
phi_psi_dihedral = physical_type("Phi/Psi pair of dihedral restraints",
                                 "phi_psi_dihedral")
sd_sd_distance = physical_type("Distance restraints 4 (SDCH-SDCH)",
                               "sd_sd_distance")
xy_distance = physical_type("Distance restraints 5 (X-Y)", "xy_distance")
nmr_distance = physical_type("NMR distance restraints 6 (X-Y)", "nmr_distance")
nmr_distance2 = physical_type("NMR distance restraints 7 (X-Y)",
                              "nmr_distance2")
min_distance = physical_type("Minimal distance restraints", "min_distance")
nonbond_spline = physical_type("Non-bonded spline restraints", "nonbond_spline")
accessibility = physical_type("Atomic accessibility restraints",
                              "accessibility")
density = physical_type("Atomic density restraints", "density")
absposition = physical_type("Absolute position restraints", "absposition")
dihedral_diff = physical_type("Dihedral angle difference restraints",
                              "dihedral_diff")
gbsa = physical_type("GBSA implicit solvent potential", "gbsa")
em_density = physical_type("EM density fitting potential", "em_density")
saxs = physical_type("SAXS restraints", "saxs")
symmetry = physical_type("Symmetry restraints", "symmetry")
