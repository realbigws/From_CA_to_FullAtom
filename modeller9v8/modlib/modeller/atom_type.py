import _modeller

class AtomType(object):
    """An atom type, as defined in the topology and other library files"""

    def __init__(self, mdl, num):
        self.mdl = mdl
        self._num = num

    def __str__(self):
        return "<%s>" % repr(self)
    def __repr__(self):
        return "Atom type %s" % self.name

    def __get_name(self):
        return _modeller.mod_atom_name_from_type(self._num,
                                                 self.mdl.env.libs.modpt)
    def __get_mass(self):
        tpl = self.mdl.env.libs.topology
        amrt = _modeller.mod_topology_amrt_get(tpl._modpt)
        return _modeller.mod_float1_get(amrt, self._num)
    def __set_mass(self, val):
        tpl = self.mdl.env.libs.topology
        amrt = _modeller.mod_topology_amrt_get(tpl._modpt)
        _modeller.mod_float1_set(amrt, self._num, val)
    def __get_element(self):
        tpl = self.mdl.env.libs.topology
        return _modeller.mod_topology_element_get(tpl._modpt, self._num)
    def __set_element(self, val):
        tpl = self.mdl.env.libs.topology
        return _modeller.mod_topology_element_set(tpl._modpt, self._num, val)
    def __get_index(self):
        return self._num + 1

    name = property(__get_name, doc="CHARMM name")
    mass = property(__get_mass, __set_mass, doc="Atomic mass")
    element = property(__get_element, __set_element, doc="Element symbol")
    index = property(__get_index, doc="Internal type index")
