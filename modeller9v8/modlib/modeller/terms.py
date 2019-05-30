"""Classes for handling user-defined energy terms"""

from modeller.util.modlist import LinkList
import _modeller

class energy_term(object):
    """Base class for user-defined energy terms. Override the L{eval} method
       and the L{_physical_type} member in your own subclasses."""

    _physical_type = None
    _cutoff = 0.0

    def _add_term(self, edat, indx):
        """Register this energy term with Modeller"""
        _modeller.mod_energy_term_new(edat, indx, self.eval, self._cutoff,
                                      self._physical_type.get_type())

    def eval(self, mdl, deriv, indats):
        raise NotImplementedError("No Python function for energy term")

    def indices_to_atoms(self, mdl, atom_indices):
        """Converts Modeller-style atom indices into atom objects"""
        return [mdl.atoms[x-1] for x in atom_indices]

class TermList(LinkList):
    """A list of L{energy_term} objects"""
    def __init__(self, edat):
        self.__edat = edat
        self.__list = []
        LinkList.__init__(self)

    def __setstate__(self, d):
        self.__dict__.update(d)
        for (indx, obj) in enumerate(self.__list):
            obj._add_term(self.__edat(), indx)

    def _insfunc(self, indx, obj):
        obj._add_term(self.__edat(), indx)
        self.__list.insert(indx, obj)
    def __len__(self):
        return len(self.__list)
    def _getfunc(self, indx):
        return self.__list[indx]
    def _delfunc(self, indx):
        del self.__list[indx]
        _modeller.mod_energy_term_del(self.__edat(), indx)
