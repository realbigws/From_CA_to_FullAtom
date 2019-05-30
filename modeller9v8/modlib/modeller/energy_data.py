"""Classes to define the type of energy function used."""

import _modeller
import terms
from modeller.error import ModellerError
from modeller.saxsdata import SAXSList
from modeller.util.modobject import modobject

__docformat__ = "epytext en"

class energy_data(modobject):
    """Defines the type of energy function to use"""

    lennard_jones_switch = [6.5, 7.5]
    coulomb_switch = [6.5, 7.5]
    contact_shell = 4.0
    relative_dielectric = 1.0
    radii_factor = 0.82
    update_dynamic = 0.39
    nonbonded_sel_atoms = 1
    covalent_cys = False
    excl_local = [True, True, True, True]
    nlogn_use = 15
    sphere_stdv = 0.05
    dynamic_pairs = False
    dynamic_sphere = True
    dynamic_coulomb = False
    dynamic_lennard = False
    dynamic_modeller = False
    dynamic_access = False
    density = None
    __proxy = None
    __termlist = None
    __saxslist = None
    __simloc = None
    __distsimloc = None

    def __init__(self, copy=None, **kwargs):
        # Use a proxy for modpt here, to avoid circular references between
        # ourselves and the termlist object:
        self.__proxy = PointerProxy(_modeller.mod_energy_data_new,
                                    _modeller.mod_energy_data_free)
        if copy:
            for member in copy.__dict__:
                if 'energy_data' not in member:
                    self.__dict__[member] = copy.__dict__[member]
            if copy.__termlist is not None:
                self.energy_terms.extend(copy.__termlist)
            if copy.__saxslist is not None:
                self.saxsdata.extend(copy.__saxslist)
            if copy.__simloc is not None:
                self.simloc = copy.simloc
            if copy.__distsimloc is not None:
                self.distsimloc = copy.distsimloc
        for key in kwargs:
            if key in dir(energy_data):
                exec("self."+key+"="+str(kwargs[key]))
            else:
                raise KeyError(str(key))

    def __getstate__(self):
        d = modobject.__getstate__(self)
        d.pop('_energy_data__saxslist', None)
        if self.__saxslist is not None and len(self.__saxslist) > 0:
            raise ModellerError("Cannot pickle SAXS data")
        return d

    def __get_energy_terms(self):
        if self.__termlist is None:
            self.__termlist = terms.TermList(self.__proxy)
        return self.__termlist

    def __get_saxsdata(self):
        if self.__saxslist is None:
            self.__saxslist = SAXSList(self.__proxy())
        return self.__saxslist

    def __get_simloc(self):
        return self.__simloc
    def __set_simloc(self, val):
        self.__simloc = val
        _modeller.mod_user_simloc_new(self.__proxy(), val)
    def __get_distsimloc(self):
        return self.__distsimloc
    def __set_distsimloc(self, val):
        self.__distsimloc = val
        _modeller.mod_user_distsimloc_new(self.__proxy(), val)

    def __get_modpt(self):
        modpt = self.__proxy()
        _modeller.mod_energy_data_set(modpt, self.contact_shell,
                                      self.relative_dielectric,
                                      self.radii_factor, self.sphere_stdv,
                                      self.update_dynamic,
                                      self.nonbonded_sel_atoms, self.nlogn_use,
                                      self.covalent_cys, self.dynamic_pairs,
                                      self.dynamic_sphere,
                                      self.dynamic_coulomb,
                                      self.dynamic_lennard,
                                      self.dynamic_modeller,
                                      self.dynamic_access,
                                      self.lennard_jones_switch,
                                      self.coulomb_switch, self.excl_local)
        if self.density:
            _modeller.mod_energy_data_density_set(modpt, self.density._modpt)
        else:
            _modeller.mod_energy_data_density_unset(modpt)
        return modpt

    modpt = property(__get_modpt)
    energy_terms = property(__get_energy_terms)
    saxsdata = property(__get_saxsdata)
    simloc = property(__get_simloc, __set_simloc)
    distsimloc = property(__get_distsimloc, __set_distsimloc)


class PointerProxy(modobject):
    """A proxy for pointers to Fortran objects. Handles their creation and
       deletion and unpickling (pickling is handled by modobject). Also useful
       to avoid circular references."""
    __modpt = None
    __newfunc = __delfunc = None

    def __init__(self, newfunc, delfunc):
        self.__newfunc = newfunc
        self.__delfunc = delfunc
        self.__modpt = newfunc()

    def __setstate__(self, d):
        self.__dict__.update(d)
        self.__modpt = self.__newfunc()

    def __call__(self):
        return self.__modpt

    def __del__(self):
        if self.__modpt:
            self.__delfunc(self.__modpt)
