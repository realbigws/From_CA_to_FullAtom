from topology import Topology
from parameters import Parameters
from modeller.util.modobject import modobject
from modeller import modfile
import _modeller

__docformat__ = "epytext en"

class Libraries(modobject):
    """All library information, such as topology and parameter files, and
       random number state."""

    __modpt = None
    __restyp_lib_file = None
    __rand_seed = None

    def __init__(self, restyp_lib, rand_seed):
        self.__modpt = _modeller.mod_libraries_new(self)
        self.__restyp_lib_file = restyp_lib
        self.__read_libs()
        self.rand_seed = rand_seed

    def __read_libs(self):
        _modeller.mod_libraries_read_libs(self.modpt, self.__restyp_lib_file)

    def __getstate__(self):
        from cStringIO import StringIO
        d = modobject.__getstate__(self)
        s = StringIO()
        self._serialize(s)
        d['__serialize'] = s.getvalue()
        return d

    def __setstate__(self, d):
        from cStringIO import StringIO
        s = StringIO(d.pop('__serialize'))
        self.__dict__.update(d)
        self.__modpt = _modeller.mod_libraries_new(self)
        self.__read_libs()
        self._unserialize(s)

    def __del__(self):
        if self.__modpt:
            _modeller.mod_libraries_free(self.__modpt)

    def random_number(self):
        """Get the next random number between 0 and 1 from Modeller's internal
           random number generator."""
        return _modeller.mod_random_number(self.modpt)

    def random_perturb(self, n):
        """Perturb the internal random number generator by some amount, n."""
        return _modeller.mod_random_perturb(self.modpt, n)

    def _serialize(self, fh):
        """Write out to a file handle"""
        fh = modfile._check_filehandle(fh)
        _modeller.mod_libraries_serialize(fh.file_pointer, self.modpt)

    def _unserialize(self, fh):
        """Read in from a file handle"""
        fh = modfile._check_filehandle(fh)
        _modeller.mod_libraries_unserialize(fh.file_pointer, self.modpt)

    def __get_modpt(self):
        return self.__modpt
    def __get_topology(self):
        return Topology(self)
    def __get_parameters(self):
        return Parameters(self)
    def __get_rand_seed(self):
        return _modeller.mod_libraries_rand_seed_get(self.modpt)
    def __set_rand_seed(self, val):
        _modeller.mod_libraries_rand_seed_set(self.modpt, val)

    modpt = property(__get_modpt)
    topology = property(__get_topology, doc="Topology information")
    parameters = property(__get_parameters, doc="Parameter information")
    rand_seed = property(__get_rand_seed, __set_rand_seed, doc="Random seed")
