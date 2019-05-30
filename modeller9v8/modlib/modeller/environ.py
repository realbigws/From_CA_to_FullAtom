from io_data import io_data
from energy_data import energy_data
from libraries import Libraries
import physical
from modeller.util import modutil, logger
from modeller.util.modobject import modobject
import _modeller

__docformat__ = "epytext en"

class environ(modobject):
    """Modeller environment (libraries etc.)"""

    #: factors for physical restraint types in scaling the schedule
    schedule_scale = physical.values(default=1.0)
    #: whether to do default NTER and CTER patching
    patch_default = True

    _rand_seed = None
    _restyp_lib_file = None
    group_restraints = None
    io = None
    edat = None
    libs = None

    def __init__(self, rand_seed=-8123, restyp_lib_file='$(LIB)/restyp.lib',
                 copy=None):
        logger.log.write_header_once()
        self.group_restraints = None
        if copy:
            self.libs = copy.libs
            self.io = io_data(copy=copy.io)
            self.edat = energy_data(copy=copy.edat)
            for member in copy.__dict__:
                if 'environ' not in member and member not in self.__dict__:
                    self.__dict__[member] = copy.__dict__[member]
        else:
            self._rand_seed = rand_seed
            self._restyp_lib_file = restyp_lib_file
            self.libs = Libraries(self._restyp_lib_file, self._rand_seed)
            self.io = io_data()
            self.edat = energy_data()

    def copy(self):
        """Returns a copy of this environment"""
        return environ(copy=self)

    def system(self, command):
        """Run a shell command."""
        return _modeller.mod_system(command, "")

    def dendrogram(self, matrix_file, cluster_cut):
        """Calculate a clustering tree."""
        return _modeller.mod_dendrogram(matrix_file, cluster_cut)

    def principal_components(self, matrix_file, file):
        """Principal components clustering."""
        return _modeller.mod_principal_components(matrix_file, file)

    def make_pssmdb(self, profile_list_file, pssmdb_name, profile_format='TEXT',
                    rr_file='$(LIB)/as1.sim.mat', matrix_offset=0.0,
                    matrix_scaling_factor=0.0069, pssm_weights_type='HH1'):
        """Create a database of PSSMs given a list of profiles"""
        return _modeller.mod_pssmdb_make(self.libs.modpt, profile_list_file,
                                         profile_format, rr_file, matrix_offset,
                                         matrix_scaling_factor, pssmdb_name,
                                         pssm_weights_type)
