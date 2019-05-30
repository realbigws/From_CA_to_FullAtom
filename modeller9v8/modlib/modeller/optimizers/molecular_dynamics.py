import _modeller
from builtin_optimizer import builtin_optimizer
import actions
from modeller.util.modutil import handle_seq_indx
from modeller.util.modlist import FixList

class molecular_dynamics(builtin_optimizer):
    """Optimize with molecular dynamics"""
    _modpt = None
    __ok_keys = ( 'edat', 'libs', 'schedule_scale', 'residue_span_range',
                  'max_iterations', 'cap_atom_shift', 'output', 'md_time_step',
                  'init_velocities', 'temperature', 'md_return',
                  'equilibrate', 'guide_factor', 'guide_time', 'friction',
                  'actions' )

    def __init__(self, output='NO_REPORT', cap_atom_shift=0.2,
                 md_time_step=4.0, init_velocities=True, temperature=293.0,
                 md_return='FINAL', equilibrate=999999, guide_factor=0.,
                 guide_time=0., friction=0., residue_span_range=(0, 99999),
                 **vars):
        builtin_optimizer.__init__(self)
        self._modpt = _modeller.mod_md_optimizer_new(self)
        self.__params = {}
        for key in ("output", "cap_atom_shift", "md_time_step",
                    "init_velocities", "temperature", "md_return",
                    "equilibrate", "guide_factor", "guide_time", "friction",
                    "residue_span_range"):
            vars[key] = eval(key)
        self._update_params(self.__params, self.__ok_keys, vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _modeller.mod_md_optimizer_new(self)

    def __del__(self):
        if self._modpt:
            _modeller.mod_md_optimizer_free(self._modpt)

    def __get_optpt(self):
        return _modeller.mod_md_optimizer_opt_get(self._modpt)

    def optimize(self, atmsel, **vars):
        """Optimize the given atom selection"""
        (mdl, libs, edat, inds, vars) = \
            self._prep_builtin_optimizer(atmsel, self.__params, self.__ok_keys,
                                         vars)
        self._prep_actions(vars)
        self.atmsel = atmsel
        func = _modeller.mod_md_optimize
        ret = func(self._modpt, mdl.modpt, edat.modpt, libs.modpt, inds, **vars)
        self.atmsel = None
        return ret

    def trace(self, out):
        """Write out information about the current optimization step"""
        if self.step == 0:
            out.write("# Molecular dynamics simulation at %.2f K\n" \
                      % self.temperature)
            out.write("# %5s %16s %8s %8s %16s %16s\n" \
                      % ('Step', 'Current energy', 'Av shift',
                         'Mx shift', 'Kinetic energy', 'Kinetic temp.'))
        log = "%7d %16.5f %8.4f %8.4f %16.5f %16.5f\n" \
              % (self.step, self.current_e, self.shiftavr, self.shiftmax,
                 self.kinetic_e, self.kinetic_temp)
        out.write(log)

    def __get_min_e(self):
        return _modeller.mod_md_optimizer_min_e_get(self._modpt)
    def __get_max_e(self):
        return _modeller.mod_md_optimizer_max_e_get(self._modpt)
    def __get_kinetic_e(self):
        return _modeller.mod_md_optimizer_kinetic_e_get(self._modpt)
    def __get_kinetic_temp(self):
        return _modeller.mod_md_optimizer_kinetic_temp_get(self._modpt)
    def __get_temperature(self):
        return _modeller.mod_md_optimizer_temperature_get(self._modpt)
    def __set_temperature(self, val):
        _modeller.mod_md_optimizer_temperature_set(self._modpt, val)
    def __get_atoms(self):
        if self.atmsel:
            return MDAtomList(self)
        else:
            raise AttributeError("'atoms' are only available during an " +
                                 "optimization - e.g. from within an 'action'")

    optpt = property(__get_optpt)
    min_e = property(__get_min_e, doc="Minimum energy")
    max_e = property(__get_max_e, doc="Maximal energy")
    kinetic_e = property(__get_kinetic_e, doc="Current kinetic energy")
    kinetic_temp = property(__get_kinetic_temp,
                            doc="Current kinetic temperature")
    temperature = property(__get_temperature, __set_temperature,
                           doc="Set temperature")
    atoms = property(__get_atoms, doc="Per-atom MD information")


class MDAtom(object):
    def __init__(self, mdopt, indx):
        self._mdopt = mdopt
        self._indx = indx
    def __get_guide_force(self):
        return _modeller.mod_md_optimizer_guide_f_get(self._mdopt._modpt,
                                                      self._indx)
    guide_force = property(__get_guide_force, doc="Guiding force")


class MDAtomList(FixList):
    def __init__(self, mdopt):
        self.__mdopt = mdopt
        FixList.__init__(self)

    def _lookup_index(self, indx, require_inrange=True):
        return handle_seq_indx(self, indx, self.__lookup_atom,
                               require_inrange=require_inrange)

    def __len__(self):
        mdl = self.__mdopt.get_selection().get_model()
        return len(mdl.atoms)

    def _getfunc(self, indx):
        return MDAtom(self.__mdopt, indx)
    def __lookup_atom(self, indx):
        mdl = self.__mdopt.get_selection().get_model()
        return mdl.atoms[indx].index - 1
