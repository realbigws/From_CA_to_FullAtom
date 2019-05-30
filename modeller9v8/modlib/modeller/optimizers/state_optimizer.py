import _modeller
from builtin_optimizer import builtin_optimizer

class state_optimizer(builtin_optimizer):
    _modpt = None
    _ok_keys = ( 'actions', 'schedule_scale', 'libs', 'edat')

    def __init__(self, **vars):
        builtin_optimizer.__init__(self)
        self._modpt = _modeller.mod_state_optimizer_new(self)
        self._params = {}
        self._update_params(self._params, self._ok_keys, vars)

    def __setstate__(self, d):
        self.__dict__.update(d)
        self._modpt = _modeller.mod_state_optimizer_new(self)

    def __del__(self):
        if self._modpt:
            _modeller.mod_state_optimizer_free(self._modpt)

    def __get_optpt(self):
        return _modeller.mod_state_optimizer_opt_get(self._modpt)

    def optimize(self, atmsel, **vars):
        """Optimize the given atom selection."""
        self.atmsel = atmsel
        (mdl, libs, edat, inds, vars) = \
            self._prep_builtin_optimizer(atmsel, self._params, self._ok_keys,
                                         vars)
        self._prep_actions(vars)
        self.__edat = edat
        self.__libs = libs
        _modeller.mod_state_optimizer_start(self._modpt, mdl.modpt, libs.modpt,
                                            edat.modpt, inds)

    def next_step(self):
        """Proceed to the next optimization step."""
        _modeller.mod_state_optimizer_next_step(self._modpt)

    def finish(self):
        """Do any cleanup at the end of optimization."""
        _modeller.mod_state_optimizer_finish(self._modpt)

    def get_state(self):
        """Converts the selection coordinates into a state vector, and returns
           it."""
        return list(_modeller.mod_state_optimizer_state_get(self._modpt))

    def set_state(self, state):
        """Converts the given state vector back into selection coordinates.
           Also calculates the atom shifts from the old state."""
        _modeller.mod_state_optimizer_state_set(self._modpt, self.__edat.modpt,
                                                state)

    def energy(self, state):
        """Given a state vector, returns the energy and first derivatives."""
        return _modeller.mod_state_optimizer_energy(self._modpt,
                                                    self.__edat.modpt,
                                                    state, self.__libs.modpt)

    def energy_only(self, state):
        """Given a state vector, returns just the energy."""
        return _modeller.mod_state_optimizer_energy_only(self._modpt,
                                                         self.__edat.modpt,
                                                         state,
                                                         self.__libs.modpt)

    def get_parameter(self, key):
        """Get an optimization keyword parameter"""
        return self._params[key]

    def get_modeller_objects(self):
        """Return Modeller optimizer, energy_data, and libraries objects,
           useful when writing the main loop in C"""
        return self._modpt, self.__edat.modpt, self.__libs.modpt

    def trace(self, out):
        """Write out information about the current optimization step"""
        if self.step == 0:
            out.write("# %5s %16s %8s %8s %7s\n" \
                       % ('Step', 'Current energy', 'Av shift',
                         'Mx shift', 'Funcs'))
        log = "%7d %16.5f %8.4f %8.4f %7d\n" \
               % (self.step, self.current_e, self.shiftavr,
                  self.shiftmax, self.funcs)
        out.write(log)

    optpt = property(__get_optpt)
