import _modeller
from optimizer import optimizer

class builtin_optimizer(optimizer):

    def __getstate__(self):
        # Python's pickle won't work with internal Modeller objects,
        # so remove them from the list of objects to be pickled (don't
        # delete the key entirely, as __getstate__ must return a non-empty
        # dict in order for our __setstate__ method to be called).
        d = self.__dict__.copy()
        for key in self.__dict__.keys():
            if key.endswith("_modpt") or key.endswith("__optpt"):
                d[key] = None
        return d

    def _prep_builtin_optimizer(self, atmsel, params, ok_keys, extravars):
        (inds, mdl) = atmsel.get_atom_indices()
        if mdl is None or len(inds) == 0:
            raise ValueError("Selection contains no atoms")
        vars = params.copy()
        self._update_params(vars, ok_keys, extravars)
        for key in ok_keys:
            if key not in vars.keys() \
               and key not in ("edat", "libs", "actions", "schedule_scale"):
                raise ValueError("a value must be given for %s" % key)

        edat = vars.pop('edat', None)
        libs = vars.pop('libs', None)
        schedule_scale = vars.pop('schedule_scale', None)

        if edat is None:
            edat = mdl.env.edat
        if libs is None:
            libs = mdl.env.libs
        if schedule_scale is None:
            schedule_scale = mdl.env.schedule_scale

        sched = _modeller.mod_model_sched_get(mdl.modpt)
        _modeller.mod_schedule_step_set(sched, 1)
        scaln = _modeller.mod_schedule_scaln_get(sched)
        for (num, val) in enumerate(schedule_scale):
            _modeller.mod_float2_set(scaln, num, 0, val)
        return (mdl, libs, edat, inds, vars)

    def _prep_actions(self, params):
        opt = self.optpt
        _modeller.mod_optimizer_actions_free(opt)
        if 'actions' in params:
            actions = params.pop('actions')
            if actions is not None:
                if not hasattr(actions, '__iter__'):
                    actions = (actions,)
                for a in actions:
                    _modeller.mod_optimizer_action_new(opt, a.__call__, a.skip,
                                                       a.first, a.last)

    def __get_step(self):
        return _modeller.mod_optimizer_step_get(self.optpt)
    def __get_funcs(self):
        return _modeller.mod_optimizer_funcs_get(self.optpt)
    def __get_init_e(self):
        return _modeller.mod_optimizer_init_e_get(self.optpt)
    def __get_current_e(self):
        return _modeller.mod_optimizer_current_e_get(self.optpt)
    def __get_shiftavr(self):
        return _modeller.mod_optimizer_shiftavr_get(self.optpt)
    def __get_shiftmax(self):
        return _modeller.mod_optimizer_shiftmax_get(self.optpt)
    step = property(__get_step, doc="Current optimization step")
    funcs = property(__get_funcs, doc="Number of calls to objective function")
    init_e = property(__get_init_e, doc="Initial energy")
    current_e = property(__get_current_e, doc="Current energy")
    shiftavr = property(__get_shiftavr, doc="Average atomic shift (angstroms)")
    shiftmax = property(__get_shiftmax, doc="Maximal atomic shift (angstroms)")
