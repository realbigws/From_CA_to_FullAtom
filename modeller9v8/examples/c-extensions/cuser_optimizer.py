from modeller.optimizers import state_optimizer

# Load the C extension module; this needs to be compiled first - see
# cuser_feat.py for suitable commands.
import _cuser_optimizer

class SteepestDescent(state_optimizer):
    """Very simple steepest descent optimizer, in C"""

    # Add options for our optimizer
    _ok_keys = state_optimizer._ok_keys + ('min_atom_shift', 'min_e_diff',
                                           'step_size', 'max_iterations')

    def __init__(self, step_size=0.0001, min_atom_shift=0.01, min_e_diff=1.0,
                 max_iterations=None, **vars):
        state_optimizer.__init__(self, step_size=step_size,
                                 min_atom_shift=min_atom_shift,
                                 min_e_diff=min_e_diff,
                                 max_iterations=max_iterations, **vars)

    def optimize(self, atmsel, **vars):
        # Do normal optimization startup
        state_optimizer.optimize(self, atmsel, **vars)

        # Get all parameters
        alpha = self.get_parameter('step_size')
        minshift = self.get_parameter('min_atom_shift')
        min_ediff = self.get_parameter('min_e_diff')
        maxit = self.get_parameter('max_iterations')

        (opt, edat, libs) = self.get_modeller_objects()
        ierr = _cuser_optimizer.mainloop(opt, edat, libs, alpha, minshift,
                                         min_ediff, maxit)
        if ierr != 0:
            raise ModellerError
