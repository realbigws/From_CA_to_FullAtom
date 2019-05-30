from modeller.optimizers import state_optimizer

class SteepestDescent(state_optimizer):
    """Very simple steepest descent optimizer, in Python"""

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

        # Main optimization loop
        state = self.get_state()
        (olde, dstate) = self.energy(state)
        while True:
            for i in range(len(state)):
                state[i] -= alpha * dstate[i]
            (newe, dstate) = self.energy(state)
            if abs(newe - olde) < min_ediff:
                print "Finished at step %d due to energy criterion" % self.step
                break
            elif self.shiftmax < minshift:
                print "Finished at step %d due to shift criterion" % self.step
                break
            elif maxit is not None and self.step >= maxit:
                print "Finished at step %d due to step criterion" % self.step
                break
            if newe < olde:
                alpha *= 2
            else:
                alpha /= 2
            olde = newe
            self.next_step()
        self.finish()
