#include <glib.h>
#include <math.h>
#include <stdlib.h>
#include "modeller.h"
#include "cuser_optimizer.h"

/** Main optimization loop */
void mainloop(struct mod_state_optimizer *opt, struct mod_energy_data *edat,
              const struct mod_libraries *libs, float alpha, float minshift,
              float min_ediff, int maxit, int *ierr)
{
  float *state, *dstate, olde;
  int n_state, n_dstate;

  *ierr = 0;

  mod_state_optimizer_state_get(opt, &state, &n_state);
  mod_state_optimizer_energy(opt, edat, state, n_state, &olde, &dstate,
                             &n_dstate, libs, ierr);
  if (*ierr != 0) {
    g_free(state);
    return;
  }

  while (*ierr == 0) {
    int i;
    float newe;

    for (i = 0; i < n_state; i++) {
      state[i] -= alpha * dstate[i];
    }
    g_free(dstate);
    mod_state_optimizer_energy(opt, edat, state, n_state, &newe, &dstate,
                               &n_dstate, libs, ierr);
    if (*ierr != 0) {
      g_free(state);
      return;
    }

    if (fabs(newe - olde) < min_ediff) {
      mod_logout("Finished at step %d due to energy criterion", opt->opt.step);
      break;
    } else if (opt->opt.shiftmax < minshift) {
      mod_logout("Finished at step %d due to shift criterion", opt->opt.step);
      break;
    } else if (opt->opt.step >= maxit) {
      mod_logout("Finished at step %d due to step criterion", opt->opt.step);
      break;
    }
    if (newe < olde) {
      alpha *= 2.0;
    } else {
      alpha /= 2.0;
    }
    olde = newe;
    mod_state_optimizer_next_step(opt, ierr);
  }
  mod_state_optimizer_finish(opt, ierr);
  g_free(state);
  g_free(dstate);
}
