#include "modeller.h"

void mainloop(struct mod_state_optimizer *opt, struct mod_energy_data *edat,
              const struct mod_libraries *libs, float alpha, float minshift,
              float min_ediff, int maxit, int *ierr);
