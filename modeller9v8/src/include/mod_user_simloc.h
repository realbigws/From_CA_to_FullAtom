/** \file mod_user_simloc.h    Functions for user-defined local similarity.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_USER_SIMLOC_H
#define MOD_USER_SIMLOC_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Callback function to return the user-defined local similarity, given the
    indices of the two alignment sequences, and the residue number. */
typedef int (*cb_simloc)(void *data, const struct mod_alignment *alignment,
                         int ires, int is1, int is2, float *val);

/** Install a new user-defined local similarity function, using the passed
    callback functions */
void mod_user_simloc_new(struct mod_energy_data *edat, cb_simloc evalfunc,
                         cb_free freefunc, void *evaldata);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_USER_SIMLOC_H */
