/** \file mod_user_distsimloc.h    Functions for user-defined distance
 *                                 similarity.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_USER_DISTSIMLOC_H
#define MOD_USER_DISTSIMLOC_H

#include "mod_types.h"
#include "mod_user_simloc.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Install a new user-defined distance local similarity function, using the
    passed callback functions */
void mod_user_distsimloc_new(struct mod_energy_data *edat, cb_simloc evalfunc,
                             cb_free freefunc, void *evaldata);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_USER_DISTSIMLOC_H */
