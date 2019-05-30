/** \file mod_user_terms.h     Functions for user-defined energy terms.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_USER_TERM_H
#define MOD_USER_TERM_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Callback function to evaluate a user-defined energy term. Given the atom
    indices, it should return the total energy and, if deriv is True, also
    the first derivatives with respect to the x/y/z atom coordinates. */
typedef int (*cb_teval)(void *data, const struct mod_model *model,
                        gboolean deriv, const int *atind, int n_atind,
                        float *e_term, float *dvx, float *dvy, float *dvz,
                        const struct mod_energy_data *enedata,
                        const struct mod_libraries *libs);

/** Create a new user-defined energy term.
    \param[in] edat energy_data object in which to create the energy term.
    \param[in] indx Index into the array of all energy terms. The new term
                    is inserted before this index.
    \param[in] evalfunc Callback to use for evaluating the energy.
    \param[in] freefunc Callback when the energy term is destroyed.
    \param[in] evaldata Pointer to state data needed by the term.
    \param[in] cutoff If non-zero, make sure the non-bonded list contains
                      atom pairs for this distance.
    \param[in] physical_type Physical restraint type in which to accumulate
                             the energies from this term.
*/
void mod_energy_term_new(struct mod_energy_data *edat, int indx,
                         cb_teval evalfunc, cb_free freefunc, void *evaldata,
                         float cutoff, int physical_type);

/** Remove the energy term at position indx */
void mod_energy_term_del(struct mod_energy_data *edat, int indx);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_USER_TERM_H */
