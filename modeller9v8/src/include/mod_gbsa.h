/** \file mod_gbsa.h           Functions for GB/SA scoring.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_GBSA_H
#define MOD_GBSA_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Create the GB/SA energy term */
void mod_gbsa_create(struct mod_energy_data *edat, int indx, int physical_type,
                     const char *library, int solvation_model, float cutoff);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_GBSA_H */
