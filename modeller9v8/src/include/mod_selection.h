/** \file mod_selection.h  Functions for handling atom selections.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_SELECTION_H
#define MOD_SELECTION_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Mutate the selected residues of the model to a specified residue
    type; only the irestyp array is changed. You can then use
    SEQUENCE_TO_ALI and WRITE_ALIGNMENT to get the result into a file. */
void mod_selection_mutate(struct mod_model *mdl, struct mod_libraries *libs,
                          const int sel1[], int n_sel1,
                          const char *residue_type, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_SELECTION_H */
