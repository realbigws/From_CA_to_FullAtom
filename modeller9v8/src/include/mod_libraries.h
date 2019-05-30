/** \file mod_libraries.h      Methods for mod_libraries objects.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_LIBRARIES_H
#define MOD_LIBRARIES_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Serialize a libraries object to a file. Return TRUE on success. */
gboolean mod_libraries_serialize(struct mod_file *fh,
                                 const struct mod_libraries *libs, int *ierr);

/** Unserialize a libraries object from a file. Return TRUE on success. */
gboolean mod_libraries_unserialize(struct mod_file *fh,
                                   struct mod_libraries *libs, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_LIBRARIES_H */
