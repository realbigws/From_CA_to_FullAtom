/** \file modeller.h       Functions for comparative model building
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 *
 *
 * Certain conventions are used for parameters and return types for
 * Modeller functions.
 *
 * Scalar parameters:
 *  - int FOO: the integer variable FOO is input-only (ditto for gboolean or
 *    float types).
 *  - const char *FOO: the null-terminated string FOO is input-only.
 *  - const struct *FOO: the pointed-to object FOO is used as an input, but is
 *    not modified.
 *  - struct *FOO: the pointed-to object FOO is used as an input, and may be
 *    modified on output.
 *  - int *FOO: the pointed-to integer FOO is filled in on output with a return
 *    value, but its initial value is not used (ditto for gboolean or float
 *    types).
 *  - char **FOO: the pointer-to char pointer FOO is filled in on output with a
 *    null-terminated string. This must be freed by the user when no longer
 *    required. The initial value is not used.
 *  - int *ierr: an error indicator from the function. The pointed-to integer
 *    is filled in with 0 on success, or non-zero on error. The initial value
 *    is unused.
 *  - int *ierr_inout: as for ierr, but the initial value is used. This is
 *    used for routines that need to release resources and do other cleanup
 *    even when an error is encountered (e.g. file closure).
 *  - const void *scriptobj: a pointer used to attach arbitrary data to a
 *    Modeller object. This is used by the Python interface to hold a reference
 *    to the Python object corresponding to a given Modeller object, so that
 *    callbacks can be passed correctly.
 *
 * Vector parameters:
 *  - const int FOO[n]: the pointed-to array FOO of n integers is used as input
 *    (ditto for gboolean or float types).
 *  - int FOO[n]: the pointed-to array FOO of n integers is filled in on output 
 *    (ditto for gboolean or float types).
 *  - const int FOO[], int n_FOO: the pointed-to array of n_FOO integers is
 *    used as input (ditto for gboolean or float types).
 *  - int **FOO, int *n_FOO: an array of integers is returned as FOO, and the
 *    size of that array is stored in n_FOO. (The two pointers should point to
 *    an int pointer, and an int, respectively.) The int array should be freed
 *    when no longer needed. Ditto for gboolean and float types.
 *  - const char FOO[n], int FOO_len: the array of n strings in Fortran-style
 *    layout (each string is of fixed length FOO_len, padded with blanks) is
 *    used as input.
 *  - const char FOO[], int FOO_len, int n_FOO: the array of n_FOO strings,
 *    each of fixed length FOO_len, is used as input.
 *
 * Return values:
 *  - const char *: a pointer to an internal null-terminated string is returned.
 *    This should NOT be freed.
 *  - char *: a pointer to a string is returned, which should be freed
 *    when no longer needed.
 */

#ifndef MOD_MODELLER_H
#define MOD_MODELLER_H

#include "mod_types.h"
#include "mod_mdt_type.h"
#include "mod_core.h"
#include "mod_selection.h"
#include "mod_file.h"
#include "mod_libs.h"
#include "mod_log.h"
#include "mod_error.h"
#include "mod_build.h"
#include "mod_info.h"
#include "mod_system.h"
#include "mod_ga341.h"
#include "mod_user_features.h"
#include "mod_user_forms.h"
#include "mod_user_simloc.h"
#include "mod_user_distsimloc.h"
#include "mod_user_terms.h"
#include "mod_optimizer_actions.h"
#include "mod_gbsa.h"
#include "mod_libraries.h"

#endif  /* MOD_MODELLER_H */
