/** \file mod_error.h      Error handling routines
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_ERROR_H
#define MOD_ERROR_H

#include <glib.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Domain for Modeller errors */
#define MOD_ERROR mod_error_quark()

/** Modeller error types */
typedef enum {
  MOD_ERROR_ZERODIV = 0,       /* Divide by zero */
  MOD_ERROR_IO,                /* Input/output error */
  MOD_ERROR_MEMORY,            /* Out of memory */
  MOD_ERROR_EOF,               /* End of file */
  MOD_ERROR_TYPE,              /* Invalid type */
  MOD_ERROR_NOTIMP,            /* Feature not implemented */
  MOD_ERROR_INDEX,             /* Index out of range */
  MOD_ERROR_VALUE,             /* Invalid value */
  MOD_ERROR_FILE_FORMAT,       /* File format error */
  MOD_ERROR_OVERFLOW,          /* Floating point overflow */
  MOD_ERROR_STATISTICS,        /* Statistics problem */
  MOD_ERROR_SEQUENCE_MISMATCH, /* Sequence mismatch between PDB and alignment */
  MOD_ERROR_FAILED             /* Generic error */
} ModellerError;

/** Domain for Modeller errors */
GQuark mod_error_quark(void);

/** Clear the error indicator */
void mod_error_clear(void);

/** Get the current error, or NULL if no error is set. (A copy of the error
    is returned, which must be freed by the user with g_error_free when
    no longer needed.) */
GError *mod_error_get(void);

/** Set the current error, using a copy of the passed-in error (which should
    thus be freed if necessary when no longer needed.) */
void mod_error_set(const GError *err);

#ifdef __cplusplus
}
#endif

#endif  /* MOD_ERROR_H */
