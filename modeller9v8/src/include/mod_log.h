/** \file mod_log.h        Logging functions.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_LOG_H
#define MOD_LOG_H

#include <stdarg.h>
#include "mod_types.h"
#include "mod_error.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Prototype for callback functions accepted by set_log_function() */
typedef int (*cb_log)(void *logdata, const char *text);

/** Write text directly to the log */
void mod_log_write(const char *str);

/** Write an 'output' message to the log */
void mod_logout(const char *format, ...);

/** Write an optional 'note' message to the log */
void mod_lognote(const char *format, ...);

/** Write an error message to the log (and set the error indicator) */
void mod_logerror(const char *routine, ModellerError type,
                  const char *format, ...);

/** Write a warning message to the log */
void mod_logwarning(const char *routine, const char *format, ...);

/** Get the current log level */
int mod_log_get(int level);

/** Set the current log level */
void mod_log_set(int level, int value);

/** Install a callback function to handle all log output */
void mod_log_function_set(cb_log logfunc, cb_free freefunc, void *logdata);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_LOG_H */
