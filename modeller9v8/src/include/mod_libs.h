/** \file mod_libs.h       Access functions for Modeller library data.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_LIBS_H
#define MOD_LIBS_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Get license string */
char *mod_license_info_get(void);

/** Get the directory holding Modeller libraries (usually "modlib") */
const char *mod_libdir_get(void);

/** Set the name of the current job */
void mod_jobname_set(const char *jobname);

/** Get the name of the current job */
const char *mod_jobname_get(void);

/** Get the filenames of the libraries used for GA341 assessment */
void mod_ga341libs_get(const char **surflib, const char **pairlib);

/** Get the directory holding Modeller binaries and TOP scripts (usually
    "bin") */
char *mod_bindir_get(void);

/** Get the filename of the configuration file used by the TOP interpreter
    (usually "modlib/top.ini") */
char *mod_inifil_get(void);

/** Return TRUE iff this Modeller binary is a debugging build */
gboolean mod_debug_get(void);

/** Set license key */
void mod_license_key_set(const char *key);

/** Set installation directory */
void mod_install_dir_set(const char *dir);

/** Set whether this Modeller build is using an embedded Python interpreter */
void mod_embedded_set(gboolean embedded);

/** Get whether this Modeller build is using an embedded Python interpreter */
gboolean mod_embedded_get(void);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_LIBS_H */
