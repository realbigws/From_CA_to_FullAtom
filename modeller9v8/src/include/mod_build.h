/** \file mod_build.h      Info functions about the Modeller build
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_BUILD_H
#define MOD_BUILD_H

#ifdef __cplusplus
extern "C" {
#endif

/** Return the date when this binary was built, as a string. */
const char *mod_build_date_get(void);

/** Return the type of this binary, as a string. */
const char *mod_exe_type_get(void);


#ifdef __cplusplus
}
#endif
#endif  /* MOD_BUILD_H */
