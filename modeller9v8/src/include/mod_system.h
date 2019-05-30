/** \file mod_system.h     Low-level Modeller system access routines
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_SYSTEM_H
#define MOD_SYSTEM_H

#ifdef __cplusplus
extern "C" {
#endif

/** Run a shell command */
void mod_system(const char *command, const char *output, int *ierr);

/** Obtain the machine hostname and OS information */
char *mod_uname(void);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_SYSTEM_H */
