/** \file mod_info.h       General information functions
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_INFO_H
#define MOD_INFO_H

#ifdef __cplusplus
extern "C" {
#endif

/** Get the current date and time */
char *mod_current_time_get(void);

/** Print wall and CPU time since this routine was last called */
void mod_time_mark(void);

/** Return the Modeller header (build and copyright info) */
char *mod_header_get(void);

/** Write the Modeller header (build and copyright info) to the log. */
void mod_header_write(void);

/** Return the short Modeller version string */
char *mod_short_version_get(void);

/** Return the long Modeller version string */
char *mod_long_version_get(void);

/** Called prior to every TOP command */
void mod_top_pre(void);

/** Called after every TOP command */
void mod_top_post(void);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_INFO_H */
