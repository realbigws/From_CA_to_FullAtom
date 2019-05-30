/** \file mod_optimizer_actions.h  Functions for instrumenting and controlling
 *                                 optimizations.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_OPT_ACT_H
#define MOD_OPT_ACT_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Prototype of an action function to be called during optimization */
typedef int (*cb_action)(void *actiondata, void *optobj);

/** Create a new action */
void mod_optimizer_action_new(struct mod_optimizer *opt, cb_action actionfunc,
                              cb_free freefunc, void *actiondata, int skip,
                              gboolean first, gboolean last);

/** Remove all of this optimizer's actions */
void mod_optimizer_actions_free(struct mod_optimizer *opt);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_OPT_ACT_H */
