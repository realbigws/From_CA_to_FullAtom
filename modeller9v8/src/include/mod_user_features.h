/** \file mod_user_features.h  Functions for user-defined feature types.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_USER_FEAT_H
#define MOD_USER_FEAT_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Callback function to evaluate a user-defined feature. Given the passed
    atom indices, return the feature value 'val'. */
typedef int (*cb_feval)(void *data, const struct mod_model *model,
                        const int *atind, int n_atind, float *val);

/** Callback function to evaluate the derivatives of a user-defined feature.
    Given the passed atom indices and feature value, return the x/y/z first
    derivatives. */
typedef int (*cb_fderiv)(void *data, const struct mod_model *model,
                         const int *atind, int n_atind, float feat,
                         float *dervx, float *dervy, float *dervz);

/** Callback function to identify whether this feature is an 'angle' (and thus
    is automatically wrapped by Modeller to ensure it remains between -pi
    and +pi). */
typedef int (*cb_fangle)(void *data, gboolean *angle);

/** Create a new user-defined feature, using the passed callback functions.
    The new feature type index is returned. */
int mod_user_feature_new(cb_feval evalfunc, void *evaldata,
                         cb_fderiv derivfunc, void *derivdata,
                         cb_fangle anglefunc, void *angledata);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_USER_FEAT_H */
