/** \file mod_user_forms.h     Functions for user-defined restraint forms.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_USER_FORM_H
#define MOD_USER_FORM_H

#include "mod_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Callback function to evaluate a user-defined restraint form. Given the
    passed features (of iftyp types, and with modal modalities) and the
    parameters pcsr, return the function value. If deriv is True, also return
    the first derivatives with respect to each feature. */
typedef int (*cb_reval)(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, gboolean deriv, float *fderv, float *val);

/** Callback function to evaluate the minimum violation of a user-defined
    restraint form, given the features, types, modalities, and parameters. */
typedef int (*cb_vmin)(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val);

/** Callback function to evaluate the heavy violation of a user-defined
    restraint form, given the features, types, modalities, and parameters. */
typedef int (*cb_vheavy)(void *data, const float *feat, const int *iftyp,
                         const int *modal, int n_feat, const float *pcsr,
                         int n_pcsr, float *val);

/** Callback function to evaluate the relative minimum violation of a
    user-defined restraint form, given the features, types, modalities, and
    parameters. */
typedef int (*cb_rvmin)(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val);

/** Callback function to evaluate the relative heavy violation of a user-defined
    restraint form, given the features, types, modalities, and parameters. */
typedef int (*cb_rvheavy)(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val);

/** Callback function to evaluate the minimum feature mean, given the
    features, types, modalities, and parameters. */
typedef int (*cb_rsrmin)(void *data, const float *feat, const int *iftyp,
                         const int *modal, int n_feat, const float *pcsr,
                         int n_pcsr, float *val);

/** Callback function to evaluate the heavy feature mean, given the
    features, types, modalities, and parameters. */
typedef int (*cb_rsrheavy)(void *data, const float *feat, const int *iftyp,
                           const int *modal, int n_feat, const float *pcsr,
                           int n_pcsr, float *val);

/** Create a new user-defined form, using the passed callback functions.
    The new form index is returned. */
int mod_user_form_new(cb_reval evalfunc, void *evaldata,
                      cb_vmin vminfunc, void *vmindata,
                      cb_vheavy vheavyfunc, void *vheavydata,
                      cb_rvmin rvminfunc, void *rvmindata,
                      cb_rvheavy rvheavyfunc, void *rvheavydata,
                      cb_rsrmin rsrminfunc, void *rsrmindata,
                      cb_rsrheavy rsrheavyfunc, void *rsrheavydata);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_USER_FORM_H */
