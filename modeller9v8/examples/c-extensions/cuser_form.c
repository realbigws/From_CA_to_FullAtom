#include <glib.h>
#include "modeller.h"
#include "cuser_form.h"

/* RT at 297.15K, in kcal/mol */
static const float RT = 0.5900991;

/* Decode parameters from Modeller parameter array */
static void get_param(const float pcsr[2], float *mean, float *stdev)
{
  *mean = pcsr[0];
  *stdev = pcsr[1];
}

/* Evaluate the harmonic form */
static int myform_eval(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, gboolean deriv, float *fderv, float *val)
{
  float mean, stdev, stdev2, delta;
  int ierr;
  get_param(pcsr, &mean, &stdev);
  stdev2 = stdev * stdev;
  delta = mod_feature_delta(feat[0], mean, iftyp[0], &ierr);
  if (ierr != 0) {
    return 1;
  }

  *val = RT * 0.5 * delta * delta / stdev2;
  if (deriv) {
    fderv[0] = RT * delta / stdev2;
  }
  return 0;
}

/* Evaluate the minimum violation of the harmonic form */
static int myform_vmin(void *data, const float *feat, const int *iftyp,
                       const int *modal, int n_feat, const float *pcsr,
                       int n_pcsr, float *val)
{
  float mean, stdev;
  int ierr;
  get_param(pcsr, &mean, &stdev);
  *val = mod_feature_delta(feat[0], mean, iftyp[0], &ierr);
  return ierr;
}

/* Evaluate the relative minimum violation of the harmonic form */
static int myform_rvmin(void *data, const float *feat, const int *iftyp,
                        const int *modal, int n_feat, const float *pcsr,
                        int n_pcsr, float *val)
{
  float mean, stdev;
  int ierr;
  get_param(pcsr, &mean, &stdev);
  *val = mod_feature_delta(feat[0], mean, iftyp[0], &ierr) / stdev;
  return ierr;
}

/* Evaluate the minimum feature mean of the harmonic form */
static int myform_minmean(void *data, const float *feat, const int *iftyp,
                          const int *modal, int n_feat, const float *pcsr,
                          int n_pcsr, float *val)
{
  float mean, stdev;
  get_param(pcsr, &mean, &stdev);
  *val = mean;
  return 0;
}

/* Create the new harmonic form, and return its identifier */
int myform_create(void)
{
  /* Only one minimum, so reuse the 'min' functions for the 'heavy' mean */
  return mod_user_form_new(myform_eval, NULL, myform_vmin, NULL, myform_vmin,
                           NULL, myform_rvmin, NULL, myform_rvmin, NULL,
                           myform_minmean, NULL, myform_minmean, NULL);
}
