#include <glib.h>
#include <math.h>
#include "modeller.h"
#include "cuser_feat.h"

/* Get relevant atomic coordinates from Modeller for the feature */
static void get_xyz(const struct mod_model *model, const int atind[2],
                    float x[2], float y[2], float z[2])
{
  int i;
  for (i = 0; i < 2; i++) {
    x[i] = mod_float1_get(&model->cd.x, atind[i] - 1);
    y[i] = mod_float1_get(&model->cd.y, atind[i] - 1);
    z[i] = mod_float1_get(&model->cd.z, atind[i] - 1);
  }
}

/* Evaluate the distance feature value */
static int myfeat_eval(void *data, const struct mod_model *model,
                       const int *atind, int n_atind, float *val)
{
  float x[2], y[2], z[2], dx, dy, dz;
  get_xyz(model, atind, x, y, z);
  dx = x[0] - x[1];
  dy = y[0] - y[1];
  dz = z[0] - z[1];
  *val = sqrt(dx * dx + dy * dy + dz * dz);
  return 0;
}

/* Evaluate the distance feature derivatives */
static int myfeat_deriv(void *data, const struct mod_model *model,
                        const int *atind, int n_atind, float feat,
                        float *dervx, float *dervy, float *dervz)
{
  float x[2], y[2], z[2];
  get_xyz(model, atind, x, y, z);
  dervx[0] = (x[0] - x[1]) / feat;
  dervy[0] = (y[0] - y[1]) / feat;
  dervz[0] = (z[0] - z[1]) / feat;
  dervx[1] = -dervx[0];
  dervy[1] = -dervy[0];
  dervz[1] = -dervz[0];
  return 0;
}

/* Distance is not an 'angle', so return False */
static int myfeat_isangle(void *data, gboolean *angle)
{
  *angle = 0;
  return 0;
}

/* Create the new distance feature, and return its identifier */
int myfeat_create(void)
{
  return mod_user_feature_new(myfeat_eval, NULL, myfeat_deriv, NULL,
                              myfeat_isangle, NULL);
}
