/** \file mod_ga341.h      The GA341 model assessment method.
 *
 *             Part of MODELLER, Copyright(c) 1989-2010 Andrej Sali
 */

#ifndef MOD_GA341_H
#define MOD_GA341_H
#ifdef __cplusplus
extern "C" {
#endif


/** Assess the given model with the GA341 method, and return the score and
    its components. */
void mod_assess_ga341(const struct mod_model *mdl,
                      const struct mod_libraries *libs,
                      float *score, float *compactness, float *e_nat_pair,
                      float *e_nat_surf, float *e_nat_comb, float *zsco_pair,
                      float *zsco_surf, float *zsco_comb, int *ierr);

/** Assess the given model with the specified assessment method (only "GA341"
    is accepted), and return the score.
    \note Deprecated: use assess_ga341() instead. */
void mod_assess_model(const struct mod_model *mdl,
                      const struct mod_libraries *libs,
                      const char *assess_method, float *molpdf, int *ierr);

#ifdef __cplusplus
}
#endif
#endif  /* MOD_GA341_H */
