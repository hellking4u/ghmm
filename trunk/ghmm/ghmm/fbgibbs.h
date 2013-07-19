#ifndef GHMM_FBGIBBS_H
#define GHMM_FBGIBBS_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif
//**uses gsl**

int sample(int seed, double* dist, int N);

void init_priors(ghmm_dmodel *mo, double ***pA, double ***pB, double **pPi);

void update(ghmm_dmodel* mo, int T, int *states, int* O, double **pA, double **pB, double *pPi);

void updateH(ghmm_dmodel* mo, int T, int *states, int* O, double **pA, double **pB, double *pPi);

int* ghmm_dmodel_fbgibbs (ghmm_dmodel * mo, int seed, int *O, int len, double **pA, double **pB, double *pPi, int burnIn); 





#ifdef __cplusplus
}
#endif
#endif
