#ifndef GHMM_FBGIBBS_H
#define GHMM_FBGIBBS_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif
//**uses gsl**

int sample(int seed, double* dist, int N);

void init_priors(double ***pA, double ***pB, double **pPi, ghmm_dmodel *mo);

void update(int seed, ghmm_dmodel* mo, int T, int *states, int* O, double **pA, double **pB, double *pPi);

void updateH(int seed, ghmm_dmodel* mo, int T, int *states, int* O, double **pA, double **pB, double *pPi);

void ghmm_dmodel_fbgibbstep (ghmm_dmodel * mo, int seed, int *O, int len, double **pA, double **pB, double *pPi, int* Q);

void ghmm_dmodel_fbgibbs (ghmm_dmodel * mo, int seed, int *O, int len, double **pA, double **pB, double *pPi, int* Q, int burnIn); 





#ifdef __cplusplus
}
#endif
#endif
