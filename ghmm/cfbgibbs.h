#ifndef GHMM_CFBGIBBS_H
#define GHMM_CFBGIBBS_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif

void precompute(int compression, ghmm_dmodel *mo, double ***mats, double**** rmats);

void ghmm_dmodel_cfbgibbstep (ghmm_dmodel *mo, int seed, int *obs, int totalobs, double **pA, double **pB, double *pPi, int* Q, int R);

void ghmm_dmodel_cfbgibbs (ghmm_dmodel *mo, int seed, int *obs, int totalobs, double **pA, double **pB, double *pPi, int* Q, int R, int burnIn);


#ifdef __cplusplus
}
#endif

#endif
