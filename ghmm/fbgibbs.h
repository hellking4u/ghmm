#ifndef GHMM_FBGIBBS_H
#define GHMM_FBGIBBS_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif
//warnining work in progress
//**uses gsl**

void ghmm_dmodel_fbgibbstep (ghmm_dmodel * mo, int seed, int *O, int len, double **pA, double **pB, double *pPi, int* Q);

void ghmm_dmodel_fbgibbs (ghmm_dmodel * mo, int seed, int *O, int len, double **pA, double **pB, double *pPi, int* Q, int burnIn); 





#ifdef __cplusplus
}
#endif
#endif
