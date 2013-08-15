#ifndef GHMM_CONTINUOUS_FBGIBBS_H
#define GHMM_CONTINUOUS_FBGIBBS_H

#include "smodel.h"
typedef struct normal_hyper{
    //mean ~ N(mue, nu)
    double mue;
    double nu;
    //var ~ gamma(a,b)
    double a;
    double b;
}normal_hyper;

int* ghmm_cmodel_fbgibbs (ghmm_cmodel * mo, ghmm_cseq* seq,
        double **pA, normal_hyper **pB, double *pPi, int burnIn);

#endif
