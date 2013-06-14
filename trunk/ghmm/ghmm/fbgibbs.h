#ifndef GHMM_FBGIBBS_H
#define GHMM_FBGIBBS_H

#include "model.h"

#ifdef __cplusplus
extern "C" {
#endif


/*Foward-Backwards Gibbs
  Caclulates Q, state sequence
  @param mo: model
  @param O:  observation sequence
  @param len: length of observation sequence
  @param m: iterations
  return state sequence 
  */
void fbgibbstep (int seed, ghmm_dmodel * mo, int *O, int len, double **A, double **B, double *Pi, double **priorA, double **priorB, double *priorPi, int steps);









#ifdef __cplusplus
}
#endif
#endif
