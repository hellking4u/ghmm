#ifndef CVITERBI_H
#define CVITERBI_H

#ifdef __cplusplus
extern "C" {
#endif


#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* #if __EXPERIMENTAL__ == 3  */

/**@name CHMM-Viterbi-Algorithmus */
/*@{ (Doc++-Group: sviterbi) */

#include "smodel.h"

/**
   Viterbi algorithm: calculation of the viterbi path (best possible
   state sequenz for a given sequenz and a given model (smo)). Also 
   calculates logp according to this path.
  @return        Viterbi-path 
  @param smo    model
  @param o       double-sequence
  @param T       sequence length
  @param log_p   log(p) of the sequence using the vitberbi path
  */
int *sviterbi(smodel *smo, double *o, int T, double *log_p);

/* #endif */ /* __EXPERIMENTAL__ == 3 */

#ifdef __cplusplus
}
#endif


#endif

/*@} (Doc++-Group: cviterbi) */


