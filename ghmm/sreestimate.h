/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : ghmm/ghmm/sreestimate.h
  created      : TIME: 17:18:06     DATE: Mon 15. November 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/
#ifndef SREESTIMATE_H
#define SREESTIMATE_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/smodel.h>

/**@name SHMM-Baum-Welch-Algorithm */
/*@{ (Doc++-Group: sreestimate) */

/** Baum-Welch-Algorithm for parameter reestimation (training) in
    a continuous (continuous output functions) HMM. Scaled version
    for multiple sequences. Sequences may carry different weights 
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285    
*/

/** @name struct smosqd_t
    structure that combines a continuous model (smo) and an integer
    sequence struct. Is used by sreestimate\_baum\_welch for 
    parameter reestimation.
 */
struct smosqd_t {
  /** pointer of continuous model*/
  smodel *smo;
  /** sequence\_d\__t pointer */
  sequence_d_t *sqd;
  /** calculated log likelihood */
  double* logp;
  /** leave reestimation loop if diff. between successive logp values 
      is smaller than eps */
  double eps;
  /** max. no of iterations */
  int max_iter;
}; 
typedef struct smosqd_t smosqd_t;

typedef struct local_store_t {
  int cos;
  double *pi_num;
  double pi_denom;
  double ***a_num;
  double **a_denom;
  double **c_num; 
  double *c_denom;
  double **mue_num;
  double **u_num;
  double **mue_u_denom; /* mue-denom. = u-denom. for sym. normal density*/
  double **sum_gt_otot; /* for truncated normal density */
  double **sum_gt_logb; /* Control according to Q-function */
} local_store_t;


/**
  Baum-Welch Algorithm for SHMMs.
  Training of model parameter with multiple double sequences (incl. scaling).
  New parameters set directly in hmm (no storage of previous values!). Matrices
  are allocated with stat_matrix_d_alloc.
  @return            0/-1 success/error
  @param cs         initial model and train sequences
  */
int sreestimate_baum_welch(smosqd_t *cs);


int sreestimate_one_step(smodel *smo, local_store_t *r, int seq_number,int *T,  double **O, double *log_p, double *seq_w);

#ifdef __cplusplus
}
#endif


#endif /* SREESTIMATE_H */

/*@} (Doc++-Group: sreestimate) */
