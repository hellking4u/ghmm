/*------------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : ghmm/ghmm/sfoba.h
  created      : TIME: 16:47:16     DATE: Mon 15. November 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/

#ifndef SFOBA_H
#define SFOBA_H

#include <ghmm/smodel.h>

#ifdef __cplusplus
extern "C" {
#endif

/**@name SHMM Forward-Backward-Algorithm */
/*@{ (Doc++-Group: sfoba) */

/** Forward-Backward-Algorithm for multiple double
    sequences with scaling.
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285
*/
	

/** Forward-Algorithm.
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param T        length of sequence
  @param alpha    alpha[t][i]
  @param scale    scale factors
  @param log_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
int sfoba_forward(smodel *smo, double *O, int T, double ***b, 
		  double **alpha, double *scale, double *log_p);

/** 
  Backward-Algorithm. 
  Calculates beta[t][i] given a double sequence and a model. Scale factors 
  given as parameter (come from sfoba\_forward).
  @param smo      model
  @param O          sequence
  @param T        length of sequence
  @param b        matrix with precalculated output probabilities. May be NULL
  @param beta     beta[t][i]
  @param scale    scale factors
  @return 0 for success, -1 for error
  */
int sfoba_backward(smodel *smo, double *O, int T, double ***b,
		   double **beta, const double *scale);

/**
  Calculation of  log( P(O|lambda) ). 
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed, alpha matrix is allocated with
  stat_matrix_d_alloc
  @param smo      model
  @param O        sequence
  @param T         length of sequence
  @param log_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
int sfoba_logp(smodel *smo, double *O, int T, double *log_p);


#ifdef __cplusplus
}
#endif


/*@} (Doc++-Group: sfoba) */

#endif
