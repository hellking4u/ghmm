/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/foba.h
  created      : TIME: 00:00:00     DATE: Tue 00. xxx 0000
  $Id$

__copyright__

*******************************************************************************/


#ifndef FOBA_H
#define FOBA_H

#include <ghmm/model.h>


#ifdef __cplusplus
extern "C" {
#endif

/**@name Forward-Backward-Algorithm 
   
   Forward-Backward-Algorithm for multiple integer
   sequences with scaling.
   For reference see:  
   Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285
       
*/

/*@{ (Doc++-Group: foba) */

/** Forward-Algorithm.
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param length: length of sequence
  @param alpha:  alpha[t][i]
  @param scale:  scale factors
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
int foba_forward(model *mo, const int *O, int length, double **alpha, 
		 double *scale, double *log_p);

/** 
  Backward-Algorithm. 
  Calculates beta[t][i] given an integer sequence and a model. Scale factors 
  given as parameter (come from sfoba\_forward).
  @param  mo      model
  @param O          sequence
  @param length   length of sequence
  @param beta     beta[t][i]emissionSequences
  @param scale    scale factors
  @return 0 for success, -1 for error
  */
int foba_backward(model *mo, const int *O, int length, double **beta, 
		  const double *scale);

/**
  Calculation of  log( P(O|lambda) ). 
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed, alpha is allocated with
  stat_matrix_d_alloc
  @param  mo      model
  @param O        sequence
  @param len       length of sequence
  @param log\_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
int foba_logp(model *mo, const int *O, int len, double *log_p);

/*@} (Doc++-Group: foba) */

#ifdef __cplusplus
}
#endif


#endif
