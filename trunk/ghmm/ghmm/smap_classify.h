/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/smap_classify.h
  created      : TIME: 13:40:31     DATE: Wed 12. January 2000
  $Id$

__copyright__

*******************************************************************************/

#ifndef SMAP_CLASSIFY_H
#define SMAP_CLASSIFY_H

#include <ghmm/smodel.h>

#ifdef __cplusplus
extern "C" {
#endif

/**@name smap functions */
/*@{ */

/**
   Maximum A Posteriori Classification Algorithm (MAPCA): 
   given a field of models and one sequence and suppose the sequence
   has been produced by one of these models. This algorithm calculates
   for each model the probability, that the seq. comes from the model.
   This bayesian approach uses a prior for the models. If none is specified
   equal prob. is assumed.
   The maps are copied into "result", which has to be of dimension "smo_number"
   Ref.: A. Kehagias: Bayesian Classification of HMM, Math. Comp. Modelling
   (1995)
   @return number of the model, that fits best to the sequence
   @param smo vector of models
   @param result gives the probability for all the models
   @param smo_number number of models
   @param O sequence
   @param T length of the sequence
 */
int smap_classify(smodel **smo, double *result, int smo_number, 
		   double *O, int T);

/**
   Alternative to MAPCA (smap_classify); calculate p[m] directly using 
   Bayes' theorem, instead of recursive over t.
   p(m | O) = p(O | m) * p(m) / (sum_i p(O | i) * p(i))
   @return number of the model, that fits best to the sequence
   @param smo vector of models
   @param result gives the probability for all the models
   @param smo_number number of models
   @param O sequence
   @param T length of the sequence
 */
int smap_bayes(smodel **smo, double *result, int smo_number, double *O, int T);

#ifdef __cplusplus
}
#endif


/*@} */

#endif /* SMAP_CLASSIFY_H */
