/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/reestimate.h
  created      : TIME: 12:39:14     DATE: Wed 18. February 1998
  $Id$

__copyright__

*******************************************************************************/


#ifndef REESTIMATE_H
#define REESTIMATE_H
#include <ghmm/sequence.h>
#include <ghmm/model.h>

#ifdef __cplusplus
extern "C" {
#endif

/**@name Baum-Welch-Algorithmus */
/*@{ (Doc++-Group: reestimate) */

/** Baum-Welch-Algorithm for parameter reestimation (training) in
    a discrete (discrete output functions) HMM. Scaled version
    for multiple sequences.  
    New parameters set directly in hmm (no storage of previous values!).
    For reference see:  
    Rabiner, L.R.: "`A Tutorial on Hidden {Markov} Models and Selected
                Applications in Speech Recognition"', Proceedings of the IEEE,
	77, no 2, 1989, pp 257--285    
  @return            0/-1 success/error
  @param mo          initial model
  @param sq          training sequences
  */

int reestimate_baum_welch(model *mo, sequence_t *sq);

#ifdef __cplusplus
}
#endif


#endif /* REESTIMATE_H */

/*@} (Doc++-Group: reestimate) */
