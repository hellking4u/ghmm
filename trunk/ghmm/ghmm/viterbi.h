/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/viterbi.h
  created      : TIME: 09:07:57     DATE: Fri 19. December 1997
  $Id$

__copyright__

*******************************************************************************/


#ifndef VITERBI_H
#define VITERBI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/model.h>

/**@name Viterbi-Algorithmus */
/*@{ (Doc++-Group: viterbi) */

/**
  Viterbi algorithm. Calculates the Viterbi path (the optimal path trough
  the model) and the Viterbi probability to a given model and a given 
  sequence.
  @return Viterbi path
  @param mo:    model
  @param o:     sequence
  @param len:   length of the sequence
  @param log_p: probability of the sequence in the Viterbi path
  */
int *viterbi(model *mo, int *o, int len, double *log_p);

/**
  Calculates the logarithmic probability to a given path through the 
  states (does not have to be the Viterbi path), given sequence and
  a model.
  @param mo:        model
  @param o:         sequence
  @param len:       length of the sequence
  @param state_seq: path through the states
  @return log P
  */
double viterbi_logp(model *mo, int *o, int len, int *state_seq);

#ifdef __cplusplus
}
#endif


#endif

/*@} (Doc++-Group: viterbi) */
