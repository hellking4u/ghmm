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
  Viterbi-Algorithmus. Berechnung des Viterbi-Pfades (bester Pfad durch den
  Zustandsraum) und der Viterbi-Wahrscheinlichkeit zu einem gegebenen Modell
  und einer gegebenen Sequenz.
  @return Viterbi Pfad
  @param mo:    Modell
  @param o:     Sequenz
  @param len:   Sequenzlaenge
  @param log_p: Wahrscheinlichkeit der Sequenz im Viterbi-Pfad
  */
int *viterbi(model *mo, int *o, int len, double *log_p);

/**
  Kontrollrechnung: logP zu gegebenem Zustandspfad (dies muss nicht
  der Viterbi Pfad sein!), gegebener Sequenz und Modell..
  @param mo:        Modell
  @param o:         Sequenz
  @param len:       Sequenzlaenge
  @param state_seq: Zustandsfolge
  @return log P
  */
double viterbi_logp(model *mo, int *o, int len, int *state_seq);

#ifdef __cplusplus
}
#endif


#endif

/*@} (Doc++-Group: viterbi) */
