#ifndef VITERBI_H
#define VITERBI_H

/**@name Viterbi-Algorithmus */
/*@{ (Doc++-Group: viterbi) */

#include "model.h"
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

#endif

/*@} (Doc++-Group: viterbi) */
