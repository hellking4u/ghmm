#ifndef CVITERBI_H
#define CVITERBI_H

#include "config.h"
#if __EXPERIMENTAL__ == 3

/**@name CHMM-Viterbi-Algorithmus */
/*@{ (Doc++-Group: cviterbi) */

#include "stdmacro.h"
#include "cmodel.h"

/**
  Viterbi-Algorithmus. Berechnung des Viterbi-Pfades (bester Pfad durch den
  Zustandsraum) und der Viterbi-"Wahrscheinlichkeit" zu einem gegebenen 
  stetigen Modell und einer gegebenen double-Sequenz.
  @return        Viterbi-Pfad (pointer auf int-Vektor)
  @param cmo     Modell
  @param o       double-Sequenz
  @param T       Sequenzlaenge
  @param log_p   "Wahrscheinlichkeit" der Sequenz im Viterbi-Pfad
  */
int *cviterbi(cmodel *cmo, double *o, int T, double *log_p);

/**
  Kontrollrechnung: logP zu geg. Zustandspfad, Sequenz und Modell..
  @return           0/-1 (Fehler)
  @param cmo        Modell
  @param o          Sequenz
  @param T          Sequenzlaenge
  @param state_seq  gegebene Zustandsfolge
  @param log_p      log(P) der Zustandsfolge (Rückgabe)
  */
int cviterbi_logp(cmodel *cmo, double *o, int T, int *state_seq, double *log_p);


#endif /* __EXPERIMENTAL__ == 3 */

#endif

/*@} (Doc++-Group: cviterbi) */


