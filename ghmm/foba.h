#ifndef FOBA_H
#define FOBA_H

#include "model.h"

/**@name Vorwaerts-Rueckwaerts-Algorithmus */
/*@{ (Doc++-Group: foba) */

/** Forward-Algorithmus.
  Berechnung von alpha[t][i], von Skalierungsfaktoren
  scale[t] und von log( P(O|lambda) ) bei gegebener Sequenz und gegebenem
  Modell.
  @param mo:     vorgebenes Modell
  @param O:      Sequenz
  @param length: Laenge der Sequenz
  @param alpha:  alpha[t][i]
  @param scale:  Skalierungsfaktoren
  @param log_p:  log( P(O|lambda) )
  */
int foba_forward(model *mo, const int *O, int length, double **alpha, 
		 double *scale, double *log_p);

/** 
  Backward-Algorithmus. Berechnung von beta[t][i], bei gegebener Sequenz und
  gegebenem Modell. Skalierungsfaktoren werden uebergeben und nicht neu
  berechnet, da schon aus forward-Alg. bekannt
  @param mo:     vorgebenes Modell
  @param O:      Sequenz
  @param length: Laenge der Sequenz
  @param beta:   beta[t][i]
  @param scale:  vorgegebene Skalierungsfaktoren
  */
int foba_backward(model *mo, const int *O, int length, double **beta, 
		  const double *scale);

/**
  Berechnung von log( P(O|lambda) ). 
  Geschieht durch Aufruf von foba\_forward.
  Funktion der Wahl, wenn es ausschliesslich um die Berechnung der
  Wahrscheinlichkeit der Sequenz geht, und nicht um alpha[t][i]. 
  @param mo:     vorgebenes Modell
  @param O:      Sequenz
  @param length: Laenge der Sequenz
  @param log_p:  log( P(O|lambda) )
  */
int foba_logp(model *mo, const int *O, int len, double *log_p);

/*@} (Doc++-Group: foba) */

#endif
