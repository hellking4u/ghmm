/*------------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : /homes/hmm/bknab/src/sfoba.h
  created      : TIME: 16:47:16     DATE: Mon 15. November 1999
  last-modified: TIME: 16:49:48     DATE: Mon 15. November 1999
------------------------------------------------------------------------------*/

#ifndef SFOBA_H
#define SFOBA_H

#include "smodel.h"

/**@name SHMM Vorwaerts-Rueckwaerts-Algorithmus */
/*@{ (Doc++-Group: sfoba) */

/** Forward-Algorithmus.
  Berechnung von alpha[t][i], Skalierungsfaktoren scale[t] und 
  log( P(O|lambda) ) bei gegebener Sequenz und gegebenem Modell.
  @param smo      vorgebenes Modell
  @param O        Sequenz
  @param T        Laenge der Sequenz
  @param b        Matrix mit bereits berechneten b's (Ausgabe"wahrsch."),
                  darf auch NULL sein!
  @param alpha    alpha[t][i]
  @param scale    Skalierungsfaktoren
  @param log_p    log( P(O|lambda) )
  */
int sfoba_forward(smodel *smo, const double *O, int T, double ***b, 
		  double **alpha, double *scale, double *log_p);

/** 
  Backward-Algorithmus. Berechnung von beta[t][i], bei gegebener Sequenz und
  gegebenem Modell. Skalierungsfaktoren werden uebergeben und nicht neu
  berechnet, da schon aus forward-Alg. bekannt
  @param smo      vorgebenes Modell
  @param O        Sequenz
  @param T        Laenge der Sequenz
  @param b        Matrix mit bereitsberechneten b's (Ausgabe"wahrsch."),
                  darf auch NULL sein!
  @param beta     beta[t][i]
  @param scale    vorgegebene Skalierungsfaktoren
  */
int sfoba_backward(smodel *smo, const double *O, int T, double ***b,
		   double **beta, const double *scale);

/**
  Berechnung von log( P(O|lambda) ). 
  Geschieht durch Aufruf von sfoba\_forward.
  Funktion der Wahl, wenn es ausschliesslich um die Berechnung der
  "Wahrscheinlichkeit" der Sequenz geht, und nicht um alpha[t][i]. 
  @param smo      vorgebenes Modell
  @param O        Sequenz
  @param T        Laenge der Sequenz
  @param log_p    log( P(O|lambda) )
  */
int sfoba_logp(smodel *smo, const double *O, int T, double *log_p);


int sfoba_forwardBS(smodel *smo, const double *O, int T, double ***b,
		    double **alpha, double *scale, double *log_p);

int sfoba_logpBS(smodel *smo, const double *O, int T, double *log_p);


/*@} (Doc++-Group: sfoba) */

#endif
