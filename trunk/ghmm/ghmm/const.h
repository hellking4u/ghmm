/******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/const.h
  created      : TIME: 13:05:23     DATE: Fri 20. February 1998
  $Id$ 

__copyright__

*******************************************************************************/

#ifndef CONST_H
#define CONST_H

/**@name Konstanten */
/*@{ (Doc++-Group: Konstanten) */

/**
  Konvergenz: Abbruch von Baum-Welch-Reeestimate, falls Differenz 
  aufeinanderfolgender Iterationen log(P) kleiner als (EPS\_ITER\_BW * log(P))..
  */
/* #define EPS_ITER_BW      0.00001 */
#define EPS_ITER_BW      0.0001

/**
  Absolute Differenz zweier Zahlen kleiner EPS\_PREC, dann: Zahlen sind gleich 
  (anstatt Abfrage auf null)
  */
#define EPS_PREC         1E-8

/**
  Gesetztes Minimum fuer U
  */
#define EPS_U            1E-4

/**
  Gesetztes Maximum fuer U (Wendepunkt bei 100 ?)
  */
#define MAX_U            10000

/**
  Maximale Zahl von Iterationen in reestimate
  */
#define MAX_ITER_BW      500

/**
  Maximale Laenge einer Sequenz
  */
#define MAX_SEQ_LEN 30
/* #define MAX_SEQ_LEN      1000000 */

/**
  Maximale Zahl von Sequenzen 
  */
#define MAX_SEQ_NUMBER   1500000

/* in float.h: DBL_EPSILON = 0.000000000000000222044604925031... 
*/

/**
   Anzahl von O-Sum-Klassen.
*/

/* #define COS 6 */

#define COS 1

/** Wert der fuer log_p zur Berechnung von Zielfunktionen
    eingesetzt wird, wenn sfoba_logp den Wert -1
    (Fehler) zurueckliefert. Der Wert von sfoba gelieferte 
    Wert - DBL_MAX ist dafuer unbrauchbar.
*/
    
#define PENALTY_LOGP -500.0
/**
   typedef density\_t fuer cmodel u. smodel.
*/
typedef enum {
  normal, 
  normal_pos, 
  normal_approx,
  density_number
} density_t;

/* linke Grenze der Tilgungssymbole */
/* #define MIN_TILG 5000 */
#define MIN_TILG 0

/**
  Linke Grenze fuer gestutzte NormalDichte
  */
#define EPS_NDT  0.1

/*@} (Doc++-Group: Konstanten) */

#endif /* CONST_H */
