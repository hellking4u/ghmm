/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : /zpr/bspk/data/hmm/ngayeva/hmm/src/sreestimate.h
  created      : TIME: 17:18:06     DATE: Mon 15. November 1999
  last-modified: TIME: 13:54:34     DATE: Wed 20. December 2000
------------------------------------------------------------------------------*/
#ifndef SREESTIMATE_H
#define SREESTIMATE_H
#include "smodel.h"

/**@name SHMM-Baum-Welch-Algorithmus */
/*@{ (Doc++-Group: sreestimate) */

/**
  Struktur zur Uebergabe an sreestimate_baum_welch: Modell, Sequenzfeld
  und berechnete Likelihood.. */
struct smosqd_t{
  /** HMM-Modell Pointer */
  smodel *smo;
  /** sequence_t Pointer */
  sequence_d_t *sqd;
  /** berechnete Likelihood */
  double *logp;
  /** Schranke, die beim sreestimate als Abbruch genommen werden kann
      eps < EPS_ITER_BW: EPS_ITER_BW verwenden */
  double eps;
  /** Max-Anzahl Iterationen fuer sreestimate 
      max_iter > MAX_ITER_BW: MAX_ITER_BW verwenden */
  int max_iter;
}; 
typedef struct smosqd_t smosqd_t;

/**
  Baum-Welch Algorithmus für SHMMs.
  Training der Modellparameter nach Baum-Welch bei
  multiplen Eingabesequenzen inklusive Scaling. Die errechneten 
  Modellparameter werden direkt im uebergebenen Modell verändert.
  @return            0/-1 (Fehler) 
  @param smo         Initialmodell
  @param sqd         Trainingssequenzen
  */
/* int sreestimate_baum_welch(smodel *smo, sequence_d_t *sqd); */
int sreestimate_baum_welch(smosqd_t *cs);

/* bereits im creestimate definiert: */
/* Fkt. zum Loesen der Nullstellengleichungen fuer abgeschn. Normaldichte */
/* double pvonmue(double mue, double A, double B); */

#endif /* SREESTIMATE_H */

/*@} (Doc++-Group: sreestimate) */
