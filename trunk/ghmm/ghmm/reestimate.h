/*******************************************************************************
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/reestimate.h
  created      : TIME: 12:39:14     DATE: Wed 18. February 1998
  last-modified: TIME: 19:08:11     DATE: Wed 25. November 1998
*******************************************************************************/
/* $Id$ */

#ifndef REESTIMATE_H
#define REESTIMATE_H
#include "model.h"

/**@name Baum-Welch-Algorithmus */
/*@{ (Doc++-Group: reestimate) */

/**
  Baum-Welch Algorithmus. Training der Modellparameter nach Baum-Welch bei
  multiplen Eingabesequenzen inklusive Scaling. NEU: Die errechneten 
  Modellparameter werden direkt im uebergebenen Modell veraendert.
  @return 0          alles ok., -1: es ist ein Fehler aufgetreten 
  @param mo          Initialmodell
  @param sq          Trainingssequenzen
  */
 
int reestimate_baum_welch(model *mo, sequence_t *sq);

#endif /* REESTIMATE_H */

/*@} (Doc++-Group: reestimate) */
