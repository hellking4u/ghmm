/*******************************************************************************
  author       : Alexander Schliep
  filename     : ghmm/ghmm/rng.c
  created      : ?
  $Id$

__copyright__

*******************************************************************************/

#ifndef RNG_H
#define RNG_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

extern gsl_rng * RNG;

void gsl_rng_init(void);

void gsl_rng_timeseed(gsl_rng * r);




#endif /* RNG_H */

