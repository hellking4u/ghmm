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

/** @name rng Initialization for the random number generator */

/**
 */
extern gsl_rng * RNG;


#ifdef __cplusplus
extern "C" {
#endif

/**
 */
void gsl_rng_init(void);

/**
 */
void gsl_rng_timeseed(gsl_rng * r);

#ifdef __cplusplus
}
#endif

#endif /* RNG_H */

