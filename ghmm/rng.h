/*******************************************************************************
  author       : Alexander Schliep
  filename     : ghmm/ghmm/rng.c
  created      : ?
  $Id$

__copyright__

*******************************************************************************/

#ifndef RNG_H
#define RNG_H

#include "ghmmutil.h"

/** @name rng Initialization for the random number generator */

/**
 */
extern GHMM_RNG * RNG;


#ifdef __cplusplus
extern "C" {
#endif

/**
 */
void ghmm_rng_init(void);

/**
 */
void ghmm_rng_timeseed(GHMM_RNG * r);

#ifdef __cplusplus
}
#endif

#endif /* RNG_H */

