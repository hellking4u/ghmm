/*******************************************************************************
  author       : Alexander Schliep
  filename     : ghmm/ghmm/rng.c
  created      : ?
  $Id$

__copyright__

*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include "ghmmutil.h"
#include "rng.h"
#include "math.h"
#include "time.h"

/* The global RNG */
GHMM_RNG * RNG;

void ghmm_rng_init(void)
{
   GHMM_RNG_ENV_SETUP();
   RNG = GHMM_RNG_ALLOC(GHMM_RNG_DEFAULT);
}

void ghmm_rng_timeseed(GHMM_RNG * r)
{
  unsigned long tm; /* Time seed */
  unsigned int timeseed;
  
  timeseed = time(NULL);
  srand(timeseed);
  tm = rand();
  GHMM_RNG_SET(r, tm);
  //printf("# using rng '%s' seed=%ld\n", GHMM_RNG_NAME(r), tm);  
  fflush(stdout);
}


/* End of: rng.c */




