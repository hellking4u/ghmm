/*******************************************************************************
  author       : Alexander Schliep
  filename     : ghmm/ghmm/rng.c
  created      : ?
  $Id$

__copyright__

*******************************************************************************/

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "rng.h"
#include "math.h"
#include "time.h"

/* The global RNG */
gsl_rng * RNG;

void gsl_rng_init(void)
{
   gsl_rng_env_setup(); /* rng is choosen according to GSL_RNG_TYPE env var */
   RNG = gsl_rng_alloc(gsl_rng_default);
}

void gsl_rng_timeseed(gsl_rng * r)
{
  unsigned long tm; /* Time seed */
  unsigned int timeseed;
  
  timeseed = time(NULL);
  srand(timeseed);
  tm = rand();
  gsl_rng_set(r, tm);
  //printf("# using GSL rng '%s' seed=%ld\n", gsl_rng_name(r), tm);  
  fflush(stdout);
}


/* End of: rng.c */




