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
#include "rng.h"
#include "math.h"
#include "time.h"

/* The global RNG */
GHMM_RNG * RNG;

#ifndef DO_WITH_GSL

ghmm_rng_state rng_state;
char rng_name[] = "random";

void ghmm_rng_set(GHMM_RNG * r, unsigned long int seed)
{
  srandom(seed);
}
 
double ghmm_rng_uniform (GHMM_RNG * r)
{
  return ((double) random())/(RAND_MAX + 1.0);
}

const char * ghmm_rng_name (GHMM_RNG * r)
{
  return rng_name;
}

void ghmm_rng_init(void)
{
  initstate(1, rng_state, sizeof(ghmm_rng_state));
  RNG = rng_state;
}

#else

void ghmm_rng_init(void)
{
   gsl_rng_env_setup();
   RNG = gsl_rng_alloc(gsl_rng_default);
}

#endif

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




