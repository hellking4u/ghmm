/*******************************************************************************
  author       : Alexander Schliep
  filename     : /homes/hmm/wichern/hmm/src/rng.c
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

/* No timeseed in GSL */
void gsl_rng_timeseed(gsl_rng * r)
{
  unsigned long tm; /* Time seed */
  
  /*  tm = time(NULL); */
  tm = rand();
  gsl_rng_set(r, tm);
  /* Printing commented out 04.02.01 by Disa */
  /* printf("# using GSL rng '%s' seed=%ld\n", gsl_rng_name(r), tm);  */
  fflush(stdout);
}


/* End of: rng.c */




