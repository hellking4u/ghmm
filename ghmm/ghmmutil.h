/*******************************************************************************
  author       : Benjamin Rich
  filename     : ghmm/ghmm/ghmmutil.h
  created      : TIME: 17:30:55     DATE: Wed 07. July 2004
  $Id$ 

__copyright__

*********************************************************************************/

#ifndef GHMMUTIL_H
#define GHMMUTIL_H

//#ifndef DO_WITH_GSL

/* Non-GSL version declared here */

//#else

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>


/* Types defined in GSL */
typedef gsl_rng GHMM_RNG;

/* Global variables defined in GSL */
#define GHMM_RNG_DEFAULT gsl_rng_default

/* Functions defined in GSL */
#define GHMM_RNG_SET gsl_rng_set
#define GHMM_RNG_UNIFORM gsl_rng_uniform
#define GHMM_RNG_NAME gsl_rng_name

/* Only used in rand.c ??? */
#define GHMM_RNG_ENV_SETUP gsl_rng_env_setup
#define GHMM_RNG_ALLOC gsl_rng_alloc

//#endif

#endif /* GHMMUTIL_H */
