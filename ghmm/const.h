/******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/const.h
  created      : TIME: 13:05:23     DATE: Fri 20. February 1998
  $Id$ 

__copyright__

*******************************************************************************/

#ifndef CONST_H
#define CONST_H

/**@name Konstanten */
/*@{ (Doc++-Group: Konstanten) */


/** 
    Convergence: Halt criterium for Baum-Welch reestimation if the difference
    of log(P) in two consecutive iterations is smaller than (EPS\_ITER\_BW * log(P))..
*/
/* #define EPS_ITER_BW      0.00001 */
#define EPS_ITER_BW      0.0001

/**
  If the absolute difference of two numbers is smaller the EPS_PREC, then the
  numbers are equal. (Instead of using zero )
  */
#define EPS_PREC         1E-8

/**
  Minimum value for U
  */
#define EPS_U            1E-4

/**
  Maximum value for U (Turning point at 100 ?)
  */
#define MAX_U            10000

/**
  Maximum number of iterations in reestimate
  */
#define MAX_ITER_BW      500

/**
  Maximum length of a sequence
  */
#define MAX_SEQ_LEN 30
/* #define MAX_SEQ_LEN      1000000 */

/**
  Maximum number of sequences 
  */
#define MAX_SEQ_NUMBER   1500000

/* in float.h: DBL_EPSILON = 0.000000000000000222044604925031... 
*/

#if 0
/* Now a member of smodel */
#define COS 1
#endif

/**
   A value that is put in for log_p in the calculation of
   the objective function if sfoba_logp returns -1 (error).
*/   
#define PENALTY_LOGP -500.0

/**
   The left limit for the normal density
*/
#define EPS_NDT  0.1

/*@} (Doc++-Group: Constants) */

#endif /* CONST_H */
