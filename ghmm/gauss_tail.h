/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/gauss_tail.h
  created      : March 2001 by Achim Gaedke from hmm/src/creestimate.h
  $Id$

__copyright__

*******************************************************************************/


#ifndef GAUSS_TAIL_H
#define GAUSS_TAIL_H


#ifdef __cplusplus
extern "C" {
#endif

/**
   @name some calculations concerning gaussian tail function
 */

/*@{ */
/**
 */
double pmue(double mue, double A, double B, double eps);

/**
 */
double pmue_umin(double mue, double A, double B, double eps);

/** 
    @name Function to find the roots of the truncated normal density function.
*/
double pmue_interpol(double mue, double A, double B, double eps);

#ifdef __cplusplus
}
#endif


/*@} */
#endif /* GAUSS_TAIL_H */
