/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/ghmm/root_finder.h
  created      : DATE: March 2001
  $Id$

__copyright__

*******************************************************************************/

#ifndef ROOT_FINDER_H
#define ROOT_FINDER_H

#ifdef __cplusplus
extern "C" {
#endif
/**
   @name root finder
 */

/*@{ */

/**
   brent root finding algorithm.
   wrapps this functioncall to the gsl 1D-root finding interface
   @author Achim G\"adke
   @param x1 lower bracket value
   @param x2 upper bracket value
   @param tol tolerance for iteration
   @param A 1st extra parameter
   @param B 2nd extra parameter
   @param eps 3rd extra parameter
   @return root
 */
double zbrent_AB(double (*func)(double, double, double, double),
		 double x1, double x2, double tol, double A, double B, 
		 double eps);

/*@} */
#ifdef __cplusplus
}
#endif

#endif
