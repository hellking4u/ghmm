/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/ghmm/root_finder.c
  created      : DATE: March 2001
  $Id$

__copyright__

*******************************************************************************/

#include "config.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


/* struct for function pointer and parameters except 1st one */
struct parameter_wrapper{
  double (*func)(double, double, double, double);
  double x2;
  double x3;
  double x4;
} parameter;

/* calls the given function during gsl root solving iteration
   the first parameter variates, all other are kept constant
*/
double function_wrapper(double x, void* p)
{
  struct parameter_wrapper* param=(struct parameter_wrapper*)p;
  return param->func(x,param->x2,param->x3,param->x4);
}


/*
  this interface is used in sreestimate.c
 */
double zbrent_AB(double (*func)(double, double, double, double),
		 double x1, double x2, double tol,
		 double A, double B, double eps)
{
  /* initialisation of wrapper structure */
  struct parameter_wrapper param;
  gsl_function f;
  gsl_interval x;
  gsl_root_fsolver* s;
  double tolerance;
  int success=0;
  double result=0;

  param.func=func; param.x2=A; param.x3=B; param.x4=eps;
  f.function=&function_wrapper; f.params=(void*)&param;
  x.lower=x1; x.upper=x2;
  tolerance=tol;

  /* initialisation */
#ifdef GSL_ROOT_FSLOVER_ALLOC_WITH_ONE_ARG
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  gsl_root_fsolver_set (s,&f,x);
#else
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent, &f, x);
#endif

  /* iteration */
  do {
    success=gsl_root_fsolver_iterate(s);
    if (success==GSL_SUCCESS)
      {
	gsl_interval new_x;
	new_x = gsl_root_fsolver_interval (s);
	success=gsl_root_test_interval(new_x,tolerance,tolerance);
      }
  } while (success==GSL_CONTINUE);

  /* result */
  if (success!=GSL_SUCCESS)
    {
      gsl_error ("solver failed", __FILE__, __LINE__, success) ;
    }
  else
    {
      result=gsl_root_fsolver_root(s);
    }

  /* destruction */
  gsl_root_fsolver_free(s);

  return result;
}




