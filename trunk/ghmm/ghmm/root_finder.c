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
#ifdef HAVE_GSL_INTERVAL /* gsl_interval vanished with version 0.9 */
  gsl_interval x;
#endif
  gsl_root_fsolver* s;
  double tolerance;
  int success=0;
  double result=0;

  param.func=func; param.x2=A; param.x3=B; param.x4=eps;
  f.function=&function_wrapper; f.params=(void*)&param;
  tolerance=tol;

  /* initialisation */
#ifdef HAVE_GSL_INTERVAL
# ifdef GSL_ROOT_FSLOVER_ALLOC_WITH_ONE_ARG
  x.lower=x1;
  x.upper=x2;
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  gsl_root_fsolver_set (s,&f,x);
# else
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent, &f, x);
# endif
#else /* gsl_interval vanished with version 0.9 */
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  gsl_root_fsolver_set (s,&f,x1,x2);
#endif

  /* iteration */
  do {
    success=gsl_root_fsolver_iterate(s);
    if (success==GSL_SUCCESS)
      {
#ifdef HAVE_GSL_INTERVAL
	gsl_interval new_x;
	new_x = gsl_root_fsolver_interval (s);
	success=gsl_root_test_interval(new_x,tolerance,tolerance);
#else /* gsl_interval vanished with version 0.9 */
	double x_up;
	double x_low;
	(void)gsl_root_fsolver_iterate (s);
	x_up=  gsl_root_fsolver_x_upper(s);
	x_low= gsl_root_fsolver_x_lower(s);
	success=gsl_root_test_interval(x_low,x_up,tolerance,tolerance);
#endif
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




