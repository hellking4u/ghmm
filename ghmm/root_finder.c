#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>


/* struct for function pointer and parameters except 1st one*/
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
};


/*
  this interface is used in sreestimate.c
 */
double zbrent_AB(double (*func)(double, double, double, double),
		 double x1, double x2, double tol,
		 double A, double B, double eps)
{
  /* initialisation of wrapper structure */
  struct parameter_wrapper param={func,A,B,eps};

  gsl_function f={&function_wrapper,&param};
  gsl_interval x={x1,x2};
  gsl_root_fsolver* s;

  double tolerance=tol;
  int success=0;
  double result=0;

  /* initialisation */
  s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent, &f, x);

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
};
