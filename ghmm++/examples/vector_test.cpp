/*******************************************************************************
  authors      : Peter Pipenbacher
  filename     : ghmm++/examples/vector_test.cpp
  $Id$

  __copyright__
*******************************************************************************/

#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


int main() {
  /* Important! initialise rng  */
  GHMM_Toolkit::gsl_rng_init();

  GHMM_DoubleVector v(10,1.0);

  v.print(stdout,"",",","\n");

  printf("normalized vector:\n");

  v.normalize();
  v.print(stdout,"",",","\n");

  printf("resize to 5 and normalize vector:\n");

  v.resize(5);
  v.normalize();
  v.print(stdout,"",",","\n");

  printf("gives all elements the value 2.5:\n");
  v.const_values(2.5); 
  v.print(stdout,"",",","\n");

  printf("gives all elements random values:\n");
  v.random_values();
  v.print(stdout,"",",","\n");

  GHMM_DoubleMatrix matrix(2,5,0.5);
  printf("2x5 matrix M:\n");
  matrix.print(stdout,"",",","\n");

  printf("M x vector:\n");
  GHMM_DoubleVector* v2 = matrix.times_vec(&v);
  v2->print(stdout,"",",","\n");
  SAFE_DELETE(v2);
  
#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return 0;
}
