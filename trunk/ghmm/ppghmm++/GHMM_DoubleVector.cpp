/* @(#)GHMM_DoubleVector.cpp created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "ghmm/vector.h"
#include "GHMM_DoubleVector.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DoubleVector::GHMM_DoubleVector() {
  c_vector = NULL;
  len      = 0;
}


GHMM_DoubleVector::GHMM_DoubleVector(double* my_c_vector, int my_len) {
  c_vector = my_c_vector;
  len      = my_len;
}


GHMM_DoubleVector::~GHMM_DoubleVector() {
  if (c_vector)
    free(c_vector);
}


const char* GHMM_DoubleVector::toString() const {
  return "GHMM_DoubleVector";
}


void GHMM_DoubleVector::print(FILE *file, char *tab, char *separator, char *ending) {
  vector_d_print(file,c_vector,len,tab,separator,ending);
}


void GHMM_DoubleVector::resize(int new_len) {
  c_vector = (double*) realloc(c_vector,new_len);

  for (int i = len; i < new_len; ++i)
    c_vector[i] = 0.0;

  len = new_len;
}
