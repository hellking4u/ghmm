/* @(#)GHMM_IntVector.cpp created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "ghmm/vector.h"
#include "GHMM_IntVector.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_IntVector::GHMM_IntVector() {
  c_vector = NULL;
  len      = 0;
}


GHMM_IntVector::GHMM_IntVector(int* my_c_vector, int my_len) {
  c_vector = my_c_vector;
  len      = my_len;
}


GHMM_IntVector::~GHMM_IntVector() {
  if (c_vector)
    free(c_vector);
}


const char* GHMM_IntVector::toString() const {
  return "GHMM_IntVector";
}


void GHMM_IntVector::print(FILE *file, char *tab, char *separator, char *ending) {
  vector_i_print(file,c_vector,len,tab,separator,ending);
}


void GHMM_IntVector::resize(int new_len, int default_value) {
  c_vector = (int*) realloc(c_vector,new_len * sizeof(int));

  for (int i = len; i < new_len; ++i)
    c_vector[i] = default_value;

  len = new_len;
}


