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


GHMM_DoubleVector::GHMM_DoubleVector(int my_len, double default_value) {
  c_vector = NULL;
  len      = 0;

  resize(my_len,default_value);
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


void GHMM_DoubleVector::resize(int new_len, double default_value) {
  if (new_len != len) {
    c_vector = (double*) realloc(c_vector,new_len * sizeof(double));
    
    for (int i = len; i < new_len; ++i)
      c_vector[i] = default_value;
    
    len = new_len;
  }
}


int GHMM_DoubleVector::normalize() {
  return vector_normalize(c_vector,len);
}


void GHMM_DoubleVector::const_values(double c) {
  vector_const_values(c_vector,len,c);
}


void GHMM_DoubleVector::const_preserve_struct(double c) {
  vector_const_preserve_struct(c_vector,len,c);
}


void GHMM_DoubleVector::random_values() {
  vector_random_values(c_vector,len);
}


void GHMM_DoubleVector::random_preserve_struct() {
  vector_random_preserve_struct(c_vector,len);
}


void GHMM_DoubleVector::print_prec(FILE *file, int width, int prec, char *tab, 
				   char *separator, char *ending) {
  vector_d_print_prec(file,c_vector,len,width,prec,tab,separator,ending);
}
