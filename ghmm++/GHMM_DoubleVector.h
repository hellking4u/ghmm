/* @(#)GHMM_DoubleVector.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_DOUBLEVECTOR_H
#define _GHMM_DOUBLEVECTOR_H 1


#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_DoubleVector;

/** */
class GHMM_DoubleVector {

 public:

  /** Constructor. */
  GHMM_DoubleVector(int len = 0, double default_value = 0.0);
  /** Construct from c type vector with given length.
      This object now owns the vector. Thus it must not
      be freed.*/
  GHMM_DoubleVector(double* my_c_vector, int len);
  /** Destructor. */
  virtual ~GHMM_DoubleVector();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Gives all elements in a vector a constant value
     @param c    given value for the elements
  */
  void const_values(double c);
  /**
     Gives all elements, not equal zero, in a vector a constant value
     @param c    given value for the elements
  */
  void const_preserve_struct(double c);
  /** 
      Scales the sum of elements in a vector to one.
      @return 0 for success; -1 for error
  */
  int normalize();
  /**
     Writes this double vector (without parenthesis)
     @param file       output stream
     @param tab        format: leading tabs
     @param separator  format: separator for columns
     @param ending     format: end of a row  
  */
  void print(FILE *file, char *tab, char *separator, char *ending);
  /**
     Writes a double vector (without parenthesis) with given number of decimal places
     @param file       output file
     @param width      format: total number of decimal places
     @param prec       format: number of decimal places after the comma
     @param tab        format: leading tabs
     @param separator  format: separator for columns
     @param ending     format: end of a row 
  */
  void print_prec(FILE *file, int width, int prec, char *tab, 
		  char *separator, char *ending);
  /** Gives all elements in a vector random values between 0 and 1. */
  void random_values();
  /** Gives all elements, not equal zero, in a vector random values between 0 and 1. */
  void random_preserve_struct();
  /** Resizes vector. When vector is enlarged, all new elements
      will be initialized with 'default_value'. */
  void resize(int new_len, double default_value = 0.0);

  /** C type double vector. */
  double* c_vector;
  /** length of C type double vector. */
  int len;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_DOUBLEVECTOR_H */
