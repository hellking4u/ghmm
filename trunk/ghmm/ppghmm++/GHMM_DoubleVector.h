/* @(#)GHMM_DoubleVector.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_DOUBLEVECTOR_H
#define _GHMM_DOUBLEVECTOR_H 1


#include <ppghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_DoubleVector;

/** */
class GHMM_DoubleVector {

 public:

  /** Constructor. */
  GHMM_DoubleVector();
  /** Construct from c type vector with given length.
      This object now owns the vector. Thus it must not
      be freed.*/
  GHMM_DoubleVector(double* my_c_vector, int my_len);
  /** Destructor. */
  virtual ~GHMM_DoubleVector();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Writes this double vector (without parenthesis)
     @param file       output stream
     @param tab        format: leading tabs
     @param separator  format: separator for columns
     @param ending     format: end of a row  
  */
  void print(FILE *file, char *tab, char *separator, char *ending);
  /** Resizes vector. */
  void resize(int new_len);

  /** C type double vector. */
  double* c_vector;
  /** length of C type double vector. */
  int len;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ppghmm++/close_code.h>

#endif /* _GHMM_DOUBLEVECTOR_H */
