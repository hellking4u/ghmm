/* @(#)GHMM_IntVector.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_INTVECTOR_H
#define _GHMM_INTVECTOR_H 1


#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_IntVector;

/** */
class GHMM_IntVector {

 public:

  /** Constructor. */
  GHMM_IntVector();
  /** Construct from c type vector with given length.
      This object now owns the vector. Thus it must not
      be freed.*/
  GHMM_IntVector(int* my_c_vector, int my_len);
  /** Destructor. */
  virtual ~GHMM_IntVector();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Writes this integer vector (without parenthesis)
     @param file       output stream
     @param tab        format: leading tabs
     @param separator  format: separator for columns
     @param ending     format: end of a row  
  */
  void print(FILE *file, char *tab, char *separator, char *ending);
  /** Resizes vector. When vector is enlarged, all new elements
      will be initialized with 'default_value'. */
  void resize(int new_len, int default_value = 0);

  /** C type integer vector. */
  int* c_vector;
  /** length of C type integer vector. */
  int len;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_INTVECTOR_H */
