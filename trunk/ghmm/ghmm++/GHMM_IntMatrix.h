/* @(#)GHMM_IntMatrix.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 * __copyright__
 */

#ifndef _GHMM_INTMATRIX_H
#define _GHMM_INTMATRIX_H 1

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_IntMatrix;
class GHMM_IntVector;

/** */
class GHMM_IntMatrix {

 public:

  /** Constructor. */
  GHMM_IntMatrix(int rows, int cols, int default_value = 0);
  /** Destructor. */
  virtual ~GHMM_IntMatrix();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Writes a double matrix (without parenthesis).
     @param file:       output file
     @param tab:        format: leading tabs
     @param separator:  format: separator for columns
     @param ending:     format: end of a row  
  */
  void print(FILE *file, char *tab, char *separator, char *ending);

  /** C style double matrix. */
  int** c_matrix;
  /** Rows of matrix. */
  int rows;
  /** Columns of matrix. */
  int cols;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_DOUBLEMATRIX_H */
