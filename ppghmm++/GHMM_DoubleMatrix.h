/* @(#)GHMM_DoubleMatrix.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_DOUBLEMATRIX_H
#define _GHMM_DOUBLEMATRIX_H 1

#include <ppghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_DoubleMatrix;

/** */
class GHMM_DoubleMatrix {

 public:

  /** Constructor. */
  GHMM_DoubleMatrix(int rows, int cols);
  /** Destructor. */
  virtual ~GHMM_DoubleMatrix();

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
  double** c_matrix;
  /** Rows of matrix. */
  int rows;
  /** Columns of matrix. */
  int cols;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ppghmm++/close_code.h>

#endif /* _GHMM_DOUBLEMATRIX_H */
