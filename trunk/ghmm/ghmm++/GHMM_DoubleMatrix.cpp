/* @(#)GHMM_DoubleMatrix.cpp created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "ghmm/matrix.h"
#include "GHMM_DoubleMatrix.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DoubleMatrix::GHMM_DoubleMatrix(int my_rows, int my_cols) {
  c_matrix = matrix_d_alloc(my_rows,my_cols);
  rows = my_rows;
  cols = my_cols;

  if (!c_matrix) {
    fprintf(stderr,"GHMM_DoubleMatrix::GHMM_DoubleMatrix(): could not alloc c_matrix.\n");
    exit(1);
  }
}


GHMM_DoubleMatrix::~GHMM_DoubleMatrix() {
  if (c_matrix)
    matrix_d_free(&c_matrix,rows);
}


const char* GHMM_DoubleMatrix::toString() const {
  return "GHMM_DoubleMatrix";
}


void GHMM_DoubleMatrix::print(FILE *file, char *tab, char *separator, char *ending) {
  matrix_d_print(file,c_matrix,rows,cols,tab,separator,ending);
}
