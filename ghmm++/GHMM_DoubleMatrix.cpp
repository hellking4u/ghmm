/* @(#)GHMM_DoubleMatrix.cpp created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "ghmm/matrix.h"
#include "ghmm/vector.h"
#include "GHMM_DoubleMatrix.h"
#include "GHMM_DoubleVector.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DoubleMatrix::GHMM_DoubleMatrix(int my_rows, int my_cols, double default_value) {
  c_matrix = matrix_d_alloc(my_rows,my_cols);
  rows = my_rows;
  cols = my_cols;

  if (!c_matrix) {
    fprintf(stderr,"GHMM_DoubleMatrix::GHMM_DoubleMatrix(): could not alloc c_matrix.\n");
    exit(1);
  }

  int i;
  int j;
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      c_matrix[i][j] = default_value;
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


GHMM_DoubleVector* GHMM_DoubleMatrix::times_vec(GHMM_DoubleVector* vec) {
  GHMM_DoubleVector* result_vector = new GHMM_DoubleVector(rows);

  vector_mat_times_vec(c_matrix,vec->c_vector,rows,cols,result_vector->c_vector);

  return result_vector;
}
