/* @(#)GHMM_IntMatrix.cpp created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "ghmm/matrix.h"
#include "ghmm/vector.h"
#include "GHMM_IntMatrix.h"
#include "GHMM_IntVector.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_IntMatrix::GHMM_IntMatrix(int my_rows, int my_cols, int default_value) {
  c_matrix = matrix_i_alloc(my_rows,my_cols);
  rows = my_rows;
  cols = my_cols;

  if (!c_matrix) {
    fprintf(stderr,"GHMM_IntMatrix::GHMM_IntMatrix(): could not alloc c_matrix.\n");
    exit(1);
  }

  int i;
  int j;
  for (i = 0; i < rows; ++i)
    for (j = 0; j < cols; ++j)
      c_matrix[i][j] = default_value;
}


GHMM_IntMatrix::~GHMM_IntMatrix() {
  if (c_matrix)
    matrix_i_free(&c_matrix,rows);
}


const char* GHMM_IntMatrix::toString() const {
  return "GHMM_IntMatrix";
}


void GHMM_IntMatrix::print(FILE *file, char *tab, char *separator, char *ending) {
  matrix_i_print(file,c_matrix,rows,cols,tab,separator,ending);
}
