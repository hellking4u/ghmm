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


void GHMM_DoubleMatrix::print_prec(FILE *file, int width, int prec, char *tab, 
				   char *separator, char *ending) {
  matrix_d_print_prec(file,c_matrix,rows,cols,width,prec,tab,separator,ending);
}


int GHMM_DoubleMatrix::notzero_columns() {
  return matrix_d_notzero_columns(c_matrix,rows,cols);
}


int GHMM_DoubleMatrix::notzero_rows() {
  return matrix_d_notzero_rows(c_matrix,cols,rows);
}


int GHMM_DoubleMatrix::normalize() {
  return matrix_d_normalize(c_matrix,rows,cols);
}


void GHMM_DoubleMatrix::random_values(double min, double max) {
  matrix_d_random_values(c_matrix,rows,cols,min,max); 
}


void GHMM_DoubleMatrix::random_const_values(double min, double max, double c) {
  matrix_d_random_const_values(c_matrix,rows,cols,min,max,c);
}


void GHMM_DoubleMatrix::const_values(double c) {
  matrix_d_const_values(c_matrix,rows,cols,c); 
}


void GHMM_DoubleMatrix::left_right_strict() {
  matrix_d_left_right_strict(c_matrix,rows,cols); 
}


void GHMM_DoubleMatrix::random_left_right(double **matrix, int rows, int cols) {
  matrix_d_random_left_right(c_matrix,rows,cols); 
}


void GHMM_DoubleMatrix::const_preserve_struct(double c) {
  matrix_d_const_preserve_struct(c_matrix,rows,cols,c);
}


void GHMM_DoubleMatrix::random_preserve_struct(double **matrix, int rows, int cols) {
  random_preserve_struct(c_matrix,rows,cols);
}


int GHMM_DoubleMatrix::gaussrows_values(GHMM_DoubleVector* mue, double u) {
  return matrix_d_gaussrows_values(c_matrix,rows,cols,&mue->c_vector,u);
}


GHMM_DoubleMatrix* GHMM_DoubleMatrix::transpose() {
  GHMM_DoubleMatrix* new_matrix = new GHMM_DoubleMatrix(cols,rows);

  matrix_d_transpose(c_matrix,rows,cols,new_matrix->c_matrix);

  return new_matrix;
}


/*
GHMM_DoubleVector* GHMM_DoubleMatrix::cholesky(double **a, double *b, int dim) {
  return matrix_cholesky(a,b,dim,x);
}
*/


/*
double GHMM_DoubleMatrix::det_symposdef() {
  return matrix_det_symposdef(c_matrix,dim,det);
}
*/

GHMM_DoubleMatrix* GHMM_DoubleMatrix::copy() {
  GHMM_DoubleMatrix* new_matrix = new GHMM_DoubleMatrix(rows,cols);

  matrix_d_copy(c_matrix,new_matrix->c_matrix,rows,cols);

  return new_matrix;
}
