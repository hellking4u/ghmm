/* @(#)GHMM_DoubleMatrix.h created by Peter Pipenbacher at 04 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_DOUBLEMATRIX_H
#define _GHMM_DOUBLEMATRIX_H 1

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_DoubleMatrix;
class GHMM_DoubleVector;

/** */
class GHMM_DoubleMatrix {

 public:

  /** Constructor. */
  GHMM_DoubleMatrix(int rows, int cols, double default_value = 0.0);
  /** Destructor. */
  virtual ~GHMM_DoubleMatrix();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Solves a linear equation system, Ax = b, for a symmetric, positiv definit matrix.
     @param a:    double matrix 
     @param b:    double vector
     @param dim:  dimension of a
     @return double vector, a solution of the system (NULL on failure)
  */
  //  GHMM_DoubleVector* cholesky(double **a, double *b, int dim);
  /** 
      Gives all elements != 0 in a matrix a constant value.
      @param c:        value for the elements
  */
  void const_preserve_struct(double c);
  /** 
      Gives all elements in a matrix a constant value.
      @param c:        value for the elements
  */  
  void const_values(double c);
  /** Copies a matrix. */
  GHMM_DoubleMatrix* copy();
  /**
     Finds the determinant of a symmetric, positiv definit matrix.
     @return 0 for success; -1 for failure
     @param dim:  dimension of a
     @return param det:  determinant of a, the returning value
  */
  //  double det_symposdef();
  /** 
      Gives each row in a matrix values according to a certain Gauss density.
      Mean values are randomly generated if mue == NULL.
      u ($\sigma^2$) must be given.
      @return 0 for success; -1 for failure
      @param mue:      pointer to the vector containing the mean values for each row
      @param u:        standard deviation, for all rows equal
  */  
  int gaussrows_values(GHMM_DoubleVector* mue, double u);
  /** 
      Scales the rowvectors of a matrix, so that they have sum 1.
      @return 0 for success; -1 for error
  */
  int normalize();
  /**
     Determines the number of entries != 0 in a row of a matrix.
     @return         number of entries
  */
  int notzero_columns();
  /**
     Determines the number of entries != 0 in a column of a matrix. 
     @return         number of entries
  */
  int notzero_rows();
  /**
     Writes a double matrix (without parenthesis).
     @param file:       output file
     @param tab:        format: leading tabs
     @param separator:  format: separator for columns
     @param ending:     format: end of a row  
  */
  void print(FILE *file, char *tab, char *separator, char *ending);
  /**
     Writes a double matrix (without parenthesis) with specifically many decimal places.
     @param file:       output file
     @param width:      format: number of places altogether
     @param prec:       format: number of decimal places
     @param tab:        format: leading tabs
     @param separator:  format: separator for columns
     @param ending:     format: end of a row
  */
  void print_prec(FILE *file, int width, int prec, char *tab, 
		  char *separator, char *ending);
  /** 
      Gives the elements in a matrix uniformly distributed values between min and max. 
      Gives all elements in the last row a constant value.
      @param min:      minimum for the random values
      @param max:      maximum for the random values
      @param c:        value for the last row
  */  
  void random_const_values(double min, double max, double c);
  /** Gives the elements in a matrix with band width 3 random values. */  
  void random_left_right(double **matrix, int rows, int cols);
  /** Gives all elements != 0 in a matrix uniformly distributed random values 
      between 0 and 1. */
  void random_preserve_struct(double **matrix, int rows, int cols);
  /** 
      Gives the elements in a matrix uniformly distributed values between min and max. 
      @param min:      minimum for the random values
      @param max:      maximum for the random values
  */  
  void random_values(double min, double max);
  /** Gives all elements on the 1. upper secondary diagonal in a matrix 
      the value 1. */  
  void left_right_strict();
  /**
     Calculates matrix x vec, where vec is a double vector
     @param vec            vector to calculate
     @return calculated vector
  */
  GHMM_DoubleVector* times_vec(GHMM_DoubleVector* vec);
  /** 
      Transposes a matrix.
      @return transposed double matrix
  */  
  GHMM_DoubleMatrix* transpose();

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

#include <ghmm++/close_code.h>

#endif /* _GHMM_DOUBLEMATRIX_H */
