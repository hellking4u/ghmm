/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/vector.h
  created      : TIME: 17:46:40     DATE: Thu 18. February 1999
  $Id$

__copyright__

*******************************************************************************/

#ifndef VECTOR_H
#define VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>

/**@name Vektor-Funktionen */
/*@{ (Doc++-Group: vector) */

/**
  Scales the sum of elements in a vector to one.
  @return 0 for success; -1 for error
  @param v    vector
  @param len  length of the vector       
*/
int vector_normalize(double *v, int len);

/**
  Gives all elements in a vector a constant value
  @param v    vector
  @param len  length of the vector
  @param c    given value for the elements
*/
void vector_const_values(double *v, int len, double c);

/**
  Gives all elements, not equal zero, in a vector a constant value
  @param v    vector
  @param len  length of the vector
  @param c    given value for the elements
*/
void vector_const_preserve_struct(double *v, int len, double c);

/**
  Gives all elements in a vector random values between 0 and 1
  @param v    vector
  @param len  length of the vector       
*/
void vector_random_values(double *v, int len);

/**
  Gives all elements, not equal zero, in a vector random values between 0 and 1
  @param v    vector
  @param len  length of the vector   
*/
void vector_random_preserve_struct(double *v, int len);

/**
  Writes a double vector (without parenthesis)
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row  
  */
void vector_d_print(FILE *file, double *vector, int len, 
		    char *tab, char *separator, char *ending);

/**
  Writes a double vector (without parenthesis) with given number of decimal places
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param width      format: total number of decimal places
  @param prec       format: number of decimal places after the comma
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row 
  */
void vector_d_print_prec(FILE *file, double *vector, int len, int width,
			 int prec, char *tab, char *separator, char *ending);

/**
  Writes an integer vector (without parenthesis)
  @param file       output file
  @param vector     vector to write
  @param len        dimension
  @param tab        format: leading tabs
  @param separator  format: separator for columns
  @param ending     format: end of a row  
  */
void vector_i_print(FILE *file, int *vector, int len,
		    char *tab, char *separator, char *ending);
/**
  Calculates Ax, where A is a double matrix and x a double vector
  @param A       n x m matrix
  @param x       vector to calculate
  @param n       number of rows
  @param m       number of columns
  @param v       calculated vector (return value)
  */
int vector_mat_times_vec(double ** A, double *x, int n, int m, double *v);

#ifdef __cplusplus
}
#endif


#endif /* VECTOR_H */

/*@} (Doc++-Group: vector) */
