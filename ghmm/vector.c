/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/vector.c
  created      : TIME: 17:48:29     DATE: Thu 18. February 1999
  $Id$

__copyright__

*******************************************************************************/

#include <float.h>
#include <stdio.h>
#include "mes.h"
#include "vector.h"
#include "rng.h"


/*============================================================================*/
/* Scales the elements of a vector to have the sum 1 */
/* PROBLEM: Entries can get very small and be rounded to 0 */
int vector_normalize(double *v, int len) {
#define CUR_PROC "vector_normalize"
  int i;
  double sum = 0.0;
  for (i = 0; i < len; i++)
    sum += v[i];
  if (sum < DBL_MIN) {
    mes_prot("Can't normalize vector. Sum eq. zero \n");
    return (-1);
  }
  for (i = 0; i < len; i++)
    v[i] /= sum;
  return 0;
#undef CUR_PROC
} /* vector_normalize */


/*============================================================================*/
void vector_const_values(double *v, int len, double c) {
  int i;
  for (i = 0; i < len; i++)
    v[i] = c;
} /* vector_const_values */


/*============================================================================*/
void vector_const_preserve_struct(double *v, int len, double c) {
  int i;
  for (i = 0; i < len; i++) {
    if (v[i] != 0.0)
      v[i] = c;
  }
} /* vector_const_preserve_struct */


/*============================================================================*/
void vector_random_values(double *v, int len) {
  int i;
  for (i = 0; i < len; i++)
    v[i] = gsl_rng_uniform(RNG);
} /* vector_random_values */


/*============================================================================*/
void vector_random_preserve_struct(double *v, int len) {
  int i;
  for (i = 0; i < len; i++) {
    if (v[i] != 0.0)
      v[i] = gsl_rng_uniform(RNG);
  }
} /* vector_random_preserve_struct */


/*============================================================================*/
void vector_d_print(FILE *file, double *vector, int len, 
		    char *tab, char *separator, char *ending) {
  int j;
  
  fprintf(file, "%s", tab);
  if (len > 0)
    fprintf(file, "%6.2f", vector[0]);

  for (j = 1; j < len; j++) 
    fprintf(file, "%s %6.2f", separator, vector[j]);
  fprintf(file, "%s\n", ending);
} /* vector_d_print */

/*============================================================================*/
void vector_d_print_prec(FILE *file, double *vector, int len, int width, 
			 int prec, char *tab, char *separator, char *ending) {
  int j;
  if (len > 0) fprintf(file, "%s%*.*f", tab, width, prec, vector[0]);
  for (j = 1; j < len; j++) 
    fprintf(file, "%s %*.*f", separator, width, prec, vector[j]);
  fprintf(file, "%s\n", ending);
} /* vector_d_print */

/*============================================================================*/
void vector_i_print(FILE *file, int *vector, int len,
		    char *tab, char *separator, char *ending) {
  int j;
  fprintf(file, "%s", tab);
  if (len > 0)
    fprintf(file, "%3d", vector[0]);
  for (j = 1; j < len; j++)
    fprintf(file, "%s %3d", separator, vector[j]);
  fprintf(file, "%s\n", ending);
} /* vector_i_print */

/*============================================================================*/
int vector_mat_times_vec(double **A, double *x, int n, int m, double *v) {
#define CUR_PROC "vector_mat_times_vec"
  int i, j;
  for (i = 0; i < n; i++){
    v[i] = 0.0;
    for (j = 0; j < m; j++)
      v[i] += A[i][j]*x[j];
  }
  return 0;
#undef CUR_PROC
} /* vector_mat_times_vec */




