/*******************************************************************************
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/matrix.c
  created      : TIME: ?            DATE: ?
  last-modified: TIME: 17:44:48     DATE: Fri 26. January 2001
*******************************************************************************/
/* $Id$ */

#include <math.h>
#include <float.h>
#include "matrix.h"
#include "vector.h"
#include "model.h"
#include "rng.h"
#include "randvar.h"
#include "mes.h"

static void lrdecomp(int dim, double **a, double *p);
static void lyequalsb(double **a, double *b, double *p, int dim, double *y);
static void ltranspxequalsy(double **a, double *y, double *p, int dim, 
			    double *x);

/*============================================================================*/

int matrix_d_read(scanner_t *s, double **matrix, int max_zeile, int max_spalte){
#define CUR_PROC "matrix_d_read"
  int len =0, zeile = 0;
  scanner_consume(s, '{'); if (s->err) return(-1);
  while( !s->eof && !s->err && s->c - '}') {
    if (zeile >= max_zeile) {
      scanner_error(s, "too many rows in matrix");
      return(-1);
    } 
    /* Speicher Alloc erfolgt in scanner_get_double_array */
    matrix[zeile] = scanner_get_double_array(s, &len);
    if (len != max_spalte) {
      scanner_error(s, "wrong number of elements in matrix");
      return(-1);
    }
    scanner_consume(s, ';'); 
    if (s->err) {
      scanner_error(s, "missing ';' or wrong number of columns");
      return(-1);
    }
    zeile++;
  }
  scanner_consume(s, '}'); if (s->err) return(-1);
  if (zeile < max_zeile){
    scanner_error(s, "rows missing in matrix");
    return(-1);
  } 
  return(0);  
#undef CUR_PROC
} /* matrix_d_read */

/*============================================================================*/

int matrix_i_read(scanner_t *s, int **matrix, int max_zeile, int max_spalte) {
#define CUR_PROC "matrix_i_read" 
  int len =0, zeile = 0;
  scanner_consume(s, '{'); if (s->err) return(-1);
  while( !s->eof && !s->err && s->c - '}') {
    if (zeile >= max_zeile) {
      scanner_error(s, "too many rows in matrix");
      return(-1);
    } 
    /* Speicher Alloc erfolgt in scanner_get_int_array */
    matrix[zeile] = scanner_get_int_array(s, &len);
    if (len != max_spalte) {
      scanner_error(s, "wrong number of elements in matrix");
      return(-1);
    }
    scanner_consume(s, ';'); 
    if (s->err) {
      scanner_error(s, "missing ';' or wrong number of columns");
      return(-1);
    }
    zeile++;
  }
  scanner_consume(s, '}'); if (s->err) return(-1);
  if (zeile < max_zeile){
    scanner_error(s, "rows missing in matrix");
    return(-1);
  } 
  return(0);  
#undef CUR_PROC
} /* matrix_i_read */

/*============================================================================*/

double** matrix_d_alloc(int zeilen, int spalten) {
#define CUR_PROC "matrix_d_alloc"
  double **matrix;
  int i;
  if (!m_calloc(matrix, zeilen)) {mes_proc(); goto STOP;}
  for (i = 0; i < zeilen; i++)
    if (!m_calloc(matrix[i], spalten)) {mes_proc(); goto STOP;}
  return matrix;
STOP:
  matrix_d_free(&matrix, zeilen);
  return NULL;
#undef CUR_PROC
} /* matrix_d_alloc */

/*============================================================================*/

double** matrix_d_alloc_copy(int zeilen, int spalten, double **copymatrix) {
#define CUR_PROC "matrix_d_alloc_copy"
  double **matrix;
  int i, j;
  if (!m_calloc(matrix, zeilen)) {mes_proc(); goto STOP;}
  for (i = 0; i < zeilen; i++) {
    if (!m_calloc(matrix[i], spalten)) {mes_proc(); goto STOP;}
    for (j = 0; j < spalten; j++)
      matrix[i][j] = copymatrix[i][j];
  }
  return matrix;
STOP:
  matrix_d_free(&matrix, zeilen);
  return NULL;
#undef CUR_PROC
} /* matrix_d_alloc_copy */

/*============================================================================*/

int** matrix_i_alloc(int zeilen, int spalten) {
#define CUR_PROC "matrix_i_alloc"
  int **matrix;
  int i;
  if (!m_calloc(matrix, zeilen)) {mes_proc(); goto STOP;}
  for (i = 0; i < zeilen; i++)
    if (!m_calloc(matrix[i], spalten)) {mes_proc(); goto STOP;}
  return matrix;
STOP:
  matrix_i_free(&matrix, zeilen);
  return NULL;
#undef CUR_PROC
} /* matrix_i_alloc */

/*============================================================================*/

int matrix_i_free(int ***matrix, long zeilen) {
# define CUR_PROC "matrix_i_free"
  long i;
  mes_check_ptr(matrix, return(-1));
  if ( !*matrix ) return(0);
  for (i = 0; i < zeilen; i++) 
    m_free((*matrix)[i]);
  m_free(*matrix);
  return (0);
# undef CUR_PROC
} /* matrix_i_free */

/*============================================================================*/

int matrix_d_free(double ***matrix, long zeilen) {
# define CUR_PROC "matrix_d_free"
  long i;
  mes_check_ptr(matrix, return(-1));
  if ( !*matrix) return(0);
  for (i = zeilen - 1; i >=  0; i--) 
    m_free((*matrix)[i]);
  m_free(*matrix);
  return (0);
# undef CUR_PROC
} /* matrix_d_free */

/*============================================================================*/

void matrix_d_print(FILE *file, double **matrix, int zeilen, int spalten, 
		    char *tab, char *separator, char *ending) {
  int i;
  for (i = 0; i < zeilen; i++)
    vector_d_print(file, matrix[i], spalten, tab, separator, ending);
} /* matrix_d_print */

/*============================================================================*/

void matrix_d_print_prec(FILE *file, double **matrix, int zeilen, int spalten, 
			 int width, int prec, char *tab, char *separator, 
			 char *ending) {
  int i;
  for (i = 0; i < zeilen; i++)
    vector_d_print_prec(file, matrix[i], spalten, width, prec, 
			tab, separator, ending);
} /* matrix_d_print_prec */

/*============================================================================*/

void matrix_i_print(FILE *file, int **matrix, int zeilen, int spalten,
		    char *tab, char *separator, char *ending) {
  int i;
  for (i = 0; i < zeilen; i++)
    vector_i_print(file, matrix[i], spalten, tab, separator, ending);
} /* matrix_i_print */


/*============================================================================*/
int matrix_d_notzero_columns(double **matrix, int row, int max_col) {
  int i, count = 0;
  for (i =0; i < max_col; i++)
    if (matrix[row][i]) count++;
  return count;
} /* matrix_d_notzero_columns */

/*============================================================================*/
int matrix_d_notzero_rows(double **matrix, int col, int max_row) {
  int i, count = 0;
  for (i =0; i < max_row; i++)
    if (matrix[i][col]) count++;
  return count;
} /* matrix_d__notzero_rows */

/*============================================================================*/
/* Normiert die Zeilenvektoren einer Matrix auf Summe 1 */
int matrix_d_normalize(double **matrix, int rows, int cols) {
#define CUR_PROC "matrix_d_normalize"
  int i;
  for (i = 0; i < rows; i++)
    if (vector_normalize(matrix[i], cols) == -1)
      mes(MES_WIN, "WARNING: sum row[%d] == 0!\n", i);
      /* return (-1); */
  return 0;
#undef CUR_PROC
} /* matrix_d_normalize */

/*============================================================================*/
void matrix_d_random_values(double **matrix, int rows, int cols,
			    double min, double max) {
  int i, j;
  double interval;
  if (max < min) {
    min = 0.0;
    max = 1.0;
  }
  interval = max - min;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = min + gsl_rng_uniform(RNG)*interval;
} /* matrix_d_random_values */

/*============================================================================*/
/* fixed value for final state */
void matrix_d_random_const_values(double **matrix, int rows, int cols,
				  double min, double max, double c) {
  int i, j;
  double interval;
  if (rows < 1) {
    mes(MES_WIN, "WARNING: rows = %d not allowed\n", rows);
    return;
  }
  if (max < min) {
    min = 0.0;
    max = 1.0;
  }
  interval = max - min;
  for (i = 0; i < rows - 1; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = min + gsl_rng_uniform(RNG)*interval;
  for (j = 0; j < cols; j++)
    matrix[rows - 1][j] = c;
} /* matrix_d_random_const_values */


/*============================================================================*/
void matrix_d_const_values(double **matrix, int rows, int cols, double c) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      matrix[i][j] = c;
} /* matrix_d_const_values */

/*============================================================================*/
void matrix_d_random_left_right(double **matrix, int rows, int cols) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      if (j == i || j == i+1) 
	matrix[i][j] = gsl_rng_uniform(RNG);
      else
	matrix[i][j] = 0.0;
} /* matrix_d_random_values */

/*============================================================================*/
void matrix_d_left_right_strict(double **matrix, int rows, int cols) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      if (j == i+1)
	matrix[i][j] = 1.0;
      else
	matrix[i][j] = 0.0;
} /* matrix_d_left_right_strict */

/*============================================================================*/
int matrix_d_gaussrows_values(double **matrix, int rows, int cols,
			    double **mue, double u) {
# define CUR_PROC "matrix_gaussrows_values"
  int res = -1;
  double *mean;
  int i, j;
  if (u <= 0.0) {mes_prot("sigma^2 <= 0.0 not allowed\n"); goto STOP;}
  if (*mue == NULL) {
    /* fuer jede Zeile zufaelliger Mittelwert mean[i] in (0, cols-1) */
    if (!m_calloc(mean, rows)) {mes_proc(); goto STOP;}
    for (i = 0; i < rows; i++)
      mean[i] = gsl_rng_uniform(RNG) * (cols-1);
    /* for (i = 0; i < rows; i++) printf("%6.4f ", mean[i]); printf("\n"); */
    *mue = mean;
  }
  else
    /* Ueberpruefen, ob vorgegebene Mittelwerte im richtigen Intervall.. */
    mean = *mue;
  for (i = 0; i < rows; i++) {
    /* fuer jeden Zustand Gaussverteilung um Mittelwert */
    for (j = 0; j < cols; j++) {
      matrix[i][j] = randvar_normal_density((double)j, mean[i], u);
      if (matrix[i][j] == -1) {mes_proc(); goto STOP;}
      /* NULL VERMEIDEN: (1.billige Version) */
      if (matrix[i][j] < 0.0001)
	matrix[i][j] = 0.0001; /* Ausgabe nur 4 Nachkommastellen! */
    }
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* matrix_gaussrows_values */

/*============================================================================*/
void matrix_d_const_preserve_struct(double **matrix, int rows, int cols,
				    double c) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (matrix[i][j] != 0)
	matrix[i][j] = c;
    }
} /* matrix_d_const_preserve_struct */

/*============================================================================*/
void matrix_d_random_preserve_struct(double **matrix, int rows, int cols) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) {
      if (matrix[i][j] != 0)
	matrix[i][j] = gsl_rng_uniform(RNG);
    }
} /* matrix_d_random_preserve_struct */

/*============================================================================*/
void matrix_d_transpose(double **A, int rows, int cols, double **A_T) {
  int i, j;
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      A_T[j][i] = A[i][j];
} /* matrix_d_transpose */


/*----------------------------------------------------------------------------*/
/*--  zerlegt pos.def.sym.A in L*L-T, d.h. bis auf Diagonalel. wird         --*/
/*--  L (=untere Dreiecksmatrix) auf A gespeichert, Diag.el. auf Vec. p     --*/
static void lrdecomp(int dim, double **a, double *p) {
  int k,i,j;
  double x;
  for (i = 0; i < dim; i++) {
    for (j = i; j < dim; j++) {
      x = a[i][j];
      for (k = i-1; k >= 0; k--)
	x = x - a[j][k]*a[i][k];
      if (i == j){
	if (x < DBL_MIN) mes(MES_WIN, "FEHLER: Pivotel.<=0!");
	p[i] = 1 / sqrt(x);
      }
      else 
	a[j][i] = x*p[i];
    }
  }
}

/*----------------------------------------------------------------------------*/
/*--  loest L*y=b, L untere Dreieckmat. auf A gespeichert, p=1/Diag.elemente -*/
static void lyequalsb(double **a, double *b, double *p, int dim, double *y) {
  int k,j;
  for (k = 0; k < dim; k++) {
    y[k] = b[k];
    for (j = 0; j < k; j++)
      y[k] = y[k] - a[k][j]*y[j];
    y[k] = y[k]*p[k];
  }
}

/*----------------------------------------------------------------------------*/
/*--  loest L-T*x=y, L-T obere Dreieckmat, ABER als L in A gespeichert !    --*/
static void ltranspxequalsy(double **a, double *y, double *p, int dim, 
			    double *x) {
  int k,j;
  for (k = dim-1; k >= 0; k--) {
    x[k] = y[k];
    for (j = k+1; j < dim; j++)
      x[k] = x[k] - a[j][k]*x[j];
    x[k] = x[k]*p[k];
  }
}


/*============================================================================*/
/*-----------------  LGS-Loesen fuer sym.,pos.def.Matrix ---------------------*/
int matrix_cholesky(double **a, double *b, int dim, double *x) {
#define CUR_PROC "matrix_cholesky"
  int res = -1;
  double *p, *y;  
  if (!m_calloc(p, dim)) {mes_proc(); goto STOP;}
  if (!m_calloc(y, dim)) {mes_proc(); goto STOP;}
  lrdecomp(dim,a,p);
  lyequalsb(a,b,p,dim,y);
  ltranspxequalsy(a,y,p,dim,x);
  res = 0;
STOP:
  return(res);
#undef CUR_PROC
}

/*-----------------  Determinante einer sym.,pos.def.Matrix ------------------*/
int matrix_det_symposdef(double **a, int dim, double *det) {
#define CUR_PROC "matrix_det_symposdef"
  int res = -1;
  int i;
  double *p, r;  
  if (!m_calloc(p, dim)) {mes_proc(); goto STOP;}
  lrdecomp(dim,a,p);
  *det = 1.0;
  for (i = 0; i < dim; i++) {
    r = 1/p[i];
    *det *= (r*r);
  }
  res = 0;
STOP:
  return(res);
#undef CUR_PROC
}

/*============================================================================*/

/* --------------- copy a matrix, alloc needs to be done outside ! -------------*/
void matrix_d_copy(double **src, double **target, int rows, int cols) {
  int i, j;

  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++)
      target [ i ][ j ] = src [ i ][ j ];
}

/*============================================================================*/

/*  calculate "target" as neighbour matrix from "src", preserve strucure 
    and row  norm */
void matrix_d_neighbour(double **src, double **target, int rows, int cols, 
			double eps) {
  int i, k;
  double sum, tmp;

  for (i = 0; i < rows; i++) {
    sum = 1.0;
    for (k = 0; k < cols; k++) {
      if (src[i][k] > 0) {
	/*  Creates a random number on the interval [0, eps] */
	tmp = schange_rand_2(eps); 
	sum += tmp;
	target[i][k] = tmp + src[i][k];
      }
      else
	target[i][k] = 0.0;
    }
    
    /* norm to 1 */
    for (k = 0; k < cols; k++) 
      target[i][k] /= sum;
  }
}

/*============================================================================*/
