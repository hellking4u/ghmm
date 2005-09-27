/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/matrix.h
*       Authors:  Bernhard Knab, Benjamin Georgi
*
*       Copyright (C) 1998-2004 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik, 
*                               Berlin
*                                   
*       Contact: schliep@ghmm.org             
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision$ 
*                       from $Date$
*             last change by $Author$.
*
*******************************************************************************/
#ifndef MATRICE_H
#define MATRICE_H

#include <stdio.h>
#include "scanner.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@name Matrix */
/*@{ (Doc++-Group: matrix) */

/**
  Allocation of a double matrix. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  double **matrix_d_alloc (int n, int m);


/**
  Allocation of a 3 dimensional double matrix. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  @param height: 3rd dimension
  */
double *** matrix3d_d_alloc(int i, int j, int k);
int matrix3d_d_free(double **** matrix, int i, int j);

int *** matrix3d_i_alloc(int i, int j, int k);
int matrix3d_i_free(int **** matrix, int i, int j);
/**
  Copying and allocation of a double matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  @param copymatrix: matrix to copy 
  */
  double **matrix_d_alloc_copy (int rows, int columns, double **copymatrix);

/**
  Free the memory of a double matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int matrix_d_free (double ***matrix, long zeilen);

/**
  Allocation of a static double matrix with a single malloc. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  double **stat_matrix_d_alloc (int n, int m);

/**
  Free the memory of a static double matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int stat_matrix_d_free (double ***matrix);

/**
  Allocation of a static int matrix with a single malloc. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  int **stat_matrix_i_alloc (int n, int m);

/**
  Free the memory of a static int matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int stat_matrix_i_free (int ***matrix);


/**
  Allocation of a integer matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
  int **matrix_i_alloc (int rows, int columns);

/**
  Free the memory of a integer matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
  int matrix_i_free (int ***matrix, long rows);


#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: matrix) */
