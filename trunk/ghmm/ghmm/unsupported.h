#ifdef GHMM_UNSUPPORTED

#ifndef GHMM_UNSUPPORTED_H
#define GHMM_UNSUPPORTED_H

/*================== matrix.h ===============================================*/
/**
  Reads in an integer matrix.
  @return 0 for succes; -1 for error
  @param s:          scanner
  @param matrix:     matrix to read
  @param max_row:    number of rows
  @param max_column: number of columns
  */
  int ighmm_dmatrix_read (scanner_t * s, int **matrix, int max_row,
                     int max_column);

/**
  Writes a double matrix (without parenthesis) with specifically many decimal places.
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param width:      format: number of places altogether
  @param prec:       format: number of decimal places
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row
  */
  void ighmm_cmatrix_print_prec (FILE * file, double **matrix, int rows,
                            int columns, int width, int prec, char *tab,
                            char *separator, char *ending);

/**
  Writes an integer matrix (without parenthesis).
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row  
  */
  void ighmm_dmatrix_print (FILE * file, int **matrix, int rows, int columns,
                       char *tab, char *separator, char *ending);

#endif /* GHMM_UNSUPPORTED_H */

#endif /* GHMM_UNSUPPORTED */
