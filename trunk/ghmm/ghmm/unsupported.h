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

/*================== mes.h ==================================================*/
  /**
   */
  void ighmm_mes_fformat (char *txt, char *filename, int line, char *proc_info);
  /**
   */
  int ighmm_mes_insert (FILE * fp, char src, int cnt);
  /**
   */
  int ighmm_mes_copy (char *oldname, char *newname);
  /**
   */
  int ighmm_mes_fgetc (FILE * fp);
  /**
   */
  int ighmm_mes_fflush (FILE * fp);
  /**
   */
  int ighmm_mes_fprintf (FILE * fp, char *format, ...);
  /**
   */
  int ighmm_mes_fputc (FILE * fp, char chr);
  /**
   */
  int ighmm_mes_fputs (FILE * fp, char *str);
  /**
   */
  int ighmm_mes_fread (FILE * fp, void *mem, int bytes);
  /**
   */
  int ighmm_mes_fread_quiet (FILE * fp, void *mem, int bytes);
  /**
   */
  int ighmm_mes_fseek (FILE * fp, long offset, int fromwhere);
#ifdef WIN32
  /**
   */
  int mes_fseek64 (FILE * fp, unsigned int uoff, unsigned int loff,
                   int fromwhere);
#endif
  /**
   */
  int ighmm_mes_ftell (FILE * fp);
  /**
   */
  int ighmm_mes_fwrite (FILE * fp, void *mem, int bytes);
  /**
   */
  int ighmm_mes_move (char *oldname, char *newname);
  /**
   */
  int ighmm_mes_remove (char *filename);
  /**
   */
  int ighmm_mes_rename (char *oldname, char *newname);
  /**
   */
  FILE *ighmm_mes_tmpfile (void);


#endif /* GHMM_UNSUPPORTED_H */

#endif /* GHMM_UNSUPPORTED */
