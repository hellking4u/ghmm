#ifndef MATRICE_H
#define MATRICE_H
#include "model.h"

/**@name Matrix */
/*@{ (Doc++-Group: matrix) */

/**
  Alloc einer double Matrix..
  @return Zeiger auf Matrix
  @param zeilen   Zeilenzahl
  @param spalten  Spaltenzahl
  */
double** matrix_d_alloc(int zeilen, int spalten);

/**
  Copy & Alloc einer double Matrix
  @return            Zeiger auf Matrix
  @param zeilen      Zeilenzahl
  @param spalten     Spaltenzahl
  @param copymatrix  zu kopierende Matrix
  */
double** matrix_d_alloc_copy(int zeilen, int spalten, double **copymatrix);

/**
  Free einer double Matrix..
  @return Erfolgsstatus: 0 succes; -1 error
  @param matrix:  freizugebende Matrix
  @param zeilen:  Zeilenzahl
  */
int matrix_d_free(double ***matrix, long zeilen);
/**
  Alloc einer int Matrix..
  @return Zeiger auf Matrix
  @param zeilen:  Zeilenzahl
  @param spalten: Spaltenzahl
  */
int** matrix_i_alloc(int zeilen, int spalten);
/**
  Free einer int Matrix..
  @return Erfolgsstatus: 0 succes; -1 error
  @param matrix:  freizugebende Matrix
  @param zeilen:  Zeilenzahl
  */
int matrix_i_free(int ***matrix, long zeilen); 

/**
  Schreiben einer double Matrix (ohne Klammerung!)..
  @param file       Ausgabedatei
  @param matrix     zu schreibende Matrix
  @param zeilen     Zeilenzahl
  @param spalten    Spaltenzahl
  @param tab        Formatierung  fuehrende Tabs
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile  
  */
void matrix_d_print(FILE *file, double **matrix, int zeilen, int spalten, 
		    char *tab, char *separator, char *ending);

/**
  Schreiben einer double Matrix (ohne Klammerung!) mit angegebenen Nachkommast.
  @param file       Ausgabedatei
  @param matrix     zu schreibende Matrix
  @param zeilen     Zeilenzahl
  @param spalten    Spaltenzahl
  @param width      Formatierung: Stellen insgesamt
  @param prec       Formatierung: Anzahl der Nachkommastellen
  @param tab        Formatierung  fuehrende Tabs/Zeichen
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile
  */
void matrix_d_print_prec(FILE *file, double **matrix, int zeilen, int spalten, 
			 int width, int prec, char *tab, char *separator, 
			 char *ending);

/**
  Schreiben einer int Matrix (ohne Klammerung!)..
  @param file       Ausgabedatei
  @param matrix     zu schreibende Matrix
  @param zeilen     Zeilenzahl
  @param spalten    Spaltenzahl
  @param tab        Formatierung: fuehrende Tabs
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile  
  */
void matrix_i_print(FILE *file, int **matrix, int zeilen, int spalten,
		    char *tab, char *separator, char *ending);
/**
  Einlesen einer double Matrix..
  @return Erfolgsstatus: 0 succes; -1 error
  @param s:       scanner
  @param matrix:  zu lesende Matrix
  @param zeilen:  Zeilenzahl
  @param spalten: Spaltenzahl
  */
int matrix_d_read(scanner_t *s, double **matrix, int max_zeile, int max_spalte);
/**
  Einlesen einer int Matrix..
  @return Erfolgsstatus: 0 succes; -1 error
  @param s:       scanner
  @param matrix:  zu lesende Matrix
  @param zeilen:  Zeilenzahl
  @param spalten: Spaltenzahl
  */
int matrix_i_read(scanner_t *s, int **matrix, int max_zeile, int max_spalte);

/**
  Anzahl der Eintraege != 0 einer Matrixzeile bestimmen..
  @return         Anzahl Eintraege
  @param matrix   Double-Matrix
  @param row      zu untersuchende Zeile
  @param max_col  Anzahl der Matrixspalten
  */
int matrix_d_notzero_columns(double **matrix, int row, int max_col);

/**
  Anzahl der Eintraege != 0 einer Matrixspalte bestimmen..
  @return         Anzahl Eintraege
  @param matrix   Double-Matrix
  @param col      zu untersuchende Spalte
  @param max_row  Anzahl der Matrixzeilen
  */
int matrix_d_notzero_rows(double **matrix, int col, int max_row);

/** 
  Normiert die Zeilenvektoren einer Matrix auf 1
  @return ??
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  */
int matrix_d_normalize(double **matrix, int rows, int cols);

/** 
  Matrix mit gleichverteilten Zufallszahlen zwischen min und max belegen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param min      Minimum für Zufallswert
  @param max      Maximum für Zufallswert
  */  
void matrix_d_random_values(double **matrix, int rows, int cols, 
			    double min, double max); 

/** 
  Matrix mit gleichverteilten Zufallszahlen zwischen min und max belegen..
  Letzte Zeile mit const. Wert c belegen.
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param min      Minimum für Zufallswert
  @param max      Maximum für Zufallswert
  @param c        Wert fuer letzte Zeile
  */  
void matrix_d_random_const_values(double **matrix, int rows, int cols,
				  double min, double max, double c);

/** 
  Matrix mit gleichverteilten Zufallszahlen zwischen min und max belegen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param c        konstanter Wert für alle Elemente
  */  
void matrix_d_const_values(double **matrix, int rows, int cols, double c); 

/** 
  Matrix mit 1en auf der 1. oberen Nebendiagonalen belegen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  */  
void matrix_d_left_right_strict(double **matrix, int rows, int cols); 

/** 
  Matrix auf der Haupt- und 1.oberen Nebendiagonalen mit Zufallszahlen belegen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  */  
void matrix_d_random_left_right(double **matrix, int rows, int cols); 


/** 
  Nicht-Null-Elemente einer Matrix auf konstanten Wert setzen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param c        vorgegebener konstanter Wert
  */
void matrix_d_const_preserve_struct(double **matrix,int rows,int cols,double c);

/** 
  Nicht-Null-Elemente einer Matrix mit gleichverteilten Zufallszahlen 
  zwischen 0 und 1 belegen..
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  */
void matrix_d_random_preserve_struct(double **matrix, int rows, int cols);

/** 
  Matrix pro Zeile nach Gaussdichte belegen.
  Mittelwerte werden zufaellig erzeugt, falls mue == NULL.
  u ($\sigma^2$) muss (global) vorgegeben werden.
  @param matrix   Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param mue      Adresse des Vektors der Mitelwerte pro Zeile (Erz.u.U.in Fkt)
  @param u    Std.abw. der Normalvert (für alle Zeilen gleich)
 */  
int matrix_d_gaussrows_values(double **matrix, int rows, int cols,
			    double **mue, double u);

/** 
  Matrix transponieren.
  @param A        Double-Matrix
  @param rows     Anzahl der Matrixzeilen
  @param cols     Anzahl der Matrixspalten
  @param A_T      Transponierte Double-Matrix (Rueckgabe)
 */  
void matrix_d_transpose(double **A, int rows, int cols, double **A_T);

int matrix_cholesky(double **a, double *b, int dim, double *x);
int matrix_det_symposdef(double **a, int dim, double *det);

/** copy a matrix. alloc needs to be done outside ! */
void matrix_d_copy(double **src, double **target, int rows, int cols);

/** calculate "target" as neighbour matrix from "src", preserve strucure and row
    norm */
void matrix_d_neighbour(double **src, double **target, int rows, int cols, double eps);


#endif

/*@} (Doc++-Group: matrix) */
