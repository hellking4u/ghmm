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
  Normierung der Summe der Elemente eines Vektors auf Eins..
  @return 0/-1
  @param v    Vektor
  @param len  Laenge des Vektors       
*/
int vector_normalize(double *v, int len);

/**
  Elemente eines Vektors mit einer Konstanten belegen
  @return -
  @param v    Vektor
  @param len  Laenge des Vektors
  @param c    konstanter Vorgabewert
*/
void vector_const_values(double *v, int len, double c);

/**
  Nicht-Null-Elemente eines Vektors mit einer Konstanten belegen
  @return -
  @param v    Vektor
  @param len  Laenge des Vektors       
  @param c    konstanter Vorgabewert
*/
void vector_const_preserve_struct(double *v, int len, double c);

/**
  Elemente eines Vektors zufällig (0..1) belegen
  @return -
  @param v    Vektor
  @param len  Laenge des Vektors       
*/
void vector_random_values(double *v, int len);

/**
  Nicht-Null-Elemente eines Vektors zufällig (0..1) belegen
  @return -
  @param v    Vektor
  @param len  Laenge des Vektors       
*/
void vector_random_preserve_struct(double *v, int len);

/**
  Schreiben eines double-Vektors (ohne Klammerung!)..
  @param file       Ausgabedatei
  @param vector     zu schreibender Vektor
  @param len        Dimension
  @param tab        Formatierung: fuehrende Tabs
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile  
  */
void vector_d_print(FILE *file, double *vector, int len, 
		    char *tab, char *separator, char *ending);

/**
  Schreiben eines double-Vektors (ohne Klammerung!) mit angegebenen Nachkommast.
  @param file       Ausgabedatei
  @param vector     zu schreibender Vektor
  @param len        Dimension
  @param width      Formatierung: Stellen insgesamt
  @param prec       Formatierung: Anzahl der Nachkommastellen
  @param tab        Formatierung: fuehrende Tabs
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile  
  */
void vector_d_print_prec(FILE *file, double *vector, int len, int width,
			 int prec, char *tab, char *separator, char *ending);

/**
  Schreiben eines int-Vektors (ohne Klammerung!)..
  @param file       Ausgabedatei
  @param vector     zu schreibender Vektor
  @param len        Dimension
  @param tab        Formatierung: fuehrende Tabs
  @param separator  Formatierung: Trennzeichen der Spalten
  @param ending     Formatierung: Endzeichen der Zeile  
  */
void vector_i_print(FILE *file, int *vector, int len,
		    char *tab, char *separator, char *ending);
/**
  Berechnen eines double-Vektors aus Matrix * Vektor
  @param A       Matrix n x m
  @param x       Vektor zur Berechnung
  @param n       Anzahl Zeilen
  @param m       Anzahl Spalten
  @param v       berechneter Vektor (Rueckgabe)
  */
int vector_mat_times_vec(double ** A, double *x, int n, int m, double *v);

#ifdef __cplusplus
}
#endif


#endif /* VECTOR_H */

/*@} (Doc++-Group: vector) */
