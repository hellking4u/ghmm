/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : ghmm/ghmm/randvar.h
  created      : TIME: 16:43:38     DATE: Wed 17. February 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/

#ifndef RANDVAR_H
#define RANDVAR_H

/**@name Hilfsfunktionen für Zufallsvariablen */
/*@{ (Doc++-Group: randvar) */


#define PI 3.141592653589793116

/**
   Berechnen der eindimensionalen Dichtefkt(mean, u) an der Stelle x..
   @return       Funktionswert
   @param x      Variablenwert
   @param mean   Mittelwert der Normalverteilung
   @param u      Varianz der Normalverteilung ($\sigma^2$)
   */
double randvar_normal_density(double x, double mean, double u);

/**
   Bestimmung der eindimensionalen Dichtefkt(mean, u) an der Stelle x.
   Der Wert wird aus einer vorberechneten Liste geholt, die erstellt wird
   beim ersten Aufruf.
   @return       Funktionswert
   @param x      Variablenwert
   @param mean   Mittelwert der Normalverteilung
   @param u      Varianz der Normalverteilung ($\sigma^2$)
   */
double randvar_normal_density_approx(double x, double mean, double u);

/**
   Bestimmung der eindimensionalen, im negativen Bereich auf 0 gesetzten
   Gauss-Dichte-Fkt(mean, u) an der Stelle x.
   @return       Funktionswert
   @param x      Variablenwert
   @param mean   Mittelwert der Normalverteilung
   @param u      Varianz der Normalverteilung ($\sigma^2$)
   */
double randvar_normal_density_pos(double x, double mean, double u);

/**
   Bestimmung der eindimensionalen, im Bereich ($-\infty$,a) auf 0 gesetzten
   Gauss-Dichte-Fkt(mean, u) an der Stelle x.
   @return       Funktionswert
   @param x      Variablenwert
   @param mean   Mittelwert der Normalverteilung
   @param u      Varianz der Normalverteilung ($sigma^2$)
   @param a     linke Grenze der Dichte
   */
double randvar_normal_density_trunc(double x, double mean, double u, double a);

/** 
   Erzeugen einer U(0,K-1) verteilten ganzzahligen Zufallszahl..
   @return          Zufallszahl
   @param seed      1: nur initialisieren des U(0,1) Generators, 
                    0: Rueckgabe der gleichverteilten ganzz. Zufallszahl 
  */
double randvar_uniform_int(int seed, int K);

/** 
   Erzeugen einer N(0,1) verteilten Zufallszahl..
   @return          Zufallszahl
   @param seed      1: nur initialisieren des U(0,1) Generators,
                    0: Rueckgabe einer standardnormalverteilen Zufallszahl 
  */
double randvar_std_normal(int seed);

/** 
   Erzeugen einer N(mue, u) verteilten Zufallszahl..
   @return          Zufallszahl
   @param seed      1: nur initialisierung
                    0: Rueckgabe einer standardnormalverteilen Zufallszahl 
  */
double randvar_normal(double mue, double u, int seed);

/** 
   Erzeugen einer positiven N(mue, u) verteilten Zufallszahl..
   @return          Zufallszahl
   @param seed      1: nur initialisieren des U(0,1) Generators,
                    0: Rueckgabe einer standardnormalverteilen Zufallszahl 
  */
double randvar_normal_pos(double mue, double u, int seed);

/**
   Bestimmen der N(0,1)-Verteilungsfunktion an der Stelle x. 
   Die Verteilung wird als Tabelle eingelesen und dann wird zwischen 
   den St""utzstellen linear interpoliert.
   @return          Wert PHI(x)
   @param x         x-Wert
   */
double randvar_get_PHI(double x);

/**
   Bestimmen des Skalierungsfaktors 1/a fuer die gestutzte Normalverteilung. 
   a entspricht dem Integral von x bis Infinity ueber die Normaldichte.
   @return          1/a
   @param x         linke Integralgrenze
   @param mean      Mittelwert der Normaldichte
   @param u         Varianz der Normaldichte
   */
double randvar_get_1durcha(double x, double mean, double u);

/**
   Bestimmen derjenigen Stuetzstelle x, fuer die PHI(x) das erste Mal 1 wird..
   @return          x mit PHI(x)==1, PHI(y) < 1 fuer alle y < x
   */
double randvar_get_xPHIless1();

/**
   Bestimmen derjenigen Stuetzstelle x, fuer die PHI(x) 
   das erste Mal gleich PHI(y) wird, wenn x und y aufeinander folgen.
   @return          x mit PHI(x)==PHI(y) mit PHI(x')<PHI(y') f.a.x',y'<x,y
   */
double randvar_get_xPHIxgleichPHIy();


/**                                 
   cumalative distribution function F(x;mean,u) der NV(mean, u)..
   @return          F(x)
   @param x         x
   @param mean      Mittelwert der Normaldichte
   @param u         Varianz der Normaldichte
*/
double randvar_normal_cdf(double x, double mean, double u);
/**                                 _
   cumalative distribution function F(x;mean,u) der gestutzten NV(mean, u)..
   @return          F(x)
   @param x         x
   @param mean      Mittelwert der Normaldichte
   @param u         Varianz der Normaldichte
*/
double randvar_normal_pos_cdf(double x, double mean, double u);

double randvar_get_xfaktphi();

double randvar_get_xstepphi();

double randvar_get_philen();

#endif /* RANDVAR_H */

/*@} (Doc++-Group: randvar) */
