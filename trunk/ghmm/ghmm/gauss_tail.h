#ifndef GAUSS_TAIL_H
#define GAUSS_TAIL_H

/**
   @name some calculations concerning gaussian tail function
 */

//@{
///
double pmue(double mue, double A, double B, double eps);

///
double pmue_umin(double mue, double A, double B, double eps);

/** 
    @name Funktion zum Loesen der Nullstellengleichungen fuer abgeschnittenen Normaldichte 
*/
double pmue_interpol(double mue, double A, double B, double eps);

//@}
#endif /* GAUSS_TAIL_H */
