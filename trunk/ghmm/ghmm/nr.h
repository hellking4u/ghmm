#ifndef NR_H
#define NR_H

/**
   @memo
   @name modified functions from numerical recipies
   should be avoided
 */

//@{

///
double zbrent_AB(double (*func)(double, double, double, double),
		 double x1, double x2, double tol, double A, double B, 
		 double eps);

///
double rtsafe(void (*funcd)(double, double *, double *), double x1, double x2,
	      double xacc);

///
void lnsrch(int n, double xold[], double fold, double g[], double p[], 
	    double x[], double *f, double stpmax, int *check, 
	    double (*func)(double []));

///
void dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	    double(*func)(double []), void (*dfunc)(double [], double []));

///
void broydn(double x[], int n, int *check,
	    void (*vecfunc)(int, double [], double []));

/** Chi-Square-Test */
double gammln(double xx);
///
void gser(double *gamser, double a, double x, double *gln);
///
void gcf(double *gammcf, double a, double x, double *gln);
///
double gammq(double a, double x);
///
void chsone(double bins[], double ebins[], int nbins, int knstrn, double *df,
	    double *chsq, double *prob);
///
void chstwo(double bins1[], double bins2[], int nbins, int knstrn, double *df,
	    double *chsq, double *prob);

/** Kolmogorov-Smirnov-Test */
void sort(unsigned long n, double arr[]);
///
void sort2(unsigned long n, double arr[], double brr[]);
///
double probks(double alam);
///
void ksone(double data[], unsigned long n, double (*func)(double), double *d,
	   double *prob);
///
void kstwo(double data1[], unsigned long n1, double data2[], unsigned long n2,
	   double *d, double *prob);

/** Integration */
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
///
double trapzd(double (*func)(double), double a, double b, int n);
///
double qromb(double (*func)(double), double a, double b);
     

/* modified functions for SHMM */
#include "smodel.h"

///
void ksone_weighted_shmm(double data[], double weights[], unsigned long n,
			 double (*func)(smodel*, int, double), smodel *smo,
			 int state, double *d, double *prob);
///
double trapzd_cdfstate_shmm(double (*func)(smodel*, int, double,double),
			    smodel* smo, int state, double c, double a, 
			    double b, int n);
///
double qromb_cdfstate_shmm(double (*func)(smodel*,int,double,double), 
			   smodel *smo,int state,double c,double a,double b);

//@}

#endif /* NR_H */
