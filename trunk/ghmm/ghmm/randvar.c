/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/randvar.c
  created      : TIME: 16:40:03     DATE: Wed 17. February 1999
  $Id$

__copyright__

*******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <math.h>
#ifdef HAVE_LIBPTHREAD
# include <pthread.h>
#endif /* HAVE_LIBPTHREAD */
#include "mes.h"
#include "mprintf.h"
#include "randvar.h"
#include "rng.h"
#include "scanner.h"
#include "float.h"
#include "const.h"

/* Liste fuer vorberechnete Werte der Dichtefkt. einer N(0,1)-Verteilung 
   Wertebereich x=0,00 - 19.99 */
#define PDFLEN 2000
#define X_STEP_PDF 0.01  /* Schrittweite */
#define X_FAKT_PDF 100   /* aequivalent zur Schrittweite */
static double pdf_stdnormal[PDFLEN];
static int pdf_stdnormal_exists = 0;

/* Liste fuer vorberechnete Werte PHI der Gaussverteilung wird eingelesen
   Wertebereich x = -9,999 - 0 */
#define X_STEP_PHI 0.001  /* Schrittweite */
#define X_FAKT_PHI 1000   /* aequivalent zur Schrittweite */

static double* PHI = NULL;
static int PHI_len = 0;
static double x_PHI_1 = -1.0;
static double x_PHI_xy = -1.0;

/*----------------------------------------------------------------------------*/
static int randvar_read_PHI() {
# define CUR_PROC "randvar_read_PHI"
  int res = -1;
  char filename[] = PHI_DATA_FILE;
  #warning "PHI_DATA_FILE depreciated!"
  scanner_t *s = NULL;
  s = scanner_alloc(filename); if(!s) {mes_proc(); goto STOP;}
  scanner_get_name(s);
  scanner_consume(s, '='); if(s->err) goto STOP; 
  if (!strcmp(s->id, "PHI")) {
    scanner_consume(s, '{'); if (s->err) goto STOP;
    PHI = scanner_get_double_earray(s, &PHI_len); if (s->err) goto STOP;
    scanner_consume(s, ';'); if (s->err) goto STOP;
    scanner_consume(s, '}'); if (s->err) goto STOP;
    scanner_consume(s, ';'); if (s->err) goto STOP;
  }
  else {
    scanner_error(s, "unknown identifier"); goto STOP;
  }
  /* printf("%.4f\n", PHI[PHI_len-1]); */

  res = 0;;
STOP:
  scanner_free(&s);
  return res;
# undef CUR_PROC
} /* randvar_read_PHI */


/*----------------------------------------------------------------------------*/
static int randvar_init_PHI() {
# define CUR_PROC "randvar_init_PHI"
#ifdef HAVE_LIBPTHREAD
  static pthread_mutex_t lock;
#endif /* HAVE_LIBPTHREAD */
  /* PHI einlesen */
  if (!PHI_len) {
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_lock(&lock); /* Setzen eines Locks, da Clustern parallel */
#endif /* HAVE_LIBPTHREAD */
    if (randvar_read_PHI() == -1) {mes_proc(); goto STOP;};
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_unlock(&lock); /* Freigabe des Locks */
#endif /* HAVE_LIBPTHREAD */
  }
  return 0;
STOP:
  return(-1);
# undef CUR_PROC
} /* randvar_init_PHI */


/*============================================================================*/
double randvar_get_xfaktphi() {
  return X_FAKT_PHI;
}

/*============================================================================*/
double randvar_get_xstepphi() {
  return X_STEP_PHI;
}

/*============================================================================*/
double randvar_get_philen() {
  return PHI_len;
}

/*============================================================================*/
double randvar_get_PHI(double x) {
# define CUR_PROC "randvar_get_PHI"
  int i;
  double phi_x;

  if (randvar_init_PHI() == -1) {mes_proc(); goto STOP;}

  /* lineare Interpolation (Alternative: Runden mit i=m_int(fabs(x)*X_FAKT))*/
  i = (int)(fabs(x) * X_FAKT_PHI); 
  if (i >= PHI_len-1) {
    i = PHI_len-1;
    phi_x = PHI[i];
  }
  else
    phi_x = PHI[i] + (fabs(x) - i*X_STEP_PHI) * (PHI[i+1]-PHI[i])/X_STEP_PHI;
  /* ACHTUNG: PHI fuer negative Werte tabelliert! */
  if (x > 0.0)
    return(1.0 - phi_x);
  else
    return(phi_x);

STOP:
  return(-1.0);
# undef CUR_PROC
} /* randvar_get_PHI */


/*============================================================================*/
/* Wann wird PHI[x,0,1] == 1? */
double randvar_get_xPHIless1() {
# define CUR_PROC "randvar_get_xPHIless1"
  double x;
  int i;
  if (x_PHI_1 == -1) {
    if (randvar_init_PHI() == -1) {mes_proc(); goto STOP;}    
    /* der letzte Wert der Tabelle hat auf jeden Fall den Wert 1 */
    for (x = (PHI_len-1)*X_STEP_PHI, i=PHI_len-1; i > 0; x -= X_STEP_PHI, i--)
      if (randvar_get_PHI(-x) > 0.0)
	break;
    /* Modifikation: x genau zwischen 2 Stuetzstellen! */
    x_PHI_1 = x - (double)X_STEP_PHI/2.0;
  }
  return(x_PHI_1);
STOP:
  return(-1.0);
# undef CUR_PROC
}

/*============================================================================*/
/* Wann wird PHI[x,0,1] ==  PHI[y,0,1] fuer hintereinanderliegende x,y ?*/
double randvar_get_xPHIxgleichPHIy() {
# define CUR_PROC "randvar_get_xPHIxgleichPHIy"
  double x,y;
  int i;
  if (x_PHI_xy == -1) {
    if (randvar_init_PHI() == -1) {mes_proc(); goto STOP;}     
    y = -1.0;
    for (x = 0.0, i=0; i < PHI_len; x += X_STEP_PHI, i++) {
      if (randvar_get_PHI(-x) == randvar_get_PHI(-y))
	break;
      y = x;
    }
    x_PHI_xy = y;
  }
  return(x_PHI_xy);
STOP:
  return(-1.0);
# undef CUR_PROC
}

/*============================================================================*/
double randvar_get_1durcha(double x, double mean, double u) {
  /* Berechnung von 1/a(x, mean, u) 
     mit a: Integral von x bis \infty ueber die Gaussdichte */
# define CUR_PROC "randvar_get_1durcha"
  int i;
  double y, z, phi_z, a;

  if (randvar_init_PHI() == -1) {mes_proc(); goto STOP;};

  y = 1/sqrt(u);
  z = (x - mean) * y; 
  /* lineare Interpolation (Alternative: Runden mit i=m_int(fabs(z)*X_FAKT))*/
  i = (int)(fabs(z)*X_FAKT_PHI); 
 
  if (i >= PHI_len-1) { 
    i = PHI_len-2;
    /* urspruenglich:
       i = PHI_len-1; letzter Wert im Table ist aber null */
    phi_z = PHI[i];
  }
  else
    phi_z = PHI[i] + (fabs(z)-i*X_STEP_PHI) * (PHI[i+1]-PHI[i])/X_STEP_PHI;
  /* ACHTUNG: PHI fuer negative Werte tabelliert! */
  if (z > 0.0) {
    if (phi_z ==  0) {
      mes_proc(); goto STOP;
    }
    else
      a = 1 / phi_z; /* PHI zwischen 0.5 und 1 */ /*??? zwischen 0.5 und 0 ! */   
  }
  else {
    a = 1 - phi_z;
    if (a > DBL_MIN) a = 1 / a;
    else {
      a = 0.0;
      mes(MES_WIN, "a ~= 0.0 critical! (mue = %.2f, u =%.2f)\n",
	  mean, u); /* goto STOP; */
    }
  }
  return a;
STOP:
  return(-1.0);
# undef CUR_PROC
} /* randvar_get_1durcha */


/*============================================================================*/
/* BEMERKUNG: 
   Berechnung dieser Dichtefunktion getestet, indem fuer beliebige
   mue und u die Integralsumme berechnet wurde:
     for (x = 0, x < ..., x += step(=0.01/0.001/0.0001)) 
       isum += step * randvar_normal_density_pos(x, mue, u);
   In allen Testfaellen "konvergierte" isum deutlich gegen 1 ! 
   (BK, 14.6.99)
   AENDERUNG:
   Bei -EPS_NDT (const.h) stutzen, damit x=0 keine Probleme bereitet
   (BK, 15.3.2000)
*/ 
double randvar_normal_density_pos(double x, double mean, double u) { 
# define CUR_PROC "randvar_normal_density_pos"
  return randvar_normal_density_trunc(x, mean, u, -EPS_NDT);
# undef CUR_PROC
} /* double randvar_normal_density_pos */


/*============================================================================*/
double randvar_normal_density_trunc(double x, double mean, double u,
				    double a) {
# define CUR_PROC "randvar_normal_density_trunc"
double c;

  if (u <= 0.0) {
    mes_prot("u <= 0.0 not allowed\n"); goto STOP;
  }
  if (x < a) return(0.0);
  
  if ((c = randvar_get_1durcha(a, mean, u)) == -1) 
    {mes_proc(); goto STOP;};
  
  return(c * randvar_normal_density(x, mean, u));
STOP:
  return(-1.0);
# undef CUR_PROC
} /* double randvar_normal_density_trunc */


/*----------------------------------------------------------------------------*/
static int randvar_init_pdf_stdnormal() {
# define CUR_PROC "randvar_init_pdf_stdnormal"
  int i;
  double x = 0.00;
  for (i = 0; i < PDFLEN; i++) {
    pdf_stdnormal[i] = 1/(sqrt(2*PI)) * exp( -1 * x * x /2);
    x += (double)X_STEP_PDF;
  }
  pdf_stdnormal_exists = 1;
  /* printf("pdf_stdnormal_exists = %d\n", pdf_stdnormal_exists); */
  return(0);
# undef CUR_PROC
} /* randvar_init_pdf_stdnormal */


/*============================================================================*/
double randvar_normal_density(double x, double mean, double u) { 
# define CUR_PROC "randvar_normal_density"
  double expo;
  if (u <= 0.0) {
    mes_prot("u <= 0.0 not allowed\n"); goto STOP;
  }
  /* evtl Nenner < EPS??? abfangen ? */
  expo = exp( -1 * m_sqr(mean - x) / (2 * u));
  return( 1/(sqrt(2*PI*u)) * expo );
STOP:
  return(-1.0);
# undef CUR_PROC
} /* double randvar_normal_density */


/*============================================================================*/
double randvar_normal_density_approx(double x, double mean, double u) { 
# define CUR_PROC "randvar_normal_density_approx"
#ifdef HAVE_LIBPTHREAD
  static pthread_mutex_t lock;
#endif /* HAVE_LIBPTHREAD */
  int i;
  double y, z, pdf_x;
  if (u <= 0.0) {
    mes_prot("u <= 0.0 not allowed\n"); goto STOP;
  }
  if (!pdf_stdnormal_exists) {
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_lock(&lock); /* Setzen eines Locks, da Clustern parallel */
#endif /* HAVE_LIBPTHREAD */
    randvar_init_pdf_stdnormal();
#ifdef HAVE_LIBPTHREAD
    pthread_mutex_unlock(&lock); /* Freigabe des Locks */
#endif /* HAVE_LIBPTHREAD */
  }
  y = 1/sqrt(u);
  z = fabs((x - mean)*y);
  i = (int)(z * X_FAKT_PDF);
  /* lineare Interpolation: */
  if (i >= PDFLEN-1) {
    i = PDFLEN-1;
    pdf_x = y * pdf_stdnormal[i];
  }
  else
    pdf_x = y * ( pdf_stdnormal[i] + 
		  (z - i*X_STEP_PDF) * 
		  (pdf_stdnormal[i+1] - pdf_stdnormal[i]) / X_STEP_PDF );
  return(pdf_x);
STOP:
  return(-1.0);
# undef CUR_PROC
} /* double randvar_normal_density_approx */


/*============================================================================*/
double randvar_std_normal(int seed) {
# define CUR_PROC "randvar_std_normal"
  static int first = 0;
  static double h1;
  static double h2;
  double U, V;
  if (seed != 0) {
    gsl_rng_set(RNG,seed);    
    return(1.0);
  }

  return( gsl_ran_gaussian(RNG, 1.0) );

# undef CUR_PROC
} /* randvar_std_normal */


/*============================================================================*/
double randvar_normal(double mue, double u, int seed) {
# define CUR_PROC "randvar_normal"
  double x;
  x = sqrt(u) * randvar_std_normal(seed) + mue;
  return(x);
# undef CUR_PROC
} /* randvar_normal */


/*============================================================================*/
#define C0 2.515517
#define C1 0.802853
#define C2 0.010328
#define D1 1.432788
#define D2 0.189269
#define D3 0.001308

double randvar_normal_pos(double mue, double u, int seed) {
# define CUR_PROC "randvar_normal_pos"
  double x = -1;
  double sigma, U, Us, Us1, Feps, Feps1, t, T;

  if (u <= 0.0) {
    mes_prot("u <= 0.0 not allowed\n"); goto STOP;
  }

  if (seed != 0) {
    gsl_rng_set(RNG,seed);    
    return(1.0);
  }

  /* Methode: solange Gauss-verteilte Zuf.Zahl erzeugen (mit GSL-Lib.), 
     bis sie im pos. Bereich liegt -> nicht effektiv, wenn mue << 0 
  while (x < 0.0) {
    x = sqrt(u) * randvar_std_normal(seed) + mue;
  } */
  
  /** Inverse Transformierung mit restricted sampling nach Fishman */
  sigma = sqrt(u);
  U = gsl_rng_uniform(RNG);               /*gsl_ran_flat(RNG,0,1) ??? */
  Feps = randvar_get_PHI(-(EPS_NDT+mue)/sigma);
  Us = Feps + (1-Feps)*U;
  /* num. besser: 1-Us = 1-Feps - (1-Feps)*U, deshalb: 
     Feps1 = 1-Feps, Us1 = 1-Us */
  Feps1 = randvar_get_PHI((EPS_NDT+mue)/sigma);
  Us1 = Feps1 - Feps1*U;
  t = m_min(Us,Us1);
  t = sqrt(-log(t*t));
  T = sigma * (t - (C0 + t*(C1 + t*C2)) / (1 + t*(D1 + t*(D2 + t*D3))) );
  if (Us-0.5 < 0) 
    x = mue - T;
  else
    x = mue + T;

STOP:
  return(x);
# undef CUR_PROC
} /* randvar_normal_pos */


/*============================================================================*/
double randvar_uniform_int(int seed, int K) {
# define CUR_PROC "randvar_uniform_int"
  if (seed != 0) {
    gsl_rng_set(RNG,seed);
    return(1.0);
  }
  else {
    double x = gsl_ran_flat(RNG, -0.5, K-0.5);
    x = m_int(x); /* m_int rundet auf Integer */
    return(x);
  }
# undef CUR_PROC
} /* randvar_uniform_int */


/*============================================================================*/
/* cumalative distribution function of N(mean, u) */
double randvar_normal_cdf(double x, double mean, double u) { 
# define CUR_PROC "randvar_normal_cdf"
  if (u <= 0.0) {
    mes_prot("u <= 0.0 not allowed\n"); goto STOP;
  }
  /* evtl Nenner < EPS abfangen ? */
  return(randvar_get_PHI((x - mean)/sqrt(u)));
STOP:
  return(-1.0);
# undef CUR_PROC
} /* double randvar_normal_cdf */

/*============================================================================*/
/* cumalative distribution function of -EPS_NDT-truncated N(mean, u) */
double randvar_normal_pos_cdf(double x, double mean, double u) { 
# define CUR_PROC "randvar_normal_pos_cdf"
  double Fx, c;
  if (x <= 0.0) return(0.0);
  if (u <= 0.0) { mes_prot("u <= 0.0 not allowed\n"); goto STOP; }
  /* evtl Nenner < EPS abfangen ? */
  Fx = randvar_get_PHI((x - mean)/sqrt(u));
  c = randvar_get_1durcha(-EPS_NDT, mean, u);
  return(c*(Fx-1)+1);
 STOP:
  return(-1.0);
# undef CUR_PROC
} /* double randvar_normal_cdf */


