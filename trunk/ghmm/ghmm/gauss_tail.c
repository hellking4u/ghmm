/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/gauss_tail.c
  created      : March 2001 by Achim Gaedke from hmm/src/creestimate.c
  $Id$

__copyright__

*******************************************************************************/

#include <float.h>
#include <math.h>
#include "const.h"

/*============================================================================*/
double pmue(double mue, double A, double B, double eps) {
  double feps, u, Atil, Btil;
  Atil = A + eps;
  Btil = B + eps*A;
  u = Btil - mue*Atil;
  //if (u < EPS_U) u = (double)EPS_U; ACHTUNG: wuerde Fkt.wert verfaelschen!
  if (u <= DBL_MIN)
    return(mue - A);
  feps = randvar_normal_density_trunc(-eps, mue, u, -eps);
  return(A - mue - u*feps);
}

/*============================================================================*/
/* pmue zur Vermeidung von numerischen Oszillationen:
   Interpolation von p(\mu) selbst zwischen 2 St""utzstellen fuer PHI
   BEACHTE: 1.Version, sehr aufwendig und um die Ecke
            -> spaeter vereinfachen! 
*/
double pmue_interpol(double mue, double A, double B, double eps) {
  double u, Atil, Btil, z,z1,z2,m1,m2,u1,u2,p1,p2,pz;
  int i1,i2;
  Atil = A + eps;
  Btil = B + eps*A;
  u = Btil - mue*Atil;
  //if (u < EPS_U) u = (double)EPS_U; ACHTUNG: wuerde Fkt.wert verfaelschen!
  if (u <= DBL_MIN)
    return(mue - A);

  /* im positiven Bereich von mue Berechnung wie gehabt */
  if (mue >= 0.0)
    return(A - mue - u*randvar_normal_density_trunc(-eps, mue, u, -eps));

  /* sonst: Interpolation der Funktion selbst zwischen 2 Stuetzstellen */
  z = (eps + mue)/sqrt(u);
    
  i1 = (int)(fabs(z) * randvar_get_xfaktphi());
  if (i1 >= randvar_get_philen()-1)
    i1 = i2 = randvar_get_philen()-1;
  else
    i2 = i1+1;
  z1 = i1/randvar_get_xfaktphi();
  z2 = i2/randvar_get_xfaktphi();

  m1 = -z1 * sqrt(Btil + eps*Atil + Atil*Atil*z1*z1*0.25) 
    - (eps + Atil*z1*z1*0.5);
  m2 = -z2 * sqrt(Btil + eps*Atil + Atil*Atil*z2*z2*0.25) 
    - (eps + Atil*z2*z2*0.5);
  u1 = Btil - m1*Atil;
  u2 = Btil - m2*Atil;

  p1 = A - m1 - u1*randvar_normal_density_trunc(-eps, m1, u1, -eps);
  p2 = A - m1 - u1*randvar_normal_density_trunc(-eps, m2, u2, -eps);

  if (i1 >= randvar_get_philen()-1)
    pz = p1;
  else {
    pz = p1 + (fabs(z)-i1*randvar_get_xstepphi()) 
      * (p2 - p1)/randvar_get_xstepphi();
    //pz = p1;
  }
  return(pz);
}

/*============================================================================*/
double pmue_umin(double mue, double A, double B, double eps) {
  double feps, u;
  u = EPS_U;
  feps = randvar_normal_density_trunc(-eps, mue, u, -eps);
  return(A - mue - u*feps);
}
