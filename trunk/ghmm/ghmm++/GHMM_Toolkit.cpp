/*
 * created: 06 Mar 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ghmm/rng.h"
#include "ghmm++/GHMM_Toolkit.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


void GHMM_Toolkit::gsl_rng_init() {
  ::gsl_rng_init();
}


string GHMM_Toolkit::toString(int var) {
  char mem[100];
  sprintf(mem,"%d",var);
  return string(mem);
}


string GHMM_Toolkit::toString(double var) {
  char mem[100];
  sprintf(mem,"%g",var);
  return string(mem);
}
