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
