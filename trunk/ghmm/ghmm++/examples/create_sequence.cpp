/* @(#)create_sequence.cpp created by Peter Pipenbacher at 25 Apr 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 * __copyright__
 */

#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


int main() {
  /* Important! initialise rng  */
  GHMM_Toolkit::gsl_rng_init();

  vector<double> v;
  v.push_back(0.5);
  v.push_back(0.35);
  v.push_back(0.75);

  GHMM_Sequence seq(GHMM_DOUBLE,v.size(),1);
  for (unsigned int i = 0; i < v.size(); ++i)
    seq.setDouble(i,v[i]);

  seq.print(stdout);
  
#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return 0;
}
