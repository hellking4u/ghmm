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
  vector<double> v;
  v.push_back(0.5);
  v.push_back(0.35);
  v.push_back(0.75);

  // init GHMM_Sequence data structure
  GHMM_Sequence seq1(GHMM_DOUBLE,v.size());
  for (unsigned int i = 0; i < v.size(); ++i)
    seq1.setDouble(i,v[i]);

  seq1.print(stdout);
  
  // create a new GHMM_Sequences data structure from 
  // the GHMM_Sequence object created above
  GHMM_Sequences seq2(&seq1);
  seq2.print(stdout);

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return 0;
}
