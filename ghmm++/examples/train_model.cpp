/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main() {
  GHMM_Document doc;

  doc.open("test3.xml","r");
  doc.readDocument();
  doc.close();

  //  GHMM_DiscreteModel* model = doc.getDiscreteModel();

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return 0;
}
