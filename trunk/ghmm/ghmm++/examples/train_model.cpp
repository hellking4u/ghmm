/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#include "ppghmm++/GHMM_DiscreteModel.h"
#include "ppghmm++/GHMM_Document.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main() {
  GHMM_Document doc;

  doc.open("test3.xml","r");
  doc.readDocument();
  doc.close();

  //  GHMM_DiscreteModel* model = doc.getDiscreteModel();

  return 0;
}
