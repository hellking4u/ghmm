/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "GHMM_AbstractModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_AbstractModel::GHMM_AbstractModel() {
}


GHMM_AbstractModel::~GHMM_AbstractModel() {
}


const char* GHMM_AbstractModel::toString() const {
  return "GHMM_AbstractModel";
}


void GHMM_AbstractModel::print(FILE *file) {
}
