/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "GHMM_Emission.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Emission::GHMM_Emission() {
  mue      = 0;
  variance = 0;
  weight   = 0;
}


GHMM_Emission::~GHMM_Emission() {
}


const char* GHMM_Emission::toString() const {
  return "GHMM_Emission";
}


XMLIO_Element* GHMM_Emission::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  bool found = false;

  if (tag == "gauss") {
    mue      = atof(attrs["mue"].c_str());
    variance = atof(attrs["variance"].c_str());
    density  = normal;
    found    = true;
  }
  if (tag == "gauss-positive") {
    mue      = atof(attrs["mue"].c_str());
    variance = atof(attrs["variance"].c_str());
    density  = normal_pos;
    found    = true;
  }
  if (tag == "gauss-approximated") {
    mue      = atof(attrs["mue"].c_str());
    variance = atof(attrs["variance"].c_str());
    density  = normal_approx;
    found    = true;
  }

  if (! found) {
    fprintf(stderr,"Tag unrecognized in <emission> element: %s\n",tag.c_str());
    exit(1);
  }

  return NULL;
}


void GHMM_Emission::XMLIO_getCharacters(const string& characters) {
  weight = atof(characters.c_str());
}
