/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ppghmm++/GHMM_Emission.h"
#include "ppghmm++/GHMM_State.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Emission::GHMM_Emission(GHMM_State* my_state) {
  mue             = 0;
  variance        = 0;
  weight          = -1;
  function_loaded = false;
  state           = my_state;
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
    fprintf(stderr,"<emission> element of state '%s' has unrecognized tag '%s'.\n",
	    state->id.c_str(),tag.c_str());
    exit(1);
  }

  if (attrs["mue"] == "") {
    fprintf(stderr,"<%s> element of state '%s' lacks mue attribute.\n",
	    tag.c_str(),state->id.c_str());
    exit(1);
  }

  if (attrs["variance"] == "") {
    fprintf(stderr,"<%s> element of state '%s' lacks variance attribute.\n",
	    tag.c_str(),state->id.c_str());
    exit(1);
  }

  function_loaded = true;

  return NULL;
}


void GHMM_Emission::XMLIO_getCharacters(const string& characters) {
  if (weight != -1) {
    fprintf(stderr,"Emission of state '%s' has multiple weights.\n",state->id.c_str());
    exit(1);
  }
  if (function_loaded) {
    fprintf(stderr,"Emission of state '%s' has character data behind density function.\n",state->id.c_str());
    exit(1);
  }
  weight = atof(characters.c_str());
}


void GHMM_Emission::XMLIO_finishedReading() {
  /* weight 1 is default value. */
  if (weight == -1)
    weight = 1;

  if (!function_loaded) {
    fprintf(stderr,"Emission of state '%s' lacks density function.\n",state->id.c_str());
    exit(1);
  }
}
