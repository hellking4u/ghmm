/*
 * created: 19 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "GHMM_Transition.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Transition::GHMM_Transition(XMLIO_Attributes &attrs) {
  reading = GHMM_TRANSITION_READING_NONE;
  source = attrs["source"];
  target = attrs["target"];
}


GHMM_Transition::~GHMM_Transition() {
}


const char* GHMM_Transition::toString() const {
  return "GHMM_Transition";
}


XMLIO_Element* GHMM_Transition::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "prob") {
    reading = GHMM_TRANSITION_READING_PROB;

    return this;
  }

  fprintf(stderr,"tag '%s' not recognized in transition element\n",tag.c_str());
  exit(1);
  
  return NULL;
}


void GHMM_Transition::XMLIO_endTag(const string& tag) {
  reading = GHMM_TRANSITION_READING_NONE;
}


void GHMM_Transition::XMLIO_getCharacters(const string& characters) {
  switch (reading) {

  case GHMM_TRANSITION_READING_PROB:
    prob = atof(characters.c_str());
    break;
    
  case GHMM_TRANSITION_READING_NONE:
    break;
  }
}
