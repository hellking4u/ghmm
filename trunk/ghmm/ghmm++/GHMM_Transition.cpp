/*
 * created: 19 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include <xmlio/XMLIO_Document.h>
#include "GHMM_Transition.h"
#include "GHMM_Toolkit.h"
#include "GHMM_State.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Transition::GHMM_Transition(XMLIO_Attributes &attrs) {
  reading           = GHMM_TRANSITION_READING_NONE;
  source            = attrs["source"];
  target            = attrs["target"];
  attributes        = attrs;
  tag               = "transition";
  xmlio_indent_type = XMLIO_INDENT_BOTH;
}


GHMM_Transition::GHMM_Transition(GHMM_State* my_source, GHMM_State* my_target, double my_prob) {
  reading              = GHMM_TRANSITION_READING_NONE;
  source               = my_source->index;
  target               = my_target->index;
  prob                 = my_prob;
  attributes["source"] = my_source->id;
  attributes["target"] = my_target->id;
  tag                  = "transition";
  xmlio_indent_type    = XMLIO_INDENT_BOTH;
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
    // printf("\tTransition prob: %g\n", prob);
    break;
    
  case GHMM_TRANSITION_READING_NONE:
    break;
  }
}


const int GHMM_Transition::XMLIO_writeContent(XMLIO_Document& writer) {
  int result = 0;

  writer.changeIndent(2);
  result += writer.writef("\n%s<prob>%f</prob>\n",writer.indent,prob);

  return result;
}
