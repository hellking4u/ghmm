/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#include <xmlio/XMLIO_Document.h>
#include "GHMM_GMLTransition.h"
#include "GHMM_Toolkit.h"
#include "GHMM_State.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLTransition::GHMM_GMLTransition(XMLIO_Attributes &attrs) : GHMM_Transition(attrs) {
  reading           = GHMM_TRANSITION_READING_NONE;
  source            = attrs["source"];
  target            = attrs["target"];
  attributes        = attrs;
  tag               = "edge";
  xmlio_indent_type = XMLIO_INDENT_BOTH;
}


GHMM_GMLTransition::GHMM_GMLTransition(GHMM_State* my_source, GHMM_State* my_target, double my_prob) : GHMM_Transition(my_source, my_target, my_prob) {
  reading              = GHMM_TRANSITION_READING_NONE;
  source               = my_source->index;
  target               = my_target->index;
  prob                 = my_prob;
  attributes["source"] = my_source->id;
  attributes["target"] = my_target->id;
  tag                  = "edge";
  xmlio_indent_type    = XMLIO_INDENT_BOTH;
}


const char* GHMM_GMLTransition::toString() const {
  return "GHMM_GMLTransition";
}


XMLIO_Element* GHMM_GMLTransition::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "data") 
    {
      if (attrs["key"] == "prob") {
	reading = GHMM_TRANSITION_READING_PROB;
	return this;
      }
    }
  fprintf(stderr,"tag '%s' not recognized in transition element\n",tag.c_str());
  exit(1);
  
  return NULL;
}


const int GHMM_GMLTransition::XMLIO_writeContent(XMLIO_Document& writer) {
  int result = 0;

  writer.changeIndent(2);
  result += writer.writef("\n%s<data key=\"prob\">%f</data>\n",writer.indent,prob);

  return result;
}
