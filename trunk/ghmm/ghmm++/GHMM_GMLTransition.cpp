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
#include "GHMM_StateT.hh"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLTransition::GHMM_GMLTransition(XMLIO_Attributes &attrs) :
GHMM_TransitionT<GHMM_GMLState>(attrs){
  reading           = GHMM_TRANSITION_READING_NONE;
  source            = attrs["source"];
  target            = attrs["target"];
  attributes        = attrs;
  tag               = "edge";
  xmlio_indent_type = XMLIO_INDENT_BOTH;
}


GHMM_GMLTransition::GHMM_GMLTransition(GHMM_GMLState* my_source, GHMM_GMLState* my_target, vector<double> my_prob) 
: GHMM_TransitionT<GHMM_GMLState>(my_source, my_target)
{
  reading              = GHMM_TRANSITION_READING_NONE;
  source               = my_source->index;
  target               = my_target->index;
  probs                = my_prob;
  attributes["source"] = my_source->id;
  attributes["target"] = my_target->id;
  tag                  = "edge";
  xmlio_indent_type    = XMLIO_INDENT_BOTH;
}

const char* GHMM_GMLTransition::toString() const {
  return "GHMM_GMLTransition";
}

void GHMM_GMLTransition::XMLIO_getCharacters(const string& characters) {
  unsigned int pos;
  switch (reading) {

  case GHMM_TRANSITION_READING_PROB:    
    for (pos = 0; pos < characters.size(); ++pos) {
      while (pos < characters.size() && isspace(characters[pos]))
		++pos;
      
      if (pos < characters.size())
		probs.push_back(atof(characters.substr(pos).c_str()));
      
      while (pos < characters.size() && !isspace(characters[pos]))
		++pos;
    }
    break;
    
  case GHMM_TRANSITION_READING_NONE:
    break;
  }
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
  int i;

  writer.changeIndent(2);
  result = writer.writeEndl();
  result += writer.writef("<data key=\"prob\">");
  for(i=0; i < (int)(probs.size())-1; i++)
    {
      result += writer.writef(" %.10f,", probs[i]); // A leading space is important
    }
  result += writer.writef(" %.10f", probs[i]); // A leading space is important
  result += writer.writef("</data>\n",prob);
  return result;
}
