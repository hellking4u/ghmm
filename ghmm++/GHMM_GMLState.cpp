/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <iostream>

#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Document.h>
#include "ghmm++/GHMM_GMLState.h"
#include "ghmm++/GHMM_Emission.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, XMLIO_Attributes& attrs):
  GHMM_State(my_model, my_index, attrs)
{
  tag               = "node"; // state
}


GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, sstate* my_state) :
  GHMM_State(my_model, my_index, my_state)
{
  tag               = "node"; // state
}

GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, state* my_state) :
  GHMM_State(my_model, my_index, my_state)
{
  tag               = "node"; // state
}


const char* GHMM_GMLState::toString() const {
  return "GHMM_GMLState";
}


XMLIO_Element* GHMM_GMLState::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  
  if (tag == "data") /* For GraphML */
    {
      if (attrs["key"] == "initial") {
	reading = GHMM_STATE_INITIAL;
	return this;
      }

      if (attrs["key"] == "label") {
	reading = GHMM_STATE_LABEL;
	return this;
      }
      if (attrs["key"] == "emissions") {
	SAFE_DELETE(emission);
	emission = new GHMM_Emission(this);
	return emission;
      }
    }
  else
    { 
      fprintf(stderr,"GMLState :: tag '%s' not recognized in state element\n",tag.c_str());
      exit(1);
    }

  return NULL;
}


void GHMM_GMLState::XMLIO_endTag(const string& tag) {
  reading = GHMM_STATE_NONE;
}


void GHMM_GMLState::XMLIO_getCharacters(const string& characters) {
  switch (reading) {

  case GHMM_STATE_INITIAL:
    initial = atof(characters.c_str());
    break;

  case GHMM_STATE_LABEL:
    label = characters.c_str();
    cout << "State label " << label << endl;
    break;
    
  case GHMM_STATE_NONE:
    break;
  }
}

XMLIO_Attributes& GHMM_GMLState::XMLIO_getAttributes() {
  attributes["id"] = id;
  return attributes;
}


const int GHMM_GMLState::XMLIO_writeContent(XMLIO_Document& writer) {
  int total_bytes = 0;
  int result;
  
  double initial = getInitial();

  writer.changeIndent(2);

  //  if (initial > 0) {
  writer.writeEndlIndent();
  
  if (result < 0)
    return result;
  total_bytes += result;
  
  result = writer.writef("<data key=\"initial\">%f</data>",initial);
  
  if (result < 0)
    return result;
  total_bytes += result;
  //  }

  /* write emission */
  result = writer.writeEndl();
  
  if (result < 0)
    return result;
  total_bytes += result;
  
  GHMM_Emission* my_emission = new GHMM_Emission(this);
  result = writer.writeElement(my_emission);
  SAFE_DELETE(my_emission);
  
  if (result < 0)
    return result;
  total_bytes += result;
  /* end write emission */

  result = writer.writeEndl();

  if (result < 0)
    return result;
  total_bytes += result;

  return total_bytes;
}








