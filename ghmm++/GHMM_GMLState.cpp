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
#include "ghmm++/GHMM_GMLDiscreteModel.h"
#include "ghmm++/GHMM_GMLDataNode.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, XMLIO_Attributes& attrs):
  GHMM_State(my_model, my_index, attrs)
{
  tag               = "node"; // state
  hasData["ngeom"]  = false;
}


GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, sstate* my_state) :
  GHMM_State(my_model, my_index, my_state)
{
  tag               = "node"; // state
  hasData["ngeom"]  = false;
}

GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, state* my_state) :
  GHMM_State(my_model, my_index, my_state)
{
  tag               = "node"; // state
  hasData["ngeom"]  = false;
}


const char* GHMM_GMLState::toString() const {
  return "GHMM_GMLState";
}


GHMM_GMLTransition* GHMM_GMLState::createTransition(int edge_index) {
  if (c_state)
    return new GHMM_GMLTransition(this,parent_model->getState(c_state->out_id[edge_index]),c_state->out_a[edge_index]);

  if (c_sstate)
    return new GHMM_GMLTransition(this,parent_model->getState(c_sstate->out_id[edge_index]),c_sstate->out_a[0][edge_index]);

  return NULL;
}


XMLIO_Element* GHMM_GMLState::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  
  if (tag == "data") /* For GraphML */
    {
      if (attrs["key"] == "initial") {
	reading = GHMM_STATE_INITIAL;
	hasData["initial"] = true;
	return this;
      }

      if (attrs["key"] == "label") {
	reading = GHMM_STATE_LABEL;
	hasData["label"] = true;
	return this;
      }
      if (attrs["key"] == "emissions") {
	hasData["emissions"] = true;
	attributes["key"] = "emissions";
	SAFE_DELETE(emission);
	emission = new GHMM_Emission(this);
	return emission;
      }

      if (attrs["key"] == "ngeom") {
	attributes["key"] = "ngeom";
	return this;
      }
    }
  else
    if (tag == "pos")
      {
	hasData["ngeom"] = true;
	vPosition[0] = atof(attrs["x"].c_str());
	vPosition[1] = atof(attrs["y"].c_str());
	return this;
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

float GHMM_GMLState::get2DPosition(int index) {
  return vPosition[index];
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

  result = writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;

  GHMM_GMLDataNode* my_label = new GHMM_GMLDataNode(this, "label");
  result = writer.writeElement(my_label);
  SAFE_DELETE(my_label);
  result += writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;

  //  }

  if ( hasData["ngeom"] == true )
    {
      GHMM_GMLDataNode* my_pos = new GHMM_GMLDataNode(this, "ngeom");
      result = writer.writeElement(my_pos);
      SAFE_DELETE(my_pos);
    }
  else
    {
      // Warning or put in some random positions
      // HMMEd expected a file to contain coordinate, or else we should change it
    }

  /* write emission */
  result += writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;
  
  GHMM_GMLDataNode* my_emission = new GHMM_GMLDataNode(this, "emissions");
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








