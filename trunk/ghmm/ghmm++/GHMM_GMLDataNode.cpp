/*
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 *
 * 
 * __copyright__
 */


#include <xmlio/XMLIO_Document.h>
#include <xmlio/XMLIO_Element.h>
#include "ghmm++/GHMM_GMLDataNode.h"
#include "ghmm++/GHMM_DiscreteModel.h"
#include "ghmm++/GHMM_GMLState.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLDataNode::GHMM_GMLDataNode(GHMM_GMLState *my_state, const char *keyvalue)
{
  state      = my_state;
  tag        = "data";
  this->keyvalue   = string(keyvalue);
  attributes["key"] = keyvalue;

  xmlio_indent_type = XMLIO_INDENT_BOTH;
}


const int GHMM_GMLDataNode::XMLIO_writeContent(XMLIO_Document& writer)
{
  int total_bytes = 0;
  int i, result;

  writer.changeIndent(2);
  
  if (keyvalue == "label")
    {
      result = writer.writef("%s",state->label.c_str());
      if (result < 0)
	return result;
      total_bytes += result;
    }

  /* discrete model */
  if (state->c_state)
  {
    if (keyvalue == "emissions")
      {
	//GHMM_Alphabet* alphabet   = state->getModel()->getAlphabet();
	GHMM_DiscreteModel* model = (GHMM_DiscreteModel*) state->getModel();
	for (i = 0; i < model->c_model->M - 1; ++i) {
	  result += writer.writef("%.2f, ",state->c_state->b[i]); // space is important	    
	}
	result += writer.writef("%.2f",state->c_state->b[i]);	   
	total_bytes += result;
      }    
  }

  if (keyvalue == "ngeom")
    {    
      XMLIO_Element    pos_element;
      char tmp[15];
      sprintf(tmp, "%5.2f", state->get2DPosition(0));
      pos_element.attributes["x"]=tmp;
      sprintf(tmp, "%5.2f", state->get2DPosition(1));
      pos_element.attributes["y"]=tmp;
      pos_element.tag = "pos";
      result = writer.writeElement(&pos_element);
      if (result < 0)
	return result;
      total_bytes += result;
    }
  return total_bytes;
}



