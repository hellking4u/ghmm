/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include "assert.h"
#include <xmlio/XMLIO_Document.h>
#include "ghmm++/GHMM_GMLDiscreteModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


void GHMM_GMLDiscreteModel::init()
{
  tag = "graph";
  // no attribute
  attributes.clear();
  setNodeTag("node");
  setTransitionTag("edge");
}

GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel() : GHMM_DiscreteModel() 
{
  init();
  gmlalphabet = NULL;
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(model* my_model) : GHMM_DiscreteModel(my_model)
{
  init();
  gmlalphabet = NULL;
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(int number_of_states, int my_M, double my_prior): GHMM_DiscreteModel(number_of_states, my_M, my_prior)
{
  init();
  gmlalphabet = NULL;
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(GHMM_GMLAlphabet* my_alphabet) : GHMM_DiscreteModel(my_alphabet)
{
  init();
  gmlalphabet = my_alphabet;
}

GHMM_GMLAlphabet* GHMM_GMLDiscreteModel::getAlphabet() const
{
  return gmlalphabet;
}

void GHMM_GMLDiscreteModel::createTransitions() {
  unsigned int i;
  int edge;

  for (i = 0; i < states.size(); ++i)
    {
      for (edge = 0; edge < states[i]->getOutEdges(); ++edge)
      {
	GHMM_GMLState *gmlstate_pt = new GHMM_GMLState(this, i+1, states[i]->c_state);
	transitions.push_back(gmlstate_pt->createTransition(edge));
	delete gmlstate_pt;
      }
    }
}

const int GHMM_GMLDiscreteModel::XMLIO_writeContent(XMLIO_Document& writer) 
{
  int total_bytes = 0;
  int result;
  unsigned int i;

  writer.changeIndent(2);
  result = writer.writeEndl();

  if (result < 0)
    return result;
  total_bytes += result;


  /* write states */
  for (i = 0; i < GHMM_DiscreteModel::states.size(); ++i) {
    if (i > 0) {
      result = writer.writeEndl();

      if (result < 0)
	return result;
      total_bytes += result;
    }

    result = writer.writeElement(GHMM_DiscreteModel::states[i]);
  
    if (result < 0)
      return result;
    total_bytes += result;

    result = writer.writeEndl();

    if (result < 0)
      return result;
    total_bytes += result;
  }

  /* write transitions */
  createTransitions();
  for (i = 0; i < transitions.size(); ++i) {
    writer.writeEndl();
    writer.writeElement(transitions[i]);
    writer.writeEndl();
  }
  cleanTransitions();

  return 0;
}




