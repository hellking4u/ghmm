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
#include "ghmm++/GHMM_GMLDiscreteModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel() : GHMM_DiscreteModel() 
{
  setNodeTag("node");
  setTransitionTag("edge");
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(model* my_model) : GHMM_DiscreteModel(my_model)
{
  setNodeTag("node");
  setTransitionTag("edge");
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(int number_of_states, int my_M, double my_prior): GHMM_DiscreteModel(number_of_states, my_M, my_prior)
{
  setNodeTag("node");
  setTransitionTag("edge");
}


GHMM_GMLDiscreteModel::GHMM_GMLDiscreteModel(GHMM_Alphabet* my_alphabet) : GHMM_DiscreteModel(my_alphabet)
{
  setNodeTag("node");
  setTransitionTag("edge");
}




