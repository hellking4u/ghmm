/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#include "ghmm++/GHMM_GMLContinuousModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLContinuousModel::GHMM_GMLContinuousModel() : GHMM_ContinuousModel() 
{ 
  setNodeTag("node");
  setTransitionTag("edge");

  /* Do not initialize anything. It will be read from xml file. */
}


GHMM_GMLContinuousModel::GHMM_GMLContinuousModel(int N, int M, int cos, density_t density,
						 double prior) 
  : GHMM_ContinuousModel(N, M, cos, density, prior) 
{
  /* Set keywords for XML-element tags */
  setNodeTag("node");
  setTransitionTag("edge");
}

