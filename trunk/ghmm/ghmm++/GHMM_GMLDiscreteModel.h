/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#ifndef _GHMM_GMLDISCRETEMODEL_H
#define _GHMM_GMLDISCRETEMODEL_H 1

#include <xmlio/XMLIO_Element.h>
#include "ghmm++/GHMM_GMLAlphabet.h"
#include "ghmm++/GHMM_GMLState.h"
#include "ghmm++/GHMM_GMLTransition.h"
#include "ghmm++/GHMM_DiscreteModel.h"


#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif


class GHMM_GMLDiscreteModel;
class GHMM_Alphabet;

/** Discrete HMM model (wrapper around model in C data structure). */
class GHMM_GMLDiscreteModel: public GHMM_DiscreteModel {

 public:
  
  /** Constructor. */
  GHMM_GMLDiscreteModel();
  /** Constructor. 
      @param number_of_states Number of the states.
      @param alphabet_size Size of the alphabet. 
      @param prior Prior for the a priori probability for the model 
                   (-1 for none). */
  GHMM_GMLDiscreteModel(int number_of_states, int alphabet_size, double prior=-1);
  /** Constructor. Construct from c model. Object now is owner of this model. 
      @param my_model model as C data structure. */
  GHMM_GMLDiscreteModel(model* my_model);
  /** Constructor. */
  GHMM_GMLDiscreteModel(GHMM_GMLAlphabet* alphabet);

 /** Returns alphabet of model or NULL, if no such alphabet exists. */
  GHMM_GMLAlphabet* getAlphabet() const;


  /** */
  const int XMLIO_writeContent(XMLIO_Document& writer);

  /** Alphabet of model. */
  GHMM_GMLAlphabet* gmlalphabet;
 
 protected:
  
  /** */
  void createTransitions();

  /** */
  vector<GHMM_GMLTransition*> transitions;

  void init();
};


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_DISCRETEMODEL_H */
