/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#ifndef _GHMM_EMISSION_H
#define _GHMM_EMISSION_H 1

#include <vector>
#include <xmlio/XMLIO_Element.h>
#include <ghmm/smodel.h>
#include <ghmm++/GHMM_Types.h>
#include <ghmm++/GHMM_State.h>
#include <ghmm++/GHMM_DiscreteModel.h>

#include <ghmm++/GHMM_StateT.hh>
#include <ghmm++/GHMM_SWDiscreteModel.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_DiscreteModel;
class GHMM_ContinuousModel;

/** */
template<class StateType, class ModelT> class GHMM_EmissionT: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_EmissionT(StateType* my_state)
    {
      state             = my_state;
      tag               = "emission";
      xmlio_indent_type = XMLIO_INDENT_BOTH;
      
      switch (getModelType()) {
	
      case GHMM_CONTINOUS:
	/* only read data if state has been initialized */
	if (my_state->c_sstate) {
	  //mue.push_back(my_state->c_sstate->mue[0]);
	  //variance.push_back(my_state->c_sstate->u[0]);
	  //weights.push_back(my_state->c_sstate->c[0]);
	  //density  = ((GHMM_ContinuousModel*)my_state->getModel())->c_model->density;
	}
	break;
	
      case GHMM_DISCRETE: break;
      }
    } // Constructor
  

  /** Returns name of class. */
  const char* toString() { return "GHMM_Emission";}

  /** Returns model type. */
  GHMM_ModelType getModelType(){ return state->getModelType(); }

  /** */
  virtual void XMLIO_finishedReading() { ; }
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters) { ; }
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) { return NULL; }
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc) { return 0; }

  /** */
  vector<double> mue;
  /** */
  vector<double> variance;
  /** */
  vector<double> weights;
  /** */
  density_t density;


 protected:

  /** Parent state. */
  StateType* state;
};

class GHMM_Emission : public GHMM_EmissionT<GHMM_State, GHMM_DiscreteModel>
{
 public:
  GHMM_Emission(GHMM_State* my_state) : GHMM_EmissionT<GHMM_State, GHMM_DiscreteModel>(my_state)
    {
      int i;
      switch (getModelType()) {
      case GHMM_DISCRETE:
	/* only read data if state has been initialized */
	if (my_state->c_state) {
	  for (i = 0; i < ((GHMM_DiscreteModel*)my_state->getModel())->c_model->M; ++i)
			weights.push_back(my_state->c_state->b[i]);
	}
	break;
      }
    }
  
  /** Destructor. */
  ~GHMM_Emission();
  
  /** Returns name of class. */
  const char* toString();
  
  /** */
  void XMLIO_finishedReading();
  /** Collects all character data. */
  void XMLIO_getCharacters(const string& characters);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  const int XMLIO_writeContent(XMLIO_Document& doc);
}; // GHMM_Emission
 
 
 class GHMM_GMLState;
 class GHMM_SWDiscreteModel;
 
class GHMM_GMLEmission : public GHMM_EmissionT<GHMM_GMLState, GHMM_SWDiscreteModel>
  {
 public:
    GHMM_GMLEmission(GHMM_GMLState* my_state) : GHMM_EmissionT<GHMM_GMLState, GHMM_SWDiscreteModel> (my_state)
      {
	int i;
	switch (getModelType()) {
	case GHMM_DISCRETE:
	  /* only read data if state has been initialized */
	  if (my_state->c_sdstate) {
	    for (i = 0; i < ((GHMM_SWDiscreteModel*)my_state->getModel())->c_model->M; ++i)
	      weights.push_back(my_state->c_state->b[i]);
	  }
	  break;
      }
      }
    
    /** Destructor. */
    ~GHMM_GMLEmission();
    
    /** Returns name of class. */
    const char* toString();
    
    /** */
    void XMLIO_finishedReading();
	/** Collects all character data. */
    void XMLIO_getCharacters(const string& characters);
    /** Called by GHMM_Document when a start tag is received. Tag and 
	attributes are passed to this function. */
    XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
    /** Writes the content (XML Spec[43]) of this element.
	You should use the public XMLIO_Document::write* functions.
	@return Returns the number of bytes written,
	but is negative when an error occured and 0 by default. */
    const int XMLIO_writeContent(XMLIO_Document& doc);
}; // GHMM_GMLEmission


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_EMISSION_H */
