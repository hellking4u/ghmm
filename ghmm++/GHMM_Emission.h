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

#include <assert.h>
#include <vector>
#include <xmlio/XMLIO_Element.h>
#include <ghmm/smodel.h>
#include <ghmm++/GHMM_Types.h>
#include <ghmm++/GHMM_State.h>
#include <ghmm++/GHMM_DiscreteModel.h>
#include <ghmm++/GHMM_ContinuousModel.h>
#include <ghmm++/GHMM_StateT.hh>
#include <ghmm++/GHMM_SWDiscreteModel.h>
#include <iostream>
#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_DiscreteModel;
class GHMM_ContinuousModel;

 template<class StateType, class ModelType> class GHMM_CEmissionT: public XMLIO_Element {
 public:
   /** Constructor. */
   GHMM_CEmissionT(StateType* my_state)
     {
       state             = my_state;
       tag               = "emission";
       xmlio_indent_type = XMLIO_INDENT_BOTH;
     }

   /** Returns name of class. */
  const char* toString() { return "GHMM_CEmissionT";}
  
  /** */
  vector<double> mue;
  /** */
  vector<double> variance;
  /** */
  vector<double> weights;
  /** */
  density_t density;

 protected:

  /** */
  virtual void XMLIO_finishedReading() = 0;
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters) = 0;
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) = 0;
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc) = 0;

  /** Parent state. */
  StateType* state;
 };

/** */
 template<class StateType, class ModelT> class GHMM_DEmissionT: public XMLIO_Element {
   
 public:
   
   /** Constructor. */
   GHMM_DEmissionT(StateType* my_state)
     {
       state             = my_state;
       tag               = "emission";
       xmlio_indent_type = XMLIO_INDENT_BOTH;      
     }
   
   /** Returns name of class. */
   const char* toString() { return "GHMM_DEmissionT";}
  
   /** Returns model type. */
   GHMM_ModelType getModelType(){ return state->getModelType(); }

   /** */
   vector<double> weights;

   /** Collects all character data. */
   void XMLIO_getCharacters(const string& characters) {
     unsigned int pos;
     for (pos = 0; pos < characters.size(); ++pos) {
       while (pos < characters.size() && isspace(characters[pos]))
	 ++pos;
       
       if (pos < characters.size())
	 weights.push_back(atof(characters.substr(pos).c_str()));
       
       while (pos < characters.size() && !isspace(characters[pos]))
	 ++pos;
     }
   }

   /** Writes the content (XML Spec[43]) of this element.
       You should use the public XMLIO_Document::write* functions.
       @return Returns the number of bytes written,
       but is negative when an error occured and 0 by default. */
   virtual const int XMLIO_writeContent(XMLIO_Document& doc) = 0;

 protected:

   /** Parent state. */
   StateType* state;
 };

 class GHMM_DEmission : public GHMM_DEmissionT<GHMM_State, GHMM_DiscreteModel>
   {
   public:
     GHMM_DEmission(GHMM_State* my_state) : GHMM_DEmissionT<GHMM_State, GHMM_DiscreteModel>(my_state)
       {
	 int i;
	 /* only read data if state has been initialized */
	 if (my_state->c_state) {
	   for (i = 0; i < ((GHMM_DiscreteModel*)my_state->getModel())->c_model->M; ++i)
	     weights.push_back(my_state->c_state->b[i]);
	 }
       }
     
     
     /** Returns name of class. */
     const char* toString() { return "GHMM_DEmission"; }
     
     /** Writes the content (XML Spec[43]) of this element.
	 You should use the public XMLIO_Document::write* functions.
	 @return Returns the number of bytes written,
	 but is negative when an error occured and 0 by default. */
     const int XMLIO_writeContent(XMLIO_Document& doc);
  
   }; // GHMM_Emission
 

class GHMM_CEmission : public GHMM_CEmissionT<GHMM_State, GHMM_ContinuousModel>
  {
  public:

    GHMM_CEmission(GHMM_State* my_state) : GHMM_CEmissionT<GHMM_State, GHMM_ContinuousModel>(my_state) {
      if (my_state->c_sstate) {  // only read data if state has been initialized 
	mue.push_back(my_state->c_sstate->mue[0]);
	variance.push_back(my_state->c_sstate->u[0]);
	weights.push_back(my_state->c_sstate->c[0]);
	// problem with getModel(), using GMLState 
	assert((smodel *)(my_state->getModel()->get_cmodel()));
	density  = ((smodel *)(my_state->getModel()->get_cmodel()))->density; 
      }
    }
    
    /** Returns name of class. */
    const char* toString() { return "GHMM_CEmission"; }
    
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
 
class GHMM_GMLEmission : public GHMM_DEmissionT<GHMM_GMLState, GHMM_SWDiscreteModel>
  {
 public:
    GHMM_GMLEmission(GHMM_GMLState* my_state) : GHMM_DEmissionT<GHMM_GMLState, GHMM_SWDiscreteModel> (my_state)
      {
	int i;
	/* only read data if state has been initialized */
	if (my_state->c_sdstate) {
	  for (i = 0; i < ((GHMM_SWDiscreteModel*)my_state->getModel())->c_model->M; ++i)
	    weights.push_back(my_state->c_sdstate->b[i]);
	}
      }
    
    /** Destructor. */
    ~GHMM_GMLEmission();
    
    /** Returns name of class. */
    const char* toString();
    

    const int XMLIO_writeContent(XMLIO_Document& doc);
}; // GHMM_GMLEmission


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_EMISSION_H */
