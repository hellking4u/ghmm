#ifndef _GHMM_STATET_H
#define _GHMM_STATET_H 1


#include <map>
#include <ghmm/sdmodel.h>

#include <ghmm++/GHMM_Toolkit.h>
#include <xmlio/XMLIO_Element.h>
#include <ghmm++/GHMM_Types.h>


#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_ContinuousModel;
class GHMM_GMLAbstractModelT;
class GHMM_Alphabet;

/** */
template<class CStateType, class TransitionT, class ModelT> class GHMM_StateT: public XMLIO_Element {

public:

  /** */
  enum GHMM_StateReadingType {GHMM_STATE_NONE,GHMM_STATE_INITIAL, GHMM_STATE_LABEL, GHMM_STATE_COUNTME};

  /** Constructor */
  GHMM_StateT(ModelT* my_model, int my_index, XMLIO_Attributes& attrs)
  {
	index             = my_index;
	c_sdstate         = NULL;
	reading           = GHMM_STATE_NONE;
	parent_model      = my_model;
	tag               = "state";
	xmlio_indent_type = XMLIO_INDENT_BOTH;

	/* by default take index as id. */
	id = attrs["id"];
	if (id == "")
		id = GHMM_Toolkit::toString(my_index);

	my_model->addStateID(id,index);
	attributes = attrs;
  }
  /** Constructor */
  GHMM_StateT(ModelT* my_model, int my_index, CStateType* my_state)
  {
    index             = my_index;
	c_sdstate         = my_state;
	reading           = GHMM_STATE_NONE;
	// emission          = NULL;
	parent_model      = my_model;
	tag               = "state";
	xmlio_indent_type = XMLIO_INDENT_BOTH;
  
	/* take index as id. */
	id = GHMM_Toolkit::toString(my_index);
	my_model->addStateID(id,index);
  }

  GHMM_StateT* getMySelf() { return this; }
  
  /** Destructor. */
  ~GHMM_StateT() {}

  /** */
  void changeOutEdge(int target, double prob);
  /** */
  void changeOutEdge(int matrix_index, int target, double prob);
  /** */
  void changeInEdge(int source, double prob);
  /** */
  void changeInEdge(int matrix_index, int source, double prob);
  
  /** Creates GHMM_Transition object of outgoing edge with given index. */
  virtual TransitionT* createTransition(int edge_index) { return NULL; }

  /** Fills given state. */
  virtual void fillState(CStateType* s) = 0;

  /** */
  void setID(const string& my_id);
  /** */
  ModelT* getModel()            const { return parent_model; }
  /** Returns model type. */
  GHMM_ModelType getModelType() const { return getModel()->getModelType(); }
  /** */
  int getOutEdges()   const { return c_sdstate->out_states; }

  /** Returns initial probability of state. */
  double getInitial() const { return c_sdstate->pi; }

  /** Sets initial probability of this state to 'prob'. */
  void setInitialProbability(double prob) { c_sdstate->pi = prob; }

  /** Sets output probability of symbol 'index' to 'prob'. */
  virtual void setOutputProbability(int index, double prob) = 0;

  /** Sets output probability of symbol to 'prob'. */
  virtual void setOutputProbability(const string& symbol, double prob) = 0;

  /** Returns name of class. */
  virtual const char* toString() const = 0;
  
  /** C type sstate. Object is not owner of this state. */
  CStateType *c_sdstate; 

  double getEmissionFrom(int i) {
    if (c_sdstate) {
      return (double) c_sdstate->b[i];
    } else {
      return -1;
    }
  }

  /** */
  string id;

  int index;


 protected:

  /** */
  void removeInEdge(int source);
  /** */
  void removeOutEdge(int target);

  /** */
  virtual void XMLIO_finishedReading() {}
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters) {}
  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  virtual void XMLIO_endTag(const string& tag) {}
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {}
  /** Returns the attributes of this element (XML Spec [40], [41]). 
      By default returns content of variable 'attributes'.*/
  virtual XMLIO_Attributes& XMLIO_getAttributes() {}
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc) {}

  /** */
  GHMM_StateReadingType reading;
  /** */
  double initial;

  /** */
  ModelT* parent_model;
};


class GHMM_GMLEmission;
class GHMM_GMLTransition;
class GHMM_SWDiscreteModel;


/** */
class GHMM_GMLState : public GHMM_StateT<sdstate, GHMM_GMLTransition, GHMM_SWDiscreteModel> {

 public:
  /** */
  GHMM_GMLState(GHMM_SWDiscreteModel* my_model, int index, XMLIO_Attributes& attrs);
  /** */
  //GHMM_GMLState(GHMM_AbstractModel* my_model, int index, sstate* my_state);
  /** */
  //GHMM_GMLState(GHMM_AbstractModel* my_model, int index, state* my_state);
  /** */
  GHMM_GMLState(GHMM_SWDiscreteModel* my_model, int index, sdstate* my_state);

  /** */
  ~GHMM_GMLState();

  /** Returns name of class. */
  const char* toString() const;

  void fillState(sdstate* s);

  GHMM_GMLTransition *createTransition(int edge_index, sdmodel *mo);

  GHMM_GMLState* getMySelf() { return this; }

  float get2DPosition(int index);

  /** Sets output probability of symbol 'index' to 'prob'. */
  void setOutputProbability(int index, double prob);

  /** Sets output probability of symbol to 'prob'. */
  void setOutputProbability(const string& symbol, double prob);

  /** */
  string label;
  int m_countme;

  GHMM_GMLEmission* emission;

 protected:

  /** Collects all character data. */
  void XMLIO_getCharacters(const string& characters);
  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  void XMLIO_endTag(const string& tag);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Returns the attributes of this element (XML Spec [40], [41]). 
      By default returns content of variable 'attributes'.*/
  XMLIO_Attributes& XMLIO_getAttributes();
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  const int XMLIO_writeContent(XMLIO_Document& doc);

  /** Position attributes */
  float vPosition[3];


 private:
  map<string, int> hasData;
};


#ifdef HAVE_NAMESPACES
}
#endif

#endif //
