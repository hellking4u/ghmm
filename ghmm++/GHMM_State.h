/*
 * created: 21 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 * 
 */

#ifndef _GHMM_STATE_H
#define _GHMM_STATE_H 1


#include <ghmm/model.h>
#include <ghmm/sdmodel.h>
#include <xmlio/XMLIO_Element.h>
#include <ghmm++/GHMM_Types.h>
#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_ContinuousModel;
class GHMM_CEmission;
class GHMM_DEmission;
class GHMM_AbstractModel;
class GHMM_Transition;
class GHMM_Alphabet;

/** */
class GHMM_State: public XMLIO_Element {

 public:

  /** */
  enum GHMM_StateReadingType {GHMM_STATE_NONE,GHMM_STATE_INITIAL, GHMM_STATE_LABEL, GHMM_STATE_COUNTME};

  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, XMLIO_Attributes& attrs);
  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, sstate* my_state);
  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, state* my_state);
  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, sdstate* my_state);
  /** Destructor. */
  virtual ~GHMM_State();

  GHMM_State* getMySelf() { return this; }

  /** */
  void changeOutEdge(int target, double prob);
  /** */
  void changeOutEdge(int matrix_index, int target, double prob);
  /** */
  void changeInEdge(int source, double prob);
  /** */
  void changeInEdge(int matrix_index, int source, double prob);
  /** Creates GHMM_Transition object of outgoing edge with given index. */
  GHMM_Transition* createTransition(int edge_index);
  /** Fills given state. */
  void fillState(sstate* s);
  /** Fills given state. */
  void fillState(state* s);
  /** Fills given state. */
  void fillState(sdstate* s);
  /** */
  void setID(const string& my_id);
  /** */
  GHMM_AbstractModel* getModel() const;
  /** Returns model type. */
  GHMM_ModelType getModelType() const;
  /** */
  int getOutEdges() const;
  /** Returns initial probability of state. */
  double getInitial() const;
  /** Sets initial probability of this state to 'prob'. */
  void setInitialProbability(double prob);
  /** Sets output probability of symbol 'index' to 'prob'. */
  void setOutputProbability(int index, double prob);
  /** Sets output probability of symbol to 'prob'. */
  void setOutputProbability(const string& symbol, double prob);

  /** Returns name of class. */
  virtual const char* toString() const;

  /** C type sstate. Object is not owner of this state. */
  sstate* c_sstate;
  /** C type state. Object is not owner of this state. */
  state* c_state;
  /** C type state. Object is not owner of this state. */
  sdstate* c_sdstate;
  /** */
  string id;
  /** */
  GHMM_CEmission* cemission;
  GHMM_DEmission* demission;
  /** */
  int index;


 protected:

  /** */
  void removeInEdge(int source);
  /** */
  void removeOutEdge(int target);

  /** */
  void XMLIO_finishedReading();
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

  /** */
  GHMM_StateReadingType reading;
  /** */
  double initial;

  /** */
  GHMM_AbstractModel* parent_model;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_STATE_H */
