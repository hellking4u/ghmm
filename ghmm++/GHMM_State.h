/*
 * created: 21 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#ifndef _GHMM_STATE_H
#define _GHMM_STATE_H 1

#include <ghmm/model.h>
#include <xmlio/XMLIO_Element.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_ContinuousModel;
class GHMM_Emission;
class GHMM_AbstractModel;

/** */
class GHMM_State: public XMLIO_Element {

 public:

  /** */
  enum GHMM_StateReadingType {GHMM_STATE_NONE,GHMM_STATE_INITIAL};

  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, XMLIO_Attributes& attrs);
  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, sstate* my_state);
  /** */
  GHMM_State(GHMM_AbstractModel* my_model, int index, state* my_state);
  /** Destructor. */
  virtual ~GHMM_State();

  /** */
  void changeOutEdge(int target, double prob);
  /** */
  void changeOutEdge(int matrix_index, int target, double prob);
  /** */
  void changeInEdge(int source, double prob);
  /** */
  void changeInEdge(int matrix_index, int source, double prob);
  /** Fills given state. */
  void fillState(sstate* s);
  /** Fills given state. */
  void fillState(state* s);
  /** Sets initial probability of this state to 'prob'. */
  void setInitialProbability(double prob);
  /** Sets output probability of symbol 'index' to 'prob'.*/
  void setOutputProbability(int index, double prob);

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  virtual void XMLIO_finishedReading();
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters);
  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  virtual void XMLIO_endTag(const string& tag);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** C type sstate. Object is not owner of this state. */
  sstate* c_sstate;
  /** C type state. Object is not owner of this state. */
  state* c_state;
  /** */
  string id;
  /** */
  GHMM_Emission* emission;
  /** */
  int index;


 private:

  /** */
  void removeInEdge(int source);
  /** */
  void removeOutEdge(int target);

  /** */
  GHMM_AbstractModel* parent_model;
  /** */
  GHMM_StateReadingType reading;
  /** */
  double initial;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_STATE_H */
