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

#include <ppghmm++/begin_code.h>

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

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  virtual void XMLIO_finishedReading();
  /** */
  virtual void XMLIO_getCharacters(const string& characters);
  /** */
  virtual void XMLIO_endTag(const string& tag);
  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** C type state. Object is not owner of this state. */
  sstate* c_sstate;
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

#include <ppghmm++/close_code.h>

#endif /* _GHMM_STATE_H */
