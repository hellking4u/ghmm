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


#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_State;
class GHMM_ContinuousModel;

/** */
class GHMM_State: public XMLIO_Element {

 public:

  /** */
  enum GHMM_StateReadingType {GHMM_STATE_NONE,GHMM_STATE_INITIAL};

  /** */
  GHMM_State(XMLIO_Attributes& attrs);
  /** Destructor. */
  virtual ~GHMM_State();

  /** */
  void fillState(GHMM_ContinuousModel* model, sstate* s);

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  virtual void XMLIO_getCharacters(const string& characters);
  /** */
  virtual void XMLIO_endTag(const string& tag);
  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** C type state. */
  state* c_state;
  /** */
  string id;


 private:

  /** */
  GHMM_StateReadingType reading;
  /** */
  int emission_weight;
  /** */
  double emission_mue;
  /** */
  double emission_variance;
  /** */
  double initial;
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* _GHMM_STATE_H */
