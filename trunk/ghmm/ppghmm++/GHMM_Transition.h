/*
 * created: 19 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_TRANSITION_H
#define _GHMM_TRANSITION_H 1

#include <xmlio/XMLIO_Element.h>


#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Transition;

/** */
class GHMM_Transition: public XMLIO_Element {

 public:

  /** */
  enum GHMM_TransitionReadingType {GHMM_TRANSITION_NONE, GHMM_TRANSITION_PROB};

  /** Constructor. */
  GHMM_Transition(XMLIO_Attributes &attrs);
  /** Destructor. */
  virtual ~GHMM_Transition();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  virtual void XMLIO_endTag(const string& tag);
  /** */
  virtual void XMLIO_getCharacters(const string& characters);
  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** */
  string source;
  /** */
  string target;
  /** */
  double prob;


 private:

  GHMM_TransitionReadingType reading;
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* _GHMM_TRANSITION_H */
