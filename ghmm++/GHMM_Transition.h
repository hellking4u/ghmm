/*
  created: 19 Feb 2002 by Peter Pipenbacher
  authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
  file   : $Source$
  $Id$
  
__copyright__

*/

#ifndef _GHMM_TRANSITION_H
#define _GHMM_TRANSITION_H 1

#include <xmlio/XMLIO_Element.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Transition;

/** Represents transition between two states. Only needed while model is constructed from
    xml file. */
class GHMM_Transition: public XMLIO_Element {

 public:

  /** For internal use of reading xml files. */
  enum GHMM_TransitionReadingType {GHMM_TRANSITION_READING_NONE, GHMM_TRANSITION_READING_PROB};

  /** Constructor. */
  GHMM_Transition(XMLIO_Attributes &attrs);
  /** Destructor. */
  virtual ~GHMM_Transition();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** ID of source state. */
  string source;
  /** ID of source target. */
  string target;
  /** Probability of transition. */
  double prob;


 private:

  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  virtual void XMLIO_endTag(const string& tag);
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** Current reading state. */
  GHMM_TransitionReadingType reading;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_TRANSITION_H */
