/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#ifndef _GHMM_GMLTRANSITION_H
#define _GHMM_GMLTRANSITION_H 1

#include <vector>
#include <ghmm++/GHMM_Transition.h>
#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GMLState;

template<class StateType> class GHMM_TransitionT : public GHMM_Transition {
 public:
  /** Constructor. */
	
  GHMM_TransitionT(XMLIO_Attributes &attrs) : GHMM_Transition(attrs)
    {;}
  
  /** Constructor. */
  GHMM_TransitionT(StateType* my_source, StateType* my_target) 
    {;}
  
  /** Destructor. */
  
  /** */
  vector<double> probs;
  
 protected:
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) = 0;
	
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc) = 0;
	
  virtual void XMLIO_getCharacters(const string& characters) = 0;
	
};

/** Represents transition between two states. Only needed while model is constructed from
    xml file. */
class GHMM_GMLTransition: public GHMM_TransitionT<GHMM_GMLState> {

 public:
  const char* toString() const;

  /** Constructor. */
  GHMM_GMLTransition(XMLIO_Attributes &attrs);
  
  GHMM_GMLTransition(GHMM_GMLState* my_source, GHMM_GMLState* my_target, vector<double> my_prob);
  
 private:

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  const int XMLIO_writeContent(XMLIO_Document& doc);

  void XMLIO_getCharacters(const string& characters);

};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_TRANSITION_H */

