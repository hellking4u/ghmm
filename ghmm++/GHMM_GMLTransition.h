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

#include <ghmm++/GHMM_Transition.h>
#include <ghmm++/GHMM_State.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GMLTransition;

/** Represents transition between two states. Only needed while model is constructed from
    xml file. */
class GHMM_GMLTransition: public GHMM_Transition {

 public:

  /** Constructor. */
  GHMM_GMLTransition(XMLIO_Attributes &attrs);
  /** Constructor. */
  GHMM_GMLTransition(GHMM_State* my_source, GHMM_State* my_target, double my_prob);
  /** Destructor. */
  /** */
  const char* toString() const;

 protected:

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  const int XMLIO_writeContent(XMLIO_Document& doc);

};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_TRANSITION_H */

