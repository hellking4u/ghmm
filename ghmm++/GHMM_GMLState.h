/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#ifndef _GHMM_GMLSTATE_H
#define _GHMM_GMLSTATE_H 1


#include <ghmm++/GHMM_Types.h>
#include <ghmm++/GHMM_State.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GMLState;

/** */
class GHMM_GMLState: public GHMM_State {

 public:
  /** */
  GHMM_GMLState(GHMM_AbstractModel* my_model, int index, XMLIO_Attributes& attrs);
  /** */
  GHMM_GMLState(GHMM_AbstractModel* my_model, int index, sstate* my_state);
  /** */
  GHMM_GMLState(GHMM_AbstractModel* my_model, int index, state* my_state);

  /** Returns name of class. */
  const char* toString() const;

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
};


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_STATE_H */
