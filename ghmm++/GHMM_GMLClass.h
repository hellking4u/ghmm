/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#ifndef _GHMM_GMLCLASS_H
#define _GHMM_GMLCLASS_H 1

#include <xmlio/XMLIO_Element.h>
#include <ghmm++/GHMM_Alphabet.h>
#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GMLClass;


/** */
class GHMM_GMLClass: public GHMM_Alphabet {

 public:

  enum READ_CHAR_TYPE { READ_SINGLE_CHAR, READ_NONE };

  /** Constructor. */
  GHMM_GMLClass();
  GHMM_GMLClass(GHMM_Alphabet *alp);  
  /** Destructor. */
  ~GHMM_GMLClass();
  const int XMLIO_writeContent(XMLIO_Document& writer);

 protected:

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  void XMLIO_endTag(const string& tag);

  void XMLIO_getCharacters(const string& characters);

 private:
  READ_CHAR_TYPE reading;

};


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_ALPHABET_H */
