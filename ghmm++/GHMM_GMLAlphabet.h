/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#ifndef _GHMM_GMLALPHABET_H
#define _GHMM_GMLALPHABET_H 1


#include <ghmm++/GHMM_Alphabet.h>
#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GMLAlphabet;


/** */
class GHMM_GMLAlphabet: public GHMM_Alphabet {

 public:

  /** Constructor. */
  GHMM_GMLAlphabet();
  /** Destructor. */
  ~GHMM_GMLAlphabet();

  
 protected:

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  // const int XMLIO_writeContent(XMLIO_Document& doc);
};


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_ALPHABET_H */
