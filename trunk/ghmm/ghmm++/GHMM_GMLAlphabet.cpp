/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#include "GHMM_GMLAlphabet.h"
#include <xmlio/XMLIO_Element.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLAlphabet::GHMM_GMLAlphabet() {
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "hmm:alphabet";
}

GHMM_GMLAlphabet::~GHMM_GMLAlphabet() {

}


XMLIO_Element* GHMM_GMLAlphabet::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "map") 
    {
      return this;
    }
  else
    if (tag == "symbol") // HACK!! We should look at the text node , not attribute
      { 
	addSymbol(attrs["code"]);
	return this;
      }
  
  return NULL;
}


