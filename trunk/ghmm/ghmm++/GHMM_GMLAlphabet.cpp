/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <iostream>
#include <xmlio/XMLIO_Document.h>
#include "ghmm++/GHMM_GMLAlphabet.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLAlphabet::GHMM_GMLAlphabet() {
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "hmm:alphabet";
  reading = READ_NONE;
  attributes["hmm:type"] = "discrete"; 
  attributes["hmm:low"]  = "0";
  attributes["hmm:high"]  = "0";
}


GHMM_GMLAlphabet::~GHMM_GMLAlphabet() {

}


XMLIO_Element* GHMM_GMLAlphabet::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  reading = READ_NONE;
  if (tag == "map") 
    {
      return this;
    }
  else
    if (tag == "symbol") // HACK!! We should look at the text node , not attribute
      { 
	reading = READ_SINGLE_CHAR;
	return this;
      }
  
  return NULL;
}

void GHMM_GMLAlphabet::XMLIO_endTag(const string& tag) {
  reading = READ_NONE;
}


void GHMM_GMLAlphabet::XMLIO_getCharacters(const string& characters) {

  switch (reading) {

  case READ_SINGLE_CHAR:
    cout << "Alphabet seen: " << characters.c_str() << endl;
    addSymbol(characters);
    char tmp[3];
    sprintf(tmp, "%d", size()-1);
    attributes["hmm:high"]  = tmp;
    break;

  case READ_NONE:
    break;
  }
}


const int GHMM_GMLAlphabet::XMLIO_writeContent(XMLIO_Document& writer) {
  int total_bytes = 0;
  int result;
  unsigned int i;

  result = writer.writeEndl();

  if (result < 0)
    return result;
  total_bytes += result;

  writer.changeIndent(2);

  result = writer.writef("%s<map>\n",writer.indent);
  total_bytes += result;

  for (i = 0; i < symbols.size(); ++i) {
    result = writer.writef("%s<symbol code=\"%2d\">%s</symbol>\n",
			   writer.indent, i, symbols[i].c_str());

    if (result < 0)
      return result;
    total_bytes += result;
  }

  result = writer.writef("%s</map>\n",writer.indent);
  total_bytes += result;

  return total_bytes;
}


