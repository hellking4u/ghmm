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
#include "ghmm++/GHMM_GMLClass.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLClass::GHMM_GMLClass() {
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "hmm:class";
  reading = READ_NONE;
  attributes["hmm:type"] = "discrete"; 
  attributes["hmm:low"]  = "0";
  attributes["hmm:high"]  = "0";
}

GHMM_GMLClass::GHMM_GMLClass(GHMM_Alphabet *alp) {
  char tmp[10];
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "hmm:class";
  reading           = READ_NONE;
  attributes["hmm:type"] = "discrete"; 
  attributes["hmm:low"]  = "0";
  sprintf(tmp, "%2d", alp->size()-1); 
  attributes["hmm:high"]  = tmp;

  for (int i = 0; i < alp->size(); ++i) {  
    addSymbol(alp->getSymbol(i));
  }
}

GHMM_GMLClass::~GHMM_GMLClass() {

}


XMLIO_Element* GHMM_GMLClass::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
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

void GHMM_GMLClass::XMLIO_endTag(const string& tag) {
  reading = READ_NONE;
}


void GHMM_GMLClass::XMLIO_getCharacters(const string& characters) {

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


const int GHMM_GMLClass::XMLIO_writeContent(XMLIO_Document& writer) {
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


