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


GHMM_GMLClassWriter::GHMM_GMLClassWriter(GHMM_Alphabet *alphabet) {
  char tmp[3];
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "hmm:class";
  attributes["hmm:low"]   = "0";
  sprintf(tmp,"%d" , alphabet->size()-1);
  attributes["hmm:high"]  = tmp;

  for (int i = 0; i < alphabet->size(); ++i) {  
    addSymbol(alphabet->getSymbol(i));
  }
}

const int GHMM_GMLClassWriter::XMLIO_writeContent(XMLIO_Document& writer) {
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
