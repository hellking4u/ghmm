/* @(#)GHMM_Alphabet.cpp created by Peter Pipenbacher at 19 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#include "GHMM_Alphabet.h"
#include "GHMM_Sequence.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Alphabet::GHMM_Alphabet() {
  alphabet_type     = GHMM_SINGLE_CHAR_ALPHABET;
  xmlio_indent_type = XMLIO_INDENT_BOTH;
  tag               = "alphabet";
}


GHMM_Alphabet::~GHMM_Alphabet() {
}


const char* GHMM_Alphabet::toString() const {
  return "GHMM_Alphabet";
}


void GHMM_Alphabet::addSymbol(const string& symbol) {
  switch (alphabet_type) {

  case GHMM_SINGLE_CHAR_ALPHABET:
    if (symbol.length() != 1) {
      fprintf(stderr,"Cannot add symbol '%s'.\nAlphabet only allow single character symbols.\n",symbol.c_str());
      exit(1);
    }
    symbols.push_back(symbol);
    symbol_map[symbol] = symbols.size() - 1;
    break;

  default: 
    fprintf(stderr,"Alphabet type does not support addSymbol() method.\n");
    exit(1);
  }
}


unsigned int GHMM_Alphabet::size() const {
  return symbols.size();
}


int GHMM_Alphabet::getIndex(const string& symbol) const {
  map<string,int>::const_iterator iter = symbol_map.find(symbol);

  if (iter != symbol_map.end())
    return iter->second;

  return -1;
}


GHMM_Sequence* GHMM_Alphabet::getSequence(const string& sequence) const {
  unsigned int i;
  GHMM_Sequence* seq = NULL;

  switch (alphabet_type) {

  case GHMM_SINGLE_CHAR_ALPHABET:
    seq = new GHMM_Sequence(GHMM_INT,sequence.length(),1);
    for (i = 0; i < sequence.length(); ++i)
      seq->c_i_sequences->seq[0][i] = getIndex(sequence.substr(i,1));

    break;

  default: 
    fprintf(stderr,"Alphabet type does not support getSequence() method.\n");
    exit(1);
  }

  return seq;
}


string GHMM_Alphabet::getSymbol(int index) const {
  return symbols[index];
}


const int GHMM_Alphabet::XMLIO_writeContent(XMLIO_Document& writer) {
  int total_bytes = 0;
  int result;
  unsigned int i;

  result = writer.writeEndl();

  if (result < 0)
    return result;
  total_bytes += result;

  writer.changeIndent(2);
  for (i = 0; i < symbols.size(); ++i) {
    result = writer.writef("%s<symbol code=\"%s\"/>\n",writer.indent,symbols[i].c_str());

    if (result < 0)
      return result;
    total_bytes += result;
  }

  return total_bytes;
}


XMLIO_Element* GHMM_Alphabet::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "symbol") // HACK!! We should look at the text node , not attribute
  { 
    addSymbol(attrs["code"]);
    return this;
  }

  return NULL;
}

