/* @(#)GHMM_Alphabet.h created by Peter Pipenbacher at 19 Mar 2002
 *
 * Authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 *
 */

#ifndef _GHMM_ALPHABET_H
#define _GHMM_ALPHABET_H 1

#include <string>
#include <vector>
#include <map>
#include <xmlio/XMLIO_Element.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Alphabet;
class GHMM_Sequence;


/** */
class GHMM_Alphabet: public XMLIO_Element {

 public:

  /** */
  enum GHMM_AlphabetType {GHMM_SINGLE_CHAR_ALPHABET};

  /** Constructor. */
  GHMM_Alphabet();
  /** Destructor. */
  virtual ~GHMM_Alphabet();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  void addSymbol(const string& symbol);
  /** */
  int getIndex(const string& symbol) const;
  /** */
  GHMM_Sequence* getSequence(const string& sequence) const;
  /** */
  string getSymbol(int index) const;
  /** */
  unsigned int size() const;

  
 private:

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc);

  /** */
  vector<string> symbols;
  /** */
  map<string,int> symbol_map;
  /** */
  GHMM_AlphabetType alphabet_type;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_ALPHABET_H */
