/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#ifndef _GHMM_GRAPHMLDOC_H
#define _GHMM_GRAPHMLDOC_H 1

#include <xmlio/XMLIO_Document.h>

#include <ghmm++/GHMM_Document.h>
#include <ghmm++/GHMM_Alphabet.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_GraphMLDoc;


class GHMM_GraphMLDoc : public GHMM_Document {

 private:
  
  GHMM_Alphabet* tmp_alphabets;

 public:

  /** Constructor. */
  GHMM_GraphMLDoc();
  /** Destructor. */
  ~GHMM_GraphMLDoc();

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  void XMLIO_endTag(const string& tag);
  /** Writes the XML prolog (XML spec [22]).
      Only XMLDecl (XML specs [23]-[26]) is supported by calling XMLIO_writeXMLDeclaration()
      @return Returns nr of bytes written or an negative error code. */ 
  int XMLIO_writeProlog();
  /** Is called when a document is closed and writes an optional trailer,
      which must be of Misc-type (XML specs [27]) after the main element 
      (XML specs [1]).
      By default this writes a newline character.
      @return Returns nr of bytes written or an negative error code. */ 
  int XMLIO_writeTrailer();


};


#ifdef HAVE_NAMESPACES
}
#endif

#endif 
