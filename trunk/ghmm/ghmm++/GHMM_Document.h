/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#ifndef _GHMM_DOCUMENT_H
#define _GHMM_DOCUMENT_H 1

#include <xmlio/XMLIO_Document.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Document;
class GHMM_DiscreteModel;
class GHMM_ContinuousModel;
class GHMM_Sequences;

/** */
class GHMM_Document: public XMLIO_Document {

 public:

  /** Constructor. */
  GHMM_Document();
  /** Destructor. */
  virtual ~GHMM_Document();

  /** Returns continuous model, which has been read from file or NULL
      if no such model exists. */
  GHMM_ContinuousModel* getContinuousModel() const;
  /** Returns discrete model, which has been read from file or NULL
      if no such model exists. */
  GHMM_DiscreteModel* getDiscreteModel() const;
  /** Returns sequences, which has been read from file or NULL
      if no such model exists. */
  GHMM_Sequences* getSequences() const;

  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Called by XMLIO_Document when a end tag is found. 
      This happens when a sub element has finished reading its
      content. By default this function does nothing. */
  virtual void XMLIO_endTag(const string& tag);
  /** Writes the XML prolog (XML spec [22]).
      Only XMLDecl (XML specs [23]-[26]) is supported by calling XMLIO_writeXMLDeclaration()
      @return Returns nr of bytes written or an negative error code. */ 
  virtual int XMLIO_writeProlog();
  /** Is called when a document is closed and writes an optional trailer,
      which must be of Misc-type (XML specs [27]) after the main element 
      (XML specs [1]).
      By default this writes a newline character.
      @return Returns nr of bytes written or an negative error code. */ 
  virtual int XMLIO_writeTrailer();

  /** Returns name of class. */
  virtual const char* toString() const;


 private:

  /** */
  GHMM_DiscreteModel* discrete_model;
  /** */
  GHMM_ContinuousModel* continuous_model;
  /** */
  GHMM_Sequences* sequences;
  /** */
  bool reading_ghmm;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_DOCUMENT_H */
