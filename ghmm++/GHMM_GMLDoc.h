/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#ifndef _GHMM_GMLDOC_H
#define _GHMM_GMLDOC_H 1

#include <xmlio/XMLIO_Document.h>

#include <ghmm++/GHMM_GMLAlphabet.h>
#include <ghmm++/GHMM_GMLClass.h>
#include <ghmm++/GHMM_ContinuousModel.h>
#include <ghmm++/GHMM_SWDiscreteModel.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif


class GHMM_GraphMLDoc : public XMLIO_Document {


 public: 

  enum enumModelType { NONE, GHMM_DISCRETE, GHMM_CONTINUOUS };

  /** Constructor. */
  GHMM_GraphMLDoc(enumModelType model_type=GHMM_DISCRETE);
  /** Destructor. */
  ~GHMM_GraphMLDoc();

  /** Returns continuous model, which has been read from file or NULL
      if no such model exists. */
  GHMM_ContinuousModel* getContinuousModel() const;
  /** Returns discrete model, which has been read from file or NULL
      if no such model exists. */
  GHMM_SWDiscreteModel* getDiscreteModel() const;
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

  GHMM_Alphabet *getClass() { return (GHMM_Alphabet*) hmmclass; }

 protected:

  /** */
  // GHMM_GMLDiscreteModel*   discrete_model;
  /** */
  GHMM_ContinuousModel* continuous_model;
  /** */
  GHMM_SWDiscreteModel*   sdiscrete_model;
  /** */
  GHMM_Sequences* sequences;
  /** */
  bool reading_ghmm;
  /** */
  enumModelType model_type;

 private:
  GHMM_GMLClass*    hmmclass;
  GHMM_GMLAlphabet* tmp_alphabets;

};


#ifdef HAVE_NAMESPACES
}
#endif

#endif 