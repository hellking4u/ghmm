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

#include <ppghmm++/begin_code.h>

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

  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** */
  virtual void XMLIO_endTag(const string& tag);

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

#include <ppghmm++/close_code.h>

#endif /* _GHMM_DOCUMENT_H */
