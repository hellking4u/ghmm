/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_EMISSION_H
#define _GHMM_EMISSION_H 1

#include <vector>
#include <xmlio/XMLIO_Element.h>
#include <ghmm/smodel.h>
#include <ghmm++/GHMM_Types.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Emission;
class GHMM_State;

/** */
class GHMM_Emission: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_Emission(GHMM_State* my_state);
  /** Destructor. */
  virtual ~GHMM_Emission();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** Returns model type. */
  GHMM_ModelType getModelType() const;

  /** */
  virtual void XMLIO_finishedReading();
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc);

  /** */
  vector<double> mue;
  /** */
  vector<double> variance;
  /** */
  vector<double> weights;
  /** */
  density_t density;


 private:

  /** Parent state. */
  GHMM_State* state;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_EMISSION_H */
