/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_EMISSION_H
#define _GHMM_EMISSION_H 1

#include <xmlio/XMLIO_Element.h>
#include <ghmm/smodel.h>

#include <ppghmm++/begin_code.h>

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

  /** */
  virtual void XMLIO_finishedReading();
  /** Collects all character data. */
  virtual void XMLIO_getCharacters(const string& characters);
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** */
  double mue;
  /** */
  double variance;
  /** */
  double weight;
  /** */
  density_t density;


 private:

  /** */
  bool function_loaded;
  /** Parent state. */
  GHMM_State* state;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ppghmm++/close_code.h>

#endif /* _GHMM_EMISSION_H */
