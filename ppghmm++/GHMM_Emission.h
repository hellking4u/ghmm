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


#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Emission;

/** */
class GHMM_Emission: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_Emission();
  /** Destructor. */
  virtual ~GHMM_Emission();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  void XMLIO_getCharacters(const string& characters);
  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** */
  double mue;
  /** */
  double variance;
  /** */
  double weight;
  /** */
  density_t density;
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* _GHMM_EMISSION_H */
