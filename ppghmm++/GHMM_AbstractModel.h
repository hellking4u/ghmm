/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_ABSTRACTMODEL_H
#define _GHMM_ABSTRACTMODEL_H 1

#include <stdio.h>
#include <xmlio/XMLIO_Element.h>


#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_AbstractModel;

/** */
class GHMM_AbstractModel: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_AbstractModel();
  /** Destructor. */
  virtual ~GHMM_AbstractModel();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Writes the model in matrix format.
     @param file: output file
  */
  virtual void print(FILE *file);
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* _GHMM_ABSTRACTMODEL_H */
