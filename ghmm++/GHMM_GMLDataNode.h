/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#ifndef _GHMM_GMLDATANODE_H
#define _GHMM_GMLDATANODE_H 1

#include <xmlio/XMLIO_Element.h>
#include <xmlio/XMLIO_Document.h>
#include <ghmm++/GHMM_Types.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif


class GHMM_GMLState;

/** */
class GHMM_GMLDataNode: public XMLIO_Element {

 public:

  /** Constructor. */
  GHMM_GMLDataNode(GHMM_GMLState* my_state, const char *keyvalue);


  virtual const int XMLIO_writeContent(XMLIO_Document& writer);

 private:

  /** Parent state. */
  GHMM_GMLState* state;

  string keyvalue;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_EMISSION_H */
