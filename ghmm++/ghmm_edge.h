/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_edge.h
  $Id$
 */

#ifndef GHMMPP_GHMM_EDGE_H
#define GHMMPP_GHMM_EDGE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xmlio/XMLIO_Object.h>
#include <gxl/gxl_edge.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     gxl_edge is extended for transition classes
   */
  class ghmm_edge: public gxl_edge, public vector<double>
    {
    public:
      ghmm_edge();
      /** constructor from XMLIO interface */
      ghmm_edge(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual const char* toString() const;
      /** */
      virtual void XMLIO_endTag(const string& tag);
      /** */
      virtual void print() const;
    };

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_H */


