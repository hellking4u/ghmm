/*
  author: Achim Gaedke
  created: 18. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_graph.h
  $Id$
 */

#ifndef GHMMPP_GHMM_GRAPH_H
#define GHMMPP_GHMM_GRAPH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <gxl/gxl_graph.h>
#include <ghmm++/ghmm_edge.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     the graph of the ghmm model
   */
  class ghmm_graph: public gxl_graph
    {
    public:
      ghmm_graph();
      /** constructor from XMLIO interface */
      ghmm_graph(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);      
      /** */
      virtual const char* toString() const;
    };

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_GRAPH_H */


