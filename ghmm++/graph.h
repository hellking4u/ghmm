/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/graph.h
  $Id$
 */

#ifndef GHMMPP_GRAPH_H
#define GHMMPP_GRAPH_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xmlio/XMLIO_Object.h>
#include <ghmm++/sequences.h>
#include <ghmm++/graph.h>
#include <ghmm++/DiscretePD.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     Node used in Graph
   */
  class Node: public XMLIO_Object
    {};

  /**
     Edge used in Graph
   */
  class Edge: public XMLIO_Object
    {};

  /**
     Graph, conform to GXL, uses Edge and Node
   */
  class Graph:public XMLIO_Object, public vector<Node*>, public vector<Edge*>
    {};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GRAPH_H */




