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

#include <set>
#include <map>
#include <vector>
#include <xmlio/XMLIO_Object.h>
#include <xmlio/XMLIO_ArrayObject.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     Node used in Graph
   */
  class Node: public XMLIO_Object
    {
    public:
      /** */
      Node(const string& tag, XMLIO_Attributes &attributes);
      /** Returns name of class. */
      virtual const char* toString() const;
      /** dumps all content to cout */
      virtual void print() const;
      /**
	 return id of this model
       */
      virtual const string& get_id() const;
    private:
      string id;
    };

  /**
     Edge used in Graph
   */
  class Edge: public XMLIO_Object
    {
    public:
      /** */
      Edge(const string& tag, XMLIO_Attributes &attributes);
      /** Returns name of class. */
      virtual const char* toString() const;
      /** dumps all content to cout */
      virtual void print() const;
      /**
	 return id of this edge
       */
      virtual const string& get_id() const;
      /**
	 return id of this edges origin
       */
      virtual const string& get_from_id() const;
      /**
	 return id of this edges destination
       */
      virtual const string& get_to_id() const;

    private:
      string id;
      string from;
      string to;
      vector<double> class_weights;
    };

  /**
     to do: not correct
   */
  class TransitionClasses: public XMLIO_ArrayObject<double>
    {};

  /**
     Graph, conform to GXL, uses Edge and Node
   */
  class Graph:public XMLIO_Object, public vector<Node*>, public vector<Edge*>
    {
    public:
      /** */
      Graph(const string& tag, XMLIO_Attributes &attributes);
      /**
	 expected elements are: Graph, Initial, Emissions and Sequences
       */
      virtual XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual void XMLIO_endTag(const string& tag);
      /** */
      virtual void XMLIO_getCharacters(const string& characters);
      /** */
      virtual void XMLIO_finishedReading();
      /** Returns name of class. */
      virtual const char* toString() const;
      /** dumps all content to cout */
      virtual void print() const;
      /**
	 return id of this model
       */
      virtual const string& get_id() const;
    private:
      /** id is mandatory for gxl graphs */
      string id;

      /** map from edge string-id  to its index */
      map<string,int> edge_idx;

      /** map from node string-id to its index */
      map<string,int> node_idx;

      /** adiascence list */
      map<int,set<int> > from_to_map;
      /** inverse adiascence list */
      map<int,set<int> > to_from_map;

    };

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GRAPH_H */




