/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_graph.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ghmm++/ghmm_graph.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm_graph::ghmm_graph():
  gxl_graph()
{}

ghmm_graph::ghmm_graph(const string& name, XMLIO_Attributes &attrs):
  gxl_graph(name,attrs)
{}

const char* ghmm_graph::toString() const
{
  return "ghmm_graph";
}


XMLIO_Element* ghmm_graph::XMLIO_startTag(const string& name, XMLIO_Attributes &attrs)
{
  if (name=="edge")
    {
      gxl_edge* new_edge=new ghmm_edge(name,attrs);
      vector<gxl_edge*>::push_back(new_edge);
      return new_edge;
    }
  else
    return gxl_graph::XMLIO_startTag(name,attrs);
}



