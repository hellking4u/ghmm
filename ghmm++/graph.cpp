/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/graph.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "graph.h"
#include "xmlio/XMLIO_IgnoreObject.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


Node::Node(const string& tag, XMLIO_Attributes &attributes)
{
  XMLIO_Attributes::const_iterator pos=attributes.find("id");
  if (pos!=attributes.end())
    {
      id=pos->second;
    }
  else
    {
      cerr<<toString()<<" id attribute missing!"<<endl;
    }
}

const char* Node::toString() const
{
  return "Node";
}

void Node::print() const
{
  cout<<toString()<<": id "<<id<<endl;
}

const string& Node::get_id() const
{
  return id;
}

Edge::Edge(const string& tag, XMLIO_Attributes &attributes)
{
  XMLIO_Attributes::const_iterator pos=attributes.find("id");
  /* looking for id */
  if (pos!=attributes.end())
    {
      id=pos->second;
    }
  else
    {
      cerr<<toString()<<" id attribute missing!"<<endl;
    }
  
  pos=attributes.find("from");
  if (pos!=attributes.end())
    {
      from=pos->second;
    }
  else
    {
      cerr<<toString()<<" from attribute missing!"<<endl;
    }

  pos=attributes.find("to");
  if (pos!=attributes.end())
    {
      to=pos->second;
    }
  else
    {
      cerr<<toString()<<" to attribute missing!"<<endl;
    }

}

const char* Edge::toString() const
{
  return "Egde";
}

void Edge::print() const
{
  cout<<toString()<<" from "<<from<<" to "<<to<<endl;
  for (vector<double>::const_iterator pos=class_weights.begin();pos!=class_weights.end();++pos)
    {
      cout<<"transition weight "<<*pos<<endl;
    }
}

const string& Edge::get_id() const
{
  return id;
}


Graph::Graph(const string& tag, XMLIO_Attributes &attributes)
{
  XMLIO_Attributes::const_iterator pos=attributes.find("id");
  if (pos!=attributes.end())
    {
      id=pos->second;
    }
  else
    {
      cerr<<toString()<<": id attribute missing!"<<endl;
    }
}

XMLIO_Object* Graph::XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes)
{
  if (tag=="Node")
    {
      Node* new_node=new Node(tag,attributes);
      vector<Node*>::push_back(new_node);
      return new_node;
    }
  else if (tag=="Edge")
    {
      Edge* new_edge=new Edge(tag,attributes);
      vector<Edge*>::push_back(new_edge);
      return new_edge;      
    }
  else
    {
      cerr<<toString()<<": found unexpected element "<<tag<<", ignoring"<<endl;
      return new XMLIO_IgnoreObject();
    }
}

void Graph::XMLIO_endTag(const string& tag)
{
  /* not needed */
}

void Graph::XMLIO_getCharacters(const string& characters)
{
  /* ignored */
}

void Graph::XMLIO_finishedReading()
{
  /* not needed now */
}

const char* Graph::toString() const
{
  return "Graph";
}

void Graph::print() const
{
  cout<<toString()<<" id: "<<id<<endl
      <<"Edges"<<endl;
  for (vector<Edge*>::const_iterator pos=vector<Edge*>::begin();pos!=vector<Edge*>::end();++pos)
    {
      if (*pos!=NULL)
	(*pos)->print();
      else
	cout<<"NULL-Element"<<endl;
    }
  cout<<"Nodes"<<endl;
  for (vector<Node*>::const_iterator pos=vector<Node*>::begin();pos!=vector<Node*>::end();++pos)
    {
      if (*pos!=NULL)
	(*pos)->print();
      else
	cout<<"NULL-Element"<<endl;
    }
}

const string& Graph::get_id() const
{
  return id;
}






