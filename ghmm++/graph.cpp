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
#include <iterator>
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

const string& Edge::get_to_id() const
{
  return to;
}

const string& Edge::get_from_id() const
{
  return from;
}


Graph::Graph(const string& tag, XMLIO_Attributes &attributes)
{
  edge_id_counter=0;
  node_id_counter=0;
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
  /* should add element to several internal structures... */
  /* add nodes to node_id map */
  if (tag=="Node")
    {
      Node* last_node=vector<Node*>::back();
      /* test whether node id is unique */
      if (node_ids.find(last_node->get_id())==node_ids.end())
	{
	  cerr<<toString()<<"double node id "<<last_node->get_id()<<" found"<<endl;
	}
      else
	{
	  node_ids.insert(pair<const string,int>(last_node->get_id(),node_id_counter));
	  node_id_counter++;
	}
      
    }
  /* add edges to edge_id map */
  else if (tag=="Edge")
    {
      Edge* last_edge=vector<Edge*>::back();
      /* test whether node id is unique */
      if (edge_ids.find(last_edge->get_id())==edge_ids.end())
	{
	  cerr<<toString()<<"double edge id "<<last_edge->get_id()<<" found"<<endl;
	}
      else
	{
	  edge_ids.insert(pair<const string,int>(last_edge->get_id(),edge_id_counter));
	  edge_id_counter++;
	}
    }
}

void Graph::XMLIO_getCharacters(const string& characters)
{
  /* ignored */
}

void Graph::XMLIO_finishedReading()
{
  /* sanity check */

  /* collect transitions, indexed by order in vector */
  vector<Edge*>::iterator edge_pos=vector<Edge*>::begin();
  while (edge_pos!=vector<Edge*>::end())
    {
      const string& from=(*edge_pos)->get_from_id();
      const string& to=(*edge_pos)->get_to_id();
      /* exist nodes in map ? what is the id? */
      int from_id;
      map<string,int>::iterator from_node_pos=node_ids.find(from);
      if (from_node_pos==node_ids.end())
	{
	  cerr<<"found unknown from-node id during transition collection"<<endl;
	  node_ids.insert(pair<const string,int>(from,node_id_counter));
	  from_id=node_id_counter;
	  node_id_counter++;
	}
      else
	{
	  from_id=from_node_pos->second;
	}
      int to_id;
      map<string,int>::iterator to_node_pos=node_ids.find(to);
      if (to_node_pos==node_ids.end())
	{
	  cerr<<"found unknown from-node id during transition collection"<<endl;
	  node_ids.insert(pair<const string,int>(to,node_id_counter));
	  to_id=node_id_counter;
	  node_id_counter++;
	}
      else
	{
	  to_id=to_node_pos->second;
	}

      /* more than one edge between two points allowed ! */

      /* adiascence list */
      map<int,set<int> >::iterator from_to_map_pos=from_to_map.find(from_id);
      if (from_to_map_pos==from_to_map.end())
	{
	  from_to_map.insert(pair<const int,set<int> >(from_id,set<int>()));
	  from_to_map_pos=from_to_map.find(from_id);
	}
      /* overkill ... need index of vector<Edge*> */
      from_to_map_pos->second.insert(distance(vector<Edge*>::begin(),edge_pos));

      /* inverse adiascence list */
      map<int,set<int> >::iterator to_from_map_pos=to_from_map.find(to_id);
      if (to_from_map_pos==to_from_map.end())
	{
	  to_from_map.insert(pair<const int,set<int> >(to_id,set<int>()));
	  to_from_map_pos=to_from_map.find(to_id);
	}
      /* overkill ... need index of vector<Edge*> */
      to_from_map_pos->second.insert(distance(vector<Edge*>::begin(),edge_pos));
      ++edge_pos;
    }
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






