/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/InitialStates.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "InitialStates.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

State::State(const string& name, XMLIO_Attributes &attrs)
{
  if (name=="State")
    {
      XMLIO_Attributes::const_iterator id_pos=attrs.find("id");
      if (id_pos==attrs.end())
	cerr<<"no id in State's attributes found!"<<endl;
      else
	id_ref=id_pos->second;
    }
  else
    {
      cerr<<name<<": this is not State, that was expected."<<endl;
    }
}

void State::print() const
{
  cout<<"State "<<id_ref<<endl;
}

const char* State::toString() const
{
  return "State";
}

const string& State::get_id() const
{
  return id_ref;
}

InitialStates::InitialStates(const string& name, XMLIO_Attributes &attrs)
{
  state_pd=NULL;
}

XMLIO_Element* InitialStates::XMLIO_startTag(const string& name, XMLIO_Attributes &attrs)
{
  if (name=="DiscretePD")
    {
      if (state_pd!=NULL)
	{
	  cerr<<"Only one DiscretePD section allowed"<<endl;
	}
      else
	{
	  state_pd=new DiscretePD<State>;
	  /* only State elements are welcome */
	  state_pd->set_element_name("State");
	  return state_pd;
	}
    }
  return NULL;
}


void InitialStates::print() const
{
  cout<<toString()<<endl;
  if (state_pd==NULL)
    {
      cout<<"empty"<<endl;
      return;
    }
  else
    {
      state_pd->print();
    }
}

