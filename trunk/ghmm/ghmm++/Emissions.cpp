/*
  author: Achim Gaedke
  created: 13. Juli 2001
  file: xmlio/examples/ghmm++/Emissions.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "Emissions.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

Emission::Emission()
{
  fix=0;
}

const char* Emission::toString() const
{
  return "Emission";
}


Emission::Emission(const string& name, XMLIO_Attributes &attrs):
  XMLIO_ContentElementArrayElement<double,XMLIO_AttributedElement>(name,attrs)
{
  fix=0;
  XMLIO_Attributes::const_iterator pos;
  /* madatory argument */  
  pos=attrs.find("State");
  if (pos==attrs.end())
    {
      cerr<<toString()<<"State id is missing!"<<endl;
    }
  else
    {
      state=pos->second;
    }
  pos=attrs.find("fix");
  if (pos!=attrs.end())
    {
      if (pos->second.empty() || pos->second=="yes")
	{
	  fix=1;
	}
      else if (pos->second=="no")
	{
	  fix=0;
	}
      else
	{
	  cerr<<toString()<<": in attribute fix yes, no or '' expected"<<endl;
	}
    }
}

Emissions::Emissions()
{
  set_element_name("Emission");
}

Emissions::Emissions(const string& name, XMLIO_Attributes &attrs)
{
  set_element_name("Emission");
}

const char* Emissions::toString() const
{
  return "Emissions";
}
