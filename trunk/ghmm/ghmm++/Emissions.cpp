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

Emission::Emission(const string& tag, XMLIO_Attributes &attributes):
  XMLIO_ContentElementArrayObject<double,XMLIO_AttributedObject>(tag,attributes)
{
  fix=0;
  XMLIO_Attributes::const_iterator pos;
  /* madatory argument */  
  pos=attributes.find("State");
  if (pos!=attributes.end())
    {
      cerr<<toString()<<"State id is missing!"<<endl;
    }
  else
    {
      state=pos->second;
    }
  pos=attributes.find("fix");
  if (pos!=attributes.end())
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

Emissions::Emissions(const string& tag, XMLIO_Attributes &attributes)
{
  set_element_name("Emission");
}
