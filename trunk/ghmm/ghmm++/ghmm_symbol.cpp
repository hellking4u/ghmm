/*
  author: Achim Gaedke
  created: 18. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_symbol.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include <xmlio/XMLIO_ArrayElement.h>
#include "ghmm_symbol.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm_symbol::ghmm_symbol():
  XMLIO_StringElement()
{
  init_members();
}

ghmm_symbol::ghmm_symbol(int my_code, const string& my_symbol):
  XMLIO_StringElement(my_symbol)
{
  init_members();
  valid=1;
  code=my_code;
}

ghmm_symbol::ghmm_symbol(const string& tag, XMLIO_Attributes &attributes):
  XMLIO_StringElement(tag,attributes)
{
  init_members();
  XMLIO_Attributes::const_iterator pos=attributes.find("code");
  if (pos!=attributes.end())
    {
      int* new_code=NULL;
      (void)XMLIO_evaluate_token(pos->second,0,pos->second.size(),new_code);
      if (new_code!=NULL)
	{
	  valid=1;
	  code=*new_code;
	  delete new_code;
	}
      else
	cerr<<toString()<<": code attribute value is no integer"<<endl;
    }
  else
    {
      cerr<<toString()<<": code attribute required!"<<endl;
    }
}

const char* ghmm_symbol::toString() const
{
  return "ghmm_symbol";
}

void ghmm_symbol::print() const
{
  if (valid)
    cout<<toString()<<": code="<<code<<" value="<<*this<<endl;
  else
    cout<<toString()<<": element not valid"<<endl;
}

void ghmm_symbol::init_members()
{
  valid=0;
  code=-1;
}


