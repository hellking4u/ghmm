/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_alphabet.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>
#include "ghmm_alphabet.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm_alphabet::ghmm_alphabet():
  XMLIO_ElementArrayElement<ghmm_symbol>()
{
  init_members();
}

ghmm_alphabet::ghmm_alphabet(const string& name, XMLIO_Attributes &attrs):
  XMLIO_ElementArrayElement<ghmm_symbol>(name,attrs)
{
  init_members();
  XMLIO_Attributes::const_iterator pos=attrs.find("id");
  if (pos!=attrs.end()) {
    id=pos->second;
  }
  else {
    cerr<<"expected an id attribute for Alphabet"<<endl;
  }
}

const XMLIO_Attributes& ghmm_alphabet::XMLIO_getAttributes() const {
  XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
  attrs["id"]=id;
  return attributes;
}

const char* ghmm_alphabet::toString() const
{
  return "ghmm_alphabet";
}

void ghmm_alphabet::print() const
{
  cout<<toString()<<":"<<endl;
  for (XMLIO_ElementArrayElement<ghmm_symbol>::const_iterator pos=begin();pos!=end();++pos)
    if (*pos!=NULL)
      (*pos)->print();

}

void ghmm_alphabet::init_members()
{
  tag="Alphabet";
  id="";
  set_element_name("Symbol");
}


