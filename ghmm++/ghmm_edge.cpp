/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_edge.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xmlio/XMLIO_ArrayElement.h>
#include <xmlio/XMLIO_StringElement.h>
#include <gxl/gxl_attr.h>
#include "ghmm_edge.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

ghmm_edge::ghmm_edge():
  gxl_edge()
{}

ghmm_edge::ghmm_edge(const string& name, XMLIO_Attributes &attrs):
  gxl_edge(name,attrs)
{}

XMLIO_Element* ghmm_edge::XMLIO_startTag(const string& name, XMLIO_Attributes &attrs)
{
  /* handle weight class information */
  return gxl_edge::XMLIO_startTag(name,attrs);
}

void ghmm_edge::XMLIO_endTag(const string& name)
{
  gxl_edge::XMLIO_endTag(name);
  if (name=="attr")
    {
      /* find out, if there are transitions*/
      gxl_attr* last_attr=my_attrs.back();
      if (last_attr->name=="transitions" && !last_attr->empty())
	{
	  /* only first attribute is tested */
	  XMLIO_ElementArrayElement<XMLIO_StringElement>* int_array;
	  int_array=dynamic_cast<XMLIO_ElementArrayElement<XMLIO_StringElement>*>(last_attr->front());
	  /* is it possible to get transition probs ? */
	  if (int_array!=NULL)
	    {
	      for(XMLIO_ElementArrayElement<XMLIO_StringElement>::const_iterator pos=int_array->begin();
		  pos!=int_array->end();
		  ++pos)
		{
		  double* new_prob=NULL;
		  (void)XMLIO_evaluate_token(**pos,0,(*pos)->size(),new_prob);
		  if (new_prob!=NULL)
		    {
		      push_back(*new_prob);
		      delete new_prob;
		    }
		}/* for pos in int_array */
	    }
	}
    }
}

const char* ghmm_edge::toString() const
{
  return "ghmm_edge";
}

void ghmm_edge::print() const
{
  cout<<toString()<<": "<<from<<"->"<<to<<endl;
  if (!empty())
    {
      int count=0;
      for (ghmm_edge::const_iterator pos=begin();pos!=end();++pos)
	cout<<"Transition probability ("<<count++<<"):"<<*pos<<endl;
    }
}


