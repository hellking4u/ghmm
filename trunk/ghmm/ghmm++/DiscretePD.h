/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/DiscretePD.h
  $Id$
 */

#ifndef GHMMPP_DISCRETEPD_H
#define GHMMPP_DISCRETEPD_H

#include <iostream>
#include <vector>
#include <xmlio/XMLIO_Object.h>
#include <xmlio/XMLIO_ObjectReader.h>
#include <xmlio/XMLIO_SkipObject.h>
#include <xmlio/XMLIO_IgnoreObject.h>
#include <xmlio/XMLIO_ArrayObject.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  template <class E>
    class DiscretePD: public XMLIO_ContentElementArrayObject<double,E>
    {
    public:
      /**
	 if no weight is given default_weight is 1.
       */
      DiscretePD()
	{
	  default_weight=1;
	  actual_weight=default_weight;
	}

      /**
	 looks for attribute default_weight
       */
      DiscretePD(const string& tag, XMLIO_Attributes &attributes)
	{
	  XMLIO_Attributes::iterator default_weight_key=attributes.find("default_weight");
	  if (default_weight_key!=attributes.end())
	    {
	      double* tmp_default;
	      (void)XMLIO_evaluate_token(pos->second,0,pos->second->size(),&tmp_default);
	      if (tmp_default==NULL)
		{
		  cerr<<toString()<<": default weight is not a floating point value"<<endl;
		  default_weight=1;
		}
	      else
		{
		  default_weight=*tmp_default;
		  delete tmp_default;
		}
	    }
	  else
	    default_weight=1;
	  actual_weight=default_weight;
	}
    private:
      double default_weight;
      double actual_weight;
    };


#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_DISCRETEPD_H */



