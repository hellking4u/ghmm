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

  /**
     a map consisting of doubles followed by elements E
   */
  template <class E>
    class DiscretePD: public map<E*,double>, public XMLIO_ContentElementArrayObject<double,E>
    {
    public:
      /**
	 if no weight is given default_weight is 1.
       */
      DiscretePD()
	{
	  default_weight=1;
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

      /**
	 makes map from ContentElementArray
       */
      virtual void XMLIO_finishedReading()
	{
	  /* parse pair vector and collect information into a map */
	  double actual_weight=default_weight;
	  for (XMLIO_ContentElementArrayObject<double,E>::iterator pos=XMLIO_ContentElementArrayObject<double,E>::begin();
	       pos!=XMLIO_ContentElementArrayObject<double,E>::end();
	       ++pos)
	    {
	      if (pos->content!=NULL)
		{
		  actual_weight=*(pos->content);
		}
	      if (pos->element!=NULL)
		{
		  map<E*,double>::insert(pair<E*,double>(pos->element,actual_weight));
		  actual_weight=default_weight;
		}
	    }
	}

      /**
       */
      virtual void print() const
	{
	  cout<<toString()<<endl;
	  for (map<E*,double>::const_iterator pos=map<E*,double>::begin();
	       pos!=map<E*,double>::end();
	       ++pos)
	    {
	      pos->first->print();
	      cout<<":"<<pos->second<<endl;
	    }
	}

    private:
      /**
       */
      double default_weight;
    };


#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_DISCRETEPD_H */
