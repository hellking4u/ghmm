/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/InitialStates.h
  $Id$
 */

#ifndef GHMMPP_INITIAL_STATES_H
#define GHMMPP_INITIAL_STATES_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <map>
#include <xmlio/XMLIO_Element.h>
#include <ghmm++/DiscretePD.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     used in InitialStates/DiscretePD for initial State distributions
   */
  class State: public XMLIO_Element
    {
    public:
      /**
       */
      State(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      void print() const;
      /**
       */
      virtual const char* toString() const;
      /**
       */
      virtual const string& get_id() const;
    private:
      /**
       */
      string id_ref;
    };

  /**
     container with all initial states probabilities
   */  
  class InitialStates: public XMLIO_Element
    {
    public:
      /**
       */
      InitialStates(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      virtual const char* toString() const {return "InitialStates";}
      /**
       */
      virtual void print() const;
      /**
       */
      virtual map<State*,double>* get_map()
	{
	  return state_pd;
	}
      
    private:
      /**
       */
      DiscretePD<State>* state_pd;
    };


#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_INITIAL_STATES_H */

