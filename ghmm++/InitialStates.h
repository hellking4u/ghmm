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

#include <xmlio/XMLIO_Object.h>
#include <ghmm++/DiscretePD.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
     used in InitialStates/DiscretePD for initial State distributions
   */
  class State: public XMLIO_Object
    {
    public:
      State(const string& tag, XMLIO_Attributes &attributes);
      void print() const;
      const char* toString() const;
    private:
      string id_ref;
    };
  
  class InitialStates: public XMLIO_Object
    {
    public:
      InitialStates(const string& tag, XMLIO_Attributes &attributes);
      virtual XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      virtual const char* toString() const {return "InitialStates";}
      virtual void print() const;
    private:
      DiscretePD<State>* state_pd;
    };


#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_INITIAL_STATES_H */

