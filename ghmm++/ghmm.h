/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm.h
  $Id$
 */

#ifndef GHMMPP_GHMM_H
#define GHMMPP_GHMM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xmlio/XMLIO_Object.h>
#include <ghmm++/sequences.h>
#include <ghmm++/graph.h>
#include <ghmm++/InitialStates.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /** 
      @memo class containing hmm
   */
  class hmm: public XMLIO_Object
    {
    public:
      /** constructor from XMLIO interface */
      hmm(const string& tag, XMLIO_Attributes &attributes);
      /**
	 expected elements are: Graph, Initial, Emissions and Sequences
       */
      virtual XMLIO_Object* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /** */
      virtual void XMLIO_endTag(const string& tag);
      /** */
      virtual void XMLIO_getCharacters(const string& characters);
      /** */
      virtual void XMLIO_finishedReading();
      /** Returns name of class. */
      virtual const char* toString() const;
      /** dumps all content to cout */
      virtual void print() const;

    protected:
      /**
	 safe the type of HMM.
	 d= discrete
	 c= continuous
	 s= switched (different classes of transition)
      */
      char type;
      /**
	 saves the initial states, as they come from xmlio
       */
      InitialStates* Initial;
      /**
	 saves a discrete hmm
       */
      model* hmm_model;
      /**
	 saves a continuous/switched hmm
       */
      smodel* shmm_model;
    };

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_H */
