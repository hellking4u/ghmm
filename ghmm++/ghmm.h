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

#include <ghmm++/ghmm_graph.h>
#include <ghmm++/ghmm_alphabet.h>
#include <ghmm++/sequences.h>
#include <ghmm++/InitialStates.h>
#include <ghmm++/Emissions.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /** 
      @memo class containing hmm
   */
  class ghmm: public XMLIO_Element
    {
    public:
      /** constructor from XMLIO interface */
      ghmm(const string& tag, XMLIO_Attributes &attributes);
      /**
	 expected elements are: Graph, Initial, Emissions and Sequences
       */
      virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
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
      /**
	 return id of this model
       */
      const string& get_id() const;
      /**
	 creates a model structure from collected data.
	 if data are not suitable, NULL is returned.
      */
      model* create_model() const;
      /**
	 creates a smodel structure from collected data.
	 if data are not suitable, NULL is returned.
	 not yet implemented
      */
      smodel* create_smodel() const;

      /**
	 generate sequences from model
       */
      sequences* generate_sequences(int number, int length);

      /**
       */
      ~ghmm();

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
	 saves the emissions informations
       */
      Emissions* ghmm_Emissions;

      /**
	 apriori probability of this model
       */
      double prior;
      /**
	 id is mandatory for each Element in gxl (!?)
       */
      string id;
      /**
       */
      ghmm_graph* my_graph;
      /**
       */
      ghmm_alphabet* my_alphabet;
    };

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_H */


