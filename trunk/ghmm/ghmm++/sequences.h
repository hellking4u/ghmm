/*
  author: Achim Gaedke
  created: 26. June 2001
  file: xmlio/examples/ghmm/sequences.h
  $Id$
 */

#ifndef GHMMPP_SEQUENCES_H
#define GHMMPP_SEQUENCES_H

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Document.h>
#include <xmlio/XMLIO_Element.h>
#include <xmlio/XMLIO_StringElement.h>
#include <xmlio/XMLIO_ArrayElement.h>
#include <ghmm++/sequence.h>
#include <ghmm/sequence.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
   */
  class sequences_DiscretePD: public XMLIO_Element
    {
    public:
      /**
       */
      sequences_DiscretePD(int** data, double* weight_data, size_t length, size_t number);
      /**
       */
      sequences_DiscretePD(double** data, double* weight_data, size_t length, size_t number);
      /**
       */
      sequences_DiscretePD(sequence_t* seq);
      /**
       */
      sequences_DiscretePD(sequence_d_t* seq);
      /**
       */
      sequences_DiscretePD(XMLIO_Attributes& attributes, const string& type);
      /**
       */
      ~sequences_DiscretePD();
      /**
       */
      XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      void XMLIO_endTag(const string& tag);
      /**
       */
      void XMLIO_finishedReading();
      /**
       */
      void XMLIO_getCharacters(const string& characters);
      /**
       */
      virtual const XMLIO_Attributes& XMLIO_getAttributes () const;
      /**
       */
      virtual const int XMLIO_writeContent (XMLIO_Document& doc) const;
      /**
       */
      const string& get_type() const {return type;}
      /**
       */
      sequence_t* create_sequence_t() const;
      /**
       */
      sequence_d_t* create_sequence_d_t() const;
      /**
       */
      virtual void print() const;
      
    private:
      /**
       */
      string type;
      /**
       */
      double actual_weight;
      /**
       */
      double default_weight;
      /**
       */
      vector<int_sequence*> int_sequence_vector;
      /**
       */
      vector<double_sequence*> double_sequence_vector;
      /**
       */
      vector<double> weight_vector;
    };

  /**
   */
  class sequences: public XMLIO_Element
    {
    public:
      /**
       */
      sequences();
      /**
       */
      sequences(int** data, double* weight, size_t length, size_t number);
      /**
       */
      sequences(double** data, double* weight, size_t length, size_t number);
      /**
       */
      sequences(sequence_t* seq);
      /**
       */
      sequences(sequence_d_t* seq);
      /**
       */
      sequences(const string& tag, XMLIO_Attributes& attributes);
      /**
       */
      ~sequences();
      /**
       */
      XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      void XMLIO_endTag(const string& tag);
      /**
       */
      void XMLIO_finishedReading();
      /**
       */
      const char* toString() const;
      /**
       */
      const string& get_type() const {return type;}
      /**
       */
      sequence_t* create_sequence_t() const;
      /**
       */
      sequence_d_t* create_sequence_d_t() const;
      /**
       */
      virtual void print() const;
      /**
       */
      virtual const XMLIO_Attributes& XMLIO_getAttributes () const;
      /**
       */
      virtual bool XMLIO_isEmpty () const;
      /**
       */
      virtual const int XMLIO_writeContent (XMLIO_Document& doc) const;

    private:
      /**
       */
      sequences_DiscretePD* sequence_array;
      /**
       */
      string type;
      /**
       */
      string coding;
    };

#ifdef HAVE_NAMESPACES
}
#endif


#endif /* GHMMPP_SEQUENCES_H */
