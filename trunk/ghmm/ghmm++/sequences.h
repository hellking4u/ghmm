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
#include <xmlio/XMLIO_Element.h>
#include <xmlio/XMLIO_StringElement.h>
#include <xmlio/XMLIO_ObjectReader.h>
#include <xmlio/XMLIO_ArrayElement.h>
#include <ghmm/sequence.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

  /**
   */
  class double_sequence: public XMLIO_ArrayElement<double>
    {
    public:
      /**
       */
      double_sequence(double* seq_data, size_t length);
      /**
       */
      double_sequence(int* seq_data, size_t length);
      /**
       */
      double_sequence(){label=NULL; id="";}
      /**
       */
      double_sequence(const string& name, XMLIO_Attributes &attributes);
      /**
       */
      ~double_sequence();
      /**
       */
      XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      double* create_double_array() const;
      /**
       */
      int* create_int_array() const;
      /**
       */
      int get_label_as_int() const;
      /**
       */
      double get_id_as_double() const;
      /**
       */
      virtual void print() const;
    private:
      /**
       */
      XMLIO_StringElement* label;
      /**
       */
      string id;
    };

  /**
   */
  class int_sequence: public XMLIO_ArrayElement<int>
    {
    public:
      /**
       */
      int_sequence(){label=NULL; id="";}
      /**
       */
      int_sequence(int* seq_data, size_t length);
      /**
       */
      int_sequence(double* seq_data, size_t length);
      /**
       */
      int_sequence(const string& name, XMLIO_Attributes &attributes);
      /**
       */
      ~int_sequence();
      /**
       */
      XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attributes);
      /**
       */
      double* create_double_array() const;
      /**
       */
      int* create_int_array() const;
      /**
       */
      int get_label_as_int() const;
      /**
       */
      double get_id_as_double() const;
      /**
       */
      virtual void print() const;
    private:
      /**
       */
      XMLIO_StringElement* label;
      /**
       */
      string id;
    };

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

  /**
   */
  class sequenceReader: public XMLIO_ObjectReader, public list<sequences*>
    {
    public:
      /**
       */
      sequenceReader();
      /**
       */
      ~sequenceReader();
      /**
       */
      size_t read_sequences(const string& filename);
      /**
       */
      XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes& attributes);
      /**
       */
      void XMLIO_endTag(const string& tag);
      
    private:
      /**
       */
      sequences* next_sequence_array;
    };

#ifdef HAVE_NAMESPACES
}
#endif


#endif /* GHMMPP_SEQUENCES_H */
