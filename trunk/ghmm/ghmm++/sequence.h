/*
  author: Achim Gaedke
  created: 2001-09-07
  file: ghmm/ghmm++/sequence.h
  $Id$
 */

#ifndef GHMMPP_SEQUENCE_H
#define GHMMPP_SEQUENCE_H

#include <string>
#include <vector>
#include <list>
#include <iostream>
#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Element.h>
#include <xmlio/XMLIO_StringElement.h>
#include <xmlio/XMLIO_ArrayElement.h>
#include <ghmm/sequence.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif


  /**
     base class for all sequences
   */
  template<typename T>
  class sequence: public XMLIO_ArrayElement<T> {
  public:
    /**
     */
    string label;

    /**
     */
    string id;

    /**
     */
    sequence<T>() {
      label = "";
      id    = "";
    }

    /**
     */
    sequence<T>(const string& name, XMLIO_Attributes &attrs) {
      tag=name;
      XMLIO_Attributes::iterator id_key=attrs.find("id");
      XMLIO_Attributes::iterator label_key=attrs.find("label");
      if (id_key != attrs.end()) {	
	id=id_key->second;
      }
      if (label_key != attrs.end()) {
	label=label_key->second;
      }
    }

    /**
     */
    virtual const XMLIO_Attributes& XMLIO_getAttributes () const {
      XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
      attrs.clear();
      if (!id.empty()) 
	attrs["id"]=id;
      if (!label.empty()) 
	attrs["label"]=label;
      return attributes;
    }

    /**
     */
    virtual int get_label_as_int() const  {
      if (label.empty()) {
	return -1;
      }
      else {
	return strtol(label.c_str(),NULL,0);
      }      
    }

    /**
     */
    virtual double get_id_as_double() const {
      if (id.empty())
	return -1.0;
      else
	return strtod(id.c_str(),NULL);
    }
  };

  /**
   */
  class double_sequence: public sequence<double> {
  public:
    /**
     */
    double_sequence(const string& name, XMLIO_Attributes &attrs);
    /**
     */
    double_sequence(double* seq_data, size_t length);
    /**
     */
    double_sequence(int* seq_data, size_t length);
    /**
     */
    double_sequence(sequence_d_t* seq, int sequence_pos);
    /**
     */
    virtual double* create_double_array() const;
    /**
     */
    virtual int* create_int_array() const;
  };

  /**
   */
  class int_sequence: public sequence<int> {
  public:
    /**
     */
    int_sequence(const string& name, XMLIO_Attributes &attrs);
    /**
     */
    int_sequence(int* seq_data, size_t length);
    /**
     */
    int_sequence(double* seq_data, size_t length);
    /**
     */
    int_sequence(sequence_t* seq, int sequence_pos); 
    /**
     */
    virtual double* create_double_array() const;
    /**
     */
    virtual int* create_int_array() const;
  };


#ifdef HAVE_NAMESPACES
}
#endif


#endif /* GHMMPP_SEQUENCES_H */
