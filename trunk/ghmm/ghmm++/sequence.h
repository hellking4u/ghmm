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

  class sequence_label: public XMLIO_StringElement{
  public:
    sequence_label(const string& name, XMLIO_Attributes &attrs) {
      tag=name;
    }
  };

  /**
     base class for all sequences
   */
  template<typename T>
  class sequence: public XMLIO_ArrayElement<T> {
  public:
    /**
     */
    sequence_label* label;

    /**
     */
    string id;

    /**
     */
    sequence<T>() {
      label = NULL;
      id    = "";
    }

    /**
     */
    sequence<T>(const string& name, XMLIO_Attributes &attrs) {
      tag=name;
      XMLIO_Attributes::iterator id_key=attrs.find("id");
      if (id_key != attrs.end())
	{
	  id=id_key->second;
	}
      label=NULL;
    }

    /**
     */
    virtual ~sequence<T>() {
      SAFE_DELETE(label);
    }
 
    /**
       gets the label of this sequence
     */
    virtual XMLIO_Element* XMLIO_startTag(const string& name, XMLIO_Attributes &attrs) {
      if (name=="label") {
	if (label==NULL) {
	  label=new sequence_label(name, attrs);
	  return label;
	}
	else {
	  cout<<"Only one label allowed in sequence"<<endl;
	}
      }
      else {
	cout<<tag<<" not supported in sequence"<<endl;
      }
      return NULL;
    }

    /**
     */
    virtual const XMLIO_Attributes& XMLIO_getAttributes () const {
      XMLIO_Attributes& attrs=(XMLIO_Attributes&)attributes;
      attrs.clear();
      if (!id.empty())
	attrs["id"]=id;
      return attributes;
    }


    /**
     */
    virtual const int XMLIO_writeContent (XMLIO_Document& doc) const {
      int result=0;
      int tmp_result;
      if (label!=NULL) {
	tmp_result=doc.writeElement(*label);
	if (tmp_result<0) return tmp_result;
      }
      tmp_result=XMLIO_ArrayElement<T>::XMLIO_writeContent(doc);
      if (tmp_result<0) return tmp_result;
      return result+tmp_result;
    }

    virtual void print() const {
      XMLIO_ArrayElement<T>::print();
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
    virtual double* create_double_array() const;
    /**
     */
    virtual int* create_int_array() const;
    /**
     */
    virtual int get_label_as_int() const;
    /**
     */
    virtual double get_id_as_double() const;
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
    virtual double* create_double_array() const;
    /**
     */
    virtual int* create_int_array() const;
    /**
     */
    virtual int get_label_as_int() const;
    /**
     */
    virtual double get_id_as_double() const;
  };


#ifdef HAVE_NAMESPACES
}
#endif


#endif /* GHMMPP_SEQUENCES_H */
