/*
  author: Achim Gaedke
  created: 2001-09-07
  file: ghmm/ghmm++/sequences_document.h
  $Id$
 */

#ifndef GHMMPP_SEQUENCES_DOCUMENT_H
#define GHMMPP_SEQUENCES_DOCUMENT_H

#include <string>
#include <list>
#include <iostream>
#include <xmlio/XMLIO_Document.h>
#include <ghmm++/sequences.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

/**
 */
class sequences_document: public XMLIO_Document, public list<sequences*> {

  /**
     should take the root element
   */
  class root_element: public XMLIO_Element{

  public:
    /** where the sequences-element should go */
    list<sequences*>* seq_list;
    /** the processed sequences-element */
    sequences* actual_sequences;

    /** for other instantiation than reading */
    root_element(list<sequences*>* the_seq_list);

    /** for read instantiation */
    root_element(const string& name, XMLIO_Attributes& attrs, list<sequences*>* the_seq_list);

    /** look out for sequences-elements */
    virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes& attrs);

    /** */
    virtual void XMLIO_endTag (const string& tag);

    /** root element writer */
    virtual const int XMLIO_writeContent (XMLIO_Document&) const;

  };

public:
  /**
     initialise reader
   */
  sequences_document();

  /**
     delete content of array
   */
  ~sequences_document();

  /**
     this is the real read procedure
   */
  size_t read_sequences(const string& filename);


  /**
     this is the real read procedure
   */
  size_t write_sequences(const string& filename);

  /**
     accept all document roots!
   */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes& attributes);

  /**
   */
  virtual void XMLIO_endTag(const string& tag);
  
private:
  /**
   */
  root_element* my_root;
};

#ifdef HAVE_NAMESPACES
}
#endif


#endif /* GHMMPP_SEQUENCES_DOCUMENT_H */
