/*
  author: Achim Gaedke
  created: 2001-09-07
  file: ghmm/ghmm++/sequences_document.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <iostream>
#include "ghmm++/sequences_document.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

sequences_document::root_element::root_element(const string& name,
					       XMLIO_Attributes& attrs, 
					       list<sequences*>* the_seq_list)
  :XMLIO_Element(name,attrs){
  tag=name;
  attributes=attrs;
  seq_list         = the_seq_list;
  actual_sequences = NULL;
}

sequences_document::root_element::root_element(list<sequences*>* the_seq_list) {
  tag="sequence_data";
  seq_list         = the_seq_list;
  actual_sequences = NULL;
}

    
XMLIO_Element*
sequences_document::root_element::XMLIO_startTag(const string& tag, XMLIO_Attributes& attrs){
  if (tag == "sequences") {
    if (actual_sequences != NULL) {
      cerr<<"root_element::XMLIO_startTag: Internal Error! actual_sequences!=NULL"<<endl;
      SAFE_DELETE(actual_sequences);
    }
    actual_sequences = new sequences(tag,attrs);
    return actual_sequences;
  }
  return NULL;
}

void 
sequences_document::root_element::XMLIO_endTag (const string& tag){
  if (tag!="sequences") return;
  if (seq_list==NULL) {
    cerr<<"expecting existent sequences list!"<<endl;
    return;
  }
  if (actual_sequences==NULL){
    cerr<<"expecting sequences object, but there is nothing!"<<endl;
  }
  seq_list->push_back(actual_sequences);
  actual_sequences=NULL;
}

const int sequences_document::root_element::XMLIO_writeContent (XMLIO_Document& doc) const{
  if (seq_list==NULL) return 0;
  int result=1;
  doc.writeEndl();
  for(list<sequences*>::const_iterator pos=seq_list->begin();
      pos != seq_list->end();
      pos++) {
    if (*pos!=NULL) {
      int tmp_result=doc.writeElement(**pos);
      if (tmp_result<0) return tmp_result;
      doc.writeEndl();
      result+=tmp_result+1;
    }
  }
  return result;
}

sequences_document::sequences_document(){
  my_root=NULL;
}

sequences_document::~sequences_document(){
  SAFE_DELETE(my_root);
  while (!empty())
    {
      SAFE_DELETE(front());
      pop_front();
    }
}

size_t
sequences_document::read_sequences(const string& filename){
  open(filename,"r");
  readDocument();
  close();
  return size();
}

size_t
sequences_document::write_sequences(const string& filename) {
  open(filename,"w");
  if (my_root!=NULL)
    writeElement(*my_root);
  else
    writeElement(root_element(this));
  close();
  return size();
}


XMLIO_Element* 
sequences_document::XMLIO_startTag(const string& tag, XMLIO_Attributes& attributes) {
  if (my_root==NULL)
    my_root=new root_element(tag,attributes,this);
  return my_root;
}

void 
sequences_document::XMLIO_endTag(const string& tag){
  SAFE_DELETE(my_root);
}
