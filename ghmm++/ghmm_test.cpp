/*
  author: Achim Gaedke
  created: 9. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_test.cpp
  $Id$
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <xmlio/XMLIO_Document.h>
#include <ghmm++/ghmm.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

class ghmm_Document: public XMLIO_Document {
  
public:
  ghmm_Document(){
    root=NULL;
  }

  virtual XMLIO_Element*
  XMLIO_startTag (const string& name, XMLIO_Attributes &attrs) {
    if (root == NULL) {
      root=new ghmm(name,attrs);
    }
    else {
      cerr<<"found a second root element!"<<endl;
    }
    return root;
  }

  virtual ghmm* get_ghmm() {
    return root;
  }

  ~ghmm_Document() {
    SAFE_DELETE(root);
  }

  ghmm* root;
};


int sequence_generation_test() {
  ghmm_Document my_model;
  my_model.open("ghmm.xml","r");
  my_model.readDocument();
  my_model.close();
  ghmm* my_ghmm=my_model.get_ghmm();
  if (my_ghmm==NULL) return 1;

  sequences* generated=my_ghmm->generate_sequences(10,10);
  if (generated!=NULL)
    generated->print();
  SAFE_DELETE(generated);
  return 0;
}

int read_test_file() {
  ghmm_Document my_model;
  if (my_model.open("ghmm.xml","r") != 0) {
    cerr<<"error while opening ghmm.xml"<<endl;
  }
  if (my_model.readDocument() != 0) {
    cerr<<"error while reading ghmm.xml"<<endl;
  }
  my_model.close();
  ghmm* my_ghmm=my_model.get_ghmm();

  if (my_ghmm != NULL) {
    XMLIO_Document output;
    output.open("/dev/stdout","w");
    output.writeElement(*my_ghmm);
    output.close();
  }
  else {
    cout<<"no ghmm document found"<<endl;
  }
  return 0;
}

int main()
{
  read_test_file();
  /*
  sequence_generation_test();
  */
  
  return 0;
}


