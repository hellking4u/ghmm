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

XMLIO_ElementReader<ghmm> my_model_reader;


int create_test()
{
  model* my_model=my_model_reader.get_element()->create_model();
  model_print(stdout,my_model);
  model_free(&my_model);
  return 0;
}

int sequence_generation_test()
{
  sequences* generated=my_model_reader.get_element()->generate_sequences(10,10);
  if (generated!=NULL)
    generated->print();
  SAFE_DELETE(generated);
  return 0;
}

int read_test_file()
{
  my_model_reader.set_doc_name("ghmm");
  my_model_reader.read_file("ghmm.xml");
  if (my_model_reader.get_element()!=NULL)
    {
      cout<<my_model_reader.toString()<<" found something: "<<endl;
      my_model_reader.get_element()->print();
    }
  return 0;
}

int main()
{
  read_test_file();
  sequence_generation_test();
  
  return 0;
}


