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
#include <xmlio/XMLIO_ObjectReader.h>
#include <ghmm++/ghmm.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main()
{
  XMLIO_ElementReader<hmm> my_model_reader;
  my_model_reader.set_doc_name("hmm");
  my_model_reader.read_file("ghmm.xml");
  if (my_model_reader.get_element()!=NULL)
    {
      cout<<my_model_reader.toString()<<"found something: "<<endl;
      my_model_reader.get_element()->print();
    }
  model* my_model=my_model_reader.get_element()->create_model();
  model_print(stdout,my_model);
  return 0;
}


