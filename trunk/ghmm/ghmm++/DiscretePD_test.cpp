#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <xmlio/XMLIO_ObjectReader.h>
#include "sequences.h"
#include "DiscretePD.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main()
{
  XMLIO_ElementArrayReader<InitialStates> ghmm_doc;
  ghmm_doc.set_doc_name("ghmm");
  ghmm_doc.set_element_name("InitialStates");
  (void)ghmm_doc.read_file("ghmm.xml");

  XMLIO_ElementArrayReader<InitialStates>::const_iterator pos=ghmm_doc.begin();
  while(pos!=ghmm_doc.end())
    {
      cout<<"element"<<endl;
      (*pos)->print();
      ++pos;
    }

  return 0;
}
