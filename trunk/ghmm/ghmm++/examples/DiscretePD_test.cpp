#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <xmlio/XMLIO_Document.h>
#include <ghmm++/InitialStates.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main()
{
  XMLIO_ElementReader<InitialStates> ghmm_doc;
  ghmm_doc.set_doc_name("ghmm");
  (void)ghmm_doc.read_file("ghmm.xml");

  return 0;
}
