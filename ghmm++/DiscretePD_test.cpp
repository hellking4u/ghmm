#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "sequences.h"
#include "DiscretePD.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main()
{
  XMLIO_ElementArrayReader<InitialStates> pd;
  pd.set_element_name("InitialStates");
  (void)pd.read_file("ghmm.xml");

  XMLIO_ElementArrayReader<InitialStates>::const_iterator pos=pd.begin();
  while(pos!=pd.end())
    {
      cout<<"element"<<endl;
      (*pos)->print();
      ++pos;
    }

  return 0;
}
