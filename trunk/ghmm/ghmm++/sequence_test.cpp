#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <xmlio/XMLIO_ArrayObject.h>

#ifdef HAVE_NAMESPACE
using namespace std;
#endif

int main_ArrayObject()
{
  XMLIO_ArrayObject<double> my_array;
  my_array.XMLIO_getCharacters("1.3 1.4 1e4 1 f2dsdf");
  my_array.push_back(0);
  my_array.print();

  XMLIO_ArrayObject<string> my2_array;
  my2_array.XMLIO_getCharacters("1.3 1.4 1e4 1 f2dsdf");
  my2_array.print();
  
  XMLIO_ArrayObject<char> my3_array;
  my3_array.XMLIO_getCharacters("1.3 1.4 1e4 1 f2dsdf");
  my3_array.print();

  return 0;
}

int main()
{
  return main_ArrayObject();
}





