#include "XMLIO_ArrayObject.h"

int main()
{
  XMLIO_Array<double> my_array;
  my_array.XMLIO_getCharacters("1 1 1 1");
  my_array.push_back(1);
  my_array.print();
  return 0;
}
