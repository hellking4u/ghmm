/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/xml_parser_test.c
  created      : DATE: May 2001
  $Id$

__copyright__

*******************************************************************************/

#include "config.h"
#if defined(__EXPERIMENTAL__) && __EXPERIMENTAL__ == 1 &&  defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT)

#include <ghmm/xmlio.h>

int main()
{
  const char* filename="simple_graph.xml";
  document_handler* my_dh;
  int result;
  my_dh=create_document_handler(filename);
  if (my_dh==NULL)
    {
      fprintf(stderr,"something wrong\n");
      return 1;
    }
  result=parse_document(my_dh);
  delete_document_handler(my_dh);
  return 0;
}

#else

#include <stdio.h>

int main()
{
  fprintf(stderr,"xml not implemented.\n");
  return 0;
}

#endif /* defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT) && __experimental__ */
