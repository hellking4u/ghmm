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
 return ghmm_xml_parse("simple_graph.xml");
}

#else

#include <stdio.h>

int main()
{
  fprintf(stderr,"xml not implemented.\n");
  return 0;
}

#endif /* defined(HAVE_EXPAT_H) && defined(HAVE_LIBEXPAT) && __experimental__ */
