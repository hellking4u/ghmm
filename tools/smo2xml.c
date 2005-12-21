#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <ghmm/ghmm.h>
#include <ghmm/xmlreader.h>
#include <ghmm/smodel.h>
#include <ghmm/obsolete.h>


/*===========================================================================*/
int main(int argc, char **argv) {

  char *docname, *writename;
  fileData_s * f;
  int i;
  int mo_number = 0;
  ghmm_cmodel ** smo;
  ghmm_set_loglevel(5+1);

  if(argc <= 1) {
    printf("Usage: %s docname.smo docname.xml", argv[0]);
    return(0);
  }

  docname = argv[1];

  smo = ghmm_c_read(docname, &mo_number);

  printf("Models %d", mo_number);
  
  /* simple test */
  if (smo) {
    writename = argv[2];
    ghmm_c_xml_write(writename,smo,mo_number);
  }


  return(0);
}
