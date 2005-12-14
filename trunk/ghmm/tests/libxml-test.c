#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <ghmm/ghmm.h>
#include <ghmm/xmlreader.h>


/*===========================================================================*/
int main(int argc, char **argv) {

  char *docname;
  fileData_s * f;
  ghmm_set_loglevel(5+1);

  if(argc <= 1) {
    printf("Usage: %s docname.xml", argv[0]);
    return(0);
  }

  docname = argv[1];
  f = parseHMMDocument(docname);

  /* simple test */
  if (f) {
    switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			    + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
    case GHMM_kContinuousHMM:
      ghmm_c_print(stdout, f->model.c);
      break;
    case GHMM_kDiscreteHMM:
      ghmm_d_print(stdout, f->model.d);
    default:
      break;
    }
  }

  return(0);
}
