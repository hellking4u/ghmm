#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <ghmm/ghmm.h>
#include <ghmm/xmlreader.h>
#include <ghmm/xmlwriter.h>


/*===========================================================================*/
int main(int argc, char **argv) {

  char *docname, *writename;
  ghmm_fileData_s * f;
  int i;
  ghmm_set_loglevel(5+1);

  if(argc <= 1) {
    printf("Usage: %s docname.xml", argv[0]);
    return(0);
  }

  docname = argv[1];
  f = ghmm_parseHMMDocument(docname);
  /* simple test */
  if (f) {
    for (i=0;i<f->noModels; i++){
      switch (f->modelType & (GHMM_kDiscreteHMM + GHMM_kTransitionClasses
			      + GHMM_kPairHMM + GHMM_kContinuousHMM)) {
      case GHMM_kContinuousHMM:
        ghmm_c_print(stdout, f->model.c[i]);
        break;
      case GHMM_kDiscreteHMM:
        ghmm_d_print(stdout, f->model.d[i]);
      default:
        break;
      }
    }
  }

  if (argc > 2) {
    writename = argv[2];
    ighmm_writeHMMDocument(f, writename);
  }

  return(0);
}
