#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/sdmodel.h>
#include <ghmm++/GHMM_convertXMLtoC.h>
#include "sdclass_change.h"

/*int cp_class_change(int *seq, int len) {
  int sum = 0;
  int i;
  for(i=0;i<=len;i++){
	  sum += seq[i];
  }
  //printf("sum = %d\n",sum);
  if (sum >= 6) {
    //printf("\n++++++++++++++++++++++++++++++++Switching class .... ");    
    return 1;
  }
  else {
    return 0;
  } 
} 		
		

void setSwitchingFunction( sdmodel *smd ) {
  smd->get_class = cp_class_change;
} */


/*========================================================================*/


/* the actual pointer to C struct model */

void mat_2d_print(double **mat, int rows, int cols) {
  int i,j;

  printf("double ** Rows = %d, Cols = %d\n", rows, cols);
  for(i=0; i < rows; i++) {
    printf("\t");
    for(j=0; j < cols; j++) {
      printf("%g ", mat[i][j]);
    }
    printf("\n");
  }
}

int __seq_d_class(const double *O, int index, double *osum) {
  return 0; 
}

sdmodel *graphmlTosdModel(char *filename)
{
  model_t   *mymodel;
  sdmodel     *cmodel;
  int nsilents=0;

  mymodel = graphmldoc_cwrapper(filename);

  if (mymodel != NULL)
    {
      if (mymodel->model_id == DISCRETE)
	{
	  cmodel = (sdmodel*)(mymodel->model_pt); /* the actual pointer to C struct model */
	  nsilents = sdmodel_initSilentStates( cmodel );
	  fprintf(stderr, "In Main: # Silent States %d\n", nsilents);

	  if ( cmodel->cos > 1 ) {
	    fprintf(stderr, "Warning: I set the function pointer sdmodel->get_class\n");
	    /* cmodel->get_class = __seq_d_class; */
	    cmodel->get_class = cp_class_change;
	  } else {
	    /* cmodel->get_class = __seq_d_class; */
	    cmodel->get_class = cp_class_change;
	  }

 	  fprintf(stderr, "In read_cxml.c\n");
	  sdmodel_topo_ordering( cmodel );
	  free((void*) mymodel);
	  return cmodel;
	}
      else
	{
	  fprintf(stderr, "graphmlToModel: I need a DISCRETE model\n");
	  free((void*) mymodel);
	  return NULL;
	}
    }
}
