/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#undef NDEBUG
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/sdmodel.h>
#include <ghmm++/GHMM_convertXMLtoC.h>

/* What this code does:

 * Testing if we can call our C++ wrapper from C
 * 1. Read an XML file in GraphML format
 * 2. Check the model (HACK!!! no checking because of silent states)
 * 3. Generate 100 sequences 
 * 4. Change the B matrix slightly and train with Baum-Welch 
 */

int seq_d_class(int N, int matchcount) {

  if (matchcount >= 2) {
    return 1;
    fprintf(stderr, "----- SWITCHING CLASS -----\n");
   } else
    return 0; 
}


int main(int argc, char* argv[]) 
{
  /* Be careful!! : these are pointers to C-struct 
   */
  model_t   *mymodel;
  sequence_t*my_output;
  sdmodel   *cmodel;
  sdstate   *state_0_pt;

  int    i;
  int    num_sequences=5;
  int    nsilents=0;
  double log_p;
  int seq[5]= {0, 1, 2, 3, 4};
  int *vipath;
  int seq_len = 5;

  if (argc < 2)
    {
      fprintf(stderr, "Please input a filename. Usage: read_cxml file.gml\n");
      exit(-1);
    }

  /* Important! initialise rng  */
  ghmm_rng_init();


  mymodel = graphmldoc_cwrapper(argv[1]);


  if (mymodel != NULL)
    {
      if (mymodel->model_id == DISCRETE)
	{
	  printf("transition matrix A:\n");
	  sdmodel_Ak_print(stdout, (sdmodel*) mymodel->model_pt, 0, " ", ",", "\n");
	  printf("Observation matrix B:\n");
	  sdmodel_B_print(stdout, (sdmodel*) mymodel->model_pt, " ", ",", "\n");

	  cmodel = (sdmodel*)(mymodel->model_pt); /* the actual pointer to C struct model */
	  state_0_pt = &(cmodel->s[0]);          /* 1st state */

	  cmodel = (sdmodel*) (mymodel->model_pt);
          cmodel->get_class = seq_d_class;
	  nsilents = sdmodel_initSilentStates(cmodel);
	  fprintf(stderr, "In Main: # Silent States %d\n", nsilents);
	  sdmodel_topo_ordering(cmodel);
	}

      assert( ((sdmodel*)(mymodel->model_pt))->silent    != NULL );
      assert( ((sdmodel*)(mymodel->model_pt))->get_class != NULL );

      /* generate sequences */
      printf("generating sequences: ...");
      
      
      /*
       * We do not specify length. A sequence will end when encounter a final state.
       */
      /***
      my_output=sdmodel_generate_sequences((sdmodel*)mymodel->model_pt,
					   1, 20,
					   5, 20);
      printf("Done\n");
      sequence_print(stdout, my_output);
      ******/

      /*** cpp_topo_ordering((sdmodel*)mymodel->model_pt); ***/
    }

  /*** testing viterbi ****/

  vipath=sdviterbi( (sdmodel*)mymodel->model_pt, seq, seq_len, &log_p);
  fprintf(stderr, "\n\n Viterbi path:");
  for(i=0; i < ((sdmodel*)mymodel->model_pt)->N*seq_len; i++) 
    fprintf(stderr, " %d ", vipath[i]);
  fprintf(stderr, "\n\n Viterbi log(prob) = %5.5f\n", log_p );

  /* reestimation */
  /*fprintf(stdout,"reestimating with Baum-Welch-algorithm...");
    reestimate_baum_welch((model*)mymodel->model_pt, my_output);
  */

  /*****
  double log_p;
  for(i = 0; i < num_sequences; i++)
    {
      viterbi(mymodel->model_pt, my_output->seq[i], my_output->seq_len[i], &log_p);
      printf( " Sequence %2d: viterbi prob = %5.5f\n", i, log_p );
    }
  ***/

  /* print the result */
  /****
  fprintf(stdout,"Done\nthe result is:\n");
  model_print(stdout,(model*)mymodel->model_pt);
  ****/

  /* free pointers */
  free(mymodel->model_pt);
  free((void*) mymodel);
  // free(my_output); 

  return 0;
}

