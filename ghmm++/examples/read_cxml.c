/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/model.h>
#include <ghmm++/GHMM_convertXMLtoC.h>

/* What this code does:

 * Testing if we can call our C++ wrapper from C
 * 1. Read an XML file in GraphML format
 * 2. Check the model (HACK!!! no checking because of silent states)
 * 3. Generate 100 sequences 
 * 4. Change the B matrix slightly and train with Baum-Welch 
 */
int main(int argc, char* argv[]) 
{
  /* Be careful!! : these are pointers to C-struct 
   */
  model_t   *mymodel;
  sequence_t*my_output;
  model     *cmodel;
  state     *state_0_pt;

  int       i;
  int      num_sequences=5;

  if (argc < 2)
    {
      fprintf(stderr, "Please input a filename. Usage: read_cxml file.gml\n");
      exit(-1);
    }

  /* Important! initialise rng  */
  gsl_rng_init();


  mymodel = graphmldoc_cwrapper(argv[1]);


  if (mymodel != NULL)
    {
      if (mymodel->model_id == DISCRETE)
	{
	  printf("transition matrix A:\n");
	  model_A_print(stdout, (model*) mymodel->model_pt, " ", ",", "\n");
	  printf("Observation matrix B:\n");
	  model_B_print(stdout, (model*) mymodel->model_pt, " ", ",", "\n");
	  cmodel = (model*)(mymodel->model_pt); /* the actual pointer to C struct model */
	  state_0_pt = &(cmodel->s[0]);         /* 1st state */
	}
    }

  /* generate sequences */
  printf("generating sequences: ...");

  /*
    We do not specify length. A sequence will end when encounter a final state.
   */

  my_output=model_generate_sequences((model*)mymodel->model_pt, /* model */
				     0, /* length of each sequence */
				     100); /* number of sequences */

  printf("Done\n");
  sequence_print(stdout, my_output);


  /* reestimation */
  /*fprintf(stdout,"reestimating with Baum-Welch-algorithm...");
    reestimate_baum_welch((model*)mymodel->model_pt, my_output);
  */

  double log_p;
  for(i = 0; i < num_sequences; i++)
    {
      viterbi(mymodel->model_pt, my_output->seq[i], my_output->seq_len[i], &log_p);
      printf( " Sequence %2d: viterbi prob = %5.5f\n", i, log_p );
    }

  /* print the result */
  fprintf(stdout,"Done\nthe result is:\n");
  model_print(stdout,(model*)mymodel->model_pt);

  /* free pointers */
  free(mymodel->model_pt);
  free((void*) mymodel);
  free(my_output);

  return 0;
}

