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
 * 2. Check the model
 * 3. Generate 100 sequences of length 100
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

  /* Important! initialise rng  */
  gsl_rng_init();

  mymodel = graphmldoc_cwrapper("gene_hmm.gml");

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
  my_output=model_generate_sequences((model*)mymodel->model_pt, /* model */
				     0,   /* random seed */
				     100, /* length of each sequence */
				     100); /* number of sequences */
  printf("Done\n");
  sequence_print(stdout, my_output);

  /* slight change of emission probabilities in state 0 */
  state_0_pt->b[0] = 0.2;
  state_0_pt->b[1] = 0.3;
  state_0_pt->b[2] = 0.4;
  state_0_pt->b[3] = 0.1;  
  state_0_pt->b[4] = 0.0;

  /* reestimation */
  fprintf(stdout,"reestimating with Baum-Welch-algorithm...");
  reestimate_baum_welch((model*)mymodel->model_pt, my_output);

  /* print the result */
  fprintf(stdout,"Done\nthe result is:\n");
  model_print(stdout,(model*)mymodel->model_pt);

  /* free pointers */
  free(mymodel->model_pt);
  free((void*) mymodel);
  free(my_output);

  return 0;
}

