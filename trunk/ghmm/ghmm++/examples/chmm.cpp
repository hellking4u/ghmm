/*******************************************************************************
  author       : Achim Gaedke, Peter Pipenbacher
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_continuous() {
  double c   = 1.0;
  double mue = 0.0;
  double u   = 1.0;
  GHMM_Sequences* my_sequences;

  /* initialise model */
  GHMM_ContinuousModel my_model(1,1,1,normal);

  /* initialise this state */
  my_model.getState(0)->setInitialProbability(1.0);
  my_model.setTransition(0,0,1.0);
  my_model.getCState(0)->c[0]   = c;  /* weight of distribution */
  my_model.getCState(0)->mue[0] = mue; /* mean */
  my_model.getCState(0)->u[0]   = u; /* variance */

  if (my_model.check() == -1)
    return 1;

  /* generate sequences */
  my_sequences = my_model.generate_sequences(1,  /* random seed */
					     10, /* length of sequences */
					     10, /* sequences */
					     0   /* label */
					     );
  /* print out sequences */
  my_sequences->print(stdout,    /* output file */
		      0          /* do not truncate to integer*/
		      );


  /* reproduce the sequence by saving and rereading the model */
  FILE *my_file;
  GHMM_Sequences* new_sequences = NULL;
  char* filename="chmm.tmp";
  GHMM_ContinuousModel new_model;

  /* write this model */
  fprintf(stdout,"printing model to file %s\n",filename);
  my_file = fopen(filename,"w");
  my_model.print(my_file);
  fclose(my_file);

  /* read this model */
  fprintf(stdout,"rereading model from file %s\n",filename);
  new_model.read(filename);

  /* generate sequences */
  fprintf(stdout,"generating sequences again\n");
  new_sequences = new_model.generate_sequences(1,  /* random seed */
					       10, /* length of sequences */
					       10, /* sequences */
					       0   /* label */
					       );

  new_sequences->print(stdout,     /* output file */
		       0           /* do not truncate to integer*/
		       );
  /* free everything */
  
  delete new_sequences;
  delete my_sequences;

  return 0;
}

int main() {
  /* Important! initialise rng  */
  GHMM_Toolkit::gsl_rng_init();

  int result = single_state_continuous();

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return result;
}
