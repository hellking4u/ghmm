/*******************************************************************************
  author       : Achim Gaedke, Peter Pipenbacher
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <ghmm/rng.h>
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_State.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_continuous() {
  double c[] = {1.0};
  double mue[] = {0.0};
  double u[] = {1.0};
  GHMM_Sequences* my_output;

  /* initialise model */
  GHMM_ContinuousModel my_model(1,1,1,normal);

  /* initialise this state */
  my_model.getState(0)->setInitialProbability(1.0);
  my_model.setTransition(0,0,1.0);
  my_model.getCState(0)->c   = c;  /* weight of distribution */
  my_model.getCState(0)->mue = mue; /* mean */
  my_model.getCState(0)->u   = u; /* variance */

  /* generate sequences */
  my_output = my_model.generate_sequences(1,  /* random seed */
					  10, /* length of sequences */
					  10, /* sequences */
					  0,  /* label */
					  0  /* maximal sequence length 0: no limit*/
					  );
  /* print out sequences */
  my_output->print(stdout,    /* output file */
		   0          /* do not truncate to integer*/
		   );


  /* reproduce the sequence by saving and rereading the model */
  {  
    FILE *my_file;
    sequence_d_t* new_output;
    char filename_buffer[]="/tmp/chmm_test.XXXXXX";
    int descriptor;
    int model_counter;
    smodel **model_array;

    descriptor=mkstemp(filename_buffer);
    /* write this model */
    my_file=fdopen(descriptor,"w+");
    fprintf(stdout,"printing model to file %s\n",filename_buffer);
    my_model.print(my_file);
    (void)fseek(my_file, 0L, SEEK_SET);
    fclose(my_file);

    /* read this model */
    fprintf(stdout,"rereading model from file %s\n",filename_buffer);
    model_array=smodel_read(filename_buffer,&model_counter);

    /* generate sequences */
    fprintf(stdout,"generating sequences again\n");
    new_output=smodel_generate_sequences(model_array[0],
					 1,  /* random seed */
					 10, /* length of sequences */
					 10, /* sequences */
					 0,  /* label */
					 0   /* maximal sequence length 0: no limit*/
					 );

    sequence_d_print(stdout,     /* output file */
		     new_output, /* sequence */
		     0           /* do not truncate to integer*/
		     );
    /* free everything */
    close(descriptor);
    unlink(filename_buffer);
    sequence_d_free(&new_output);
    while(model_counter>0)
      {
	smodel_free(&(model_array[model_counter-1]));
	model_counter-=1;
      }
  }

  delete my_output;

  return 0;
}

int main()
{
  /* Important! initialise rng  */
  gsl_rng_init();

  return single_state_continuous();
}
