/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <ghmm/matrix.h>
#include <ghmm/rng.h>
#include <ghmm/smodel.h>

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_continuous()
{
  sstate single_state;
  smodel my_model;

  double trans_prob_single_state[]={1.0};
  double trans_prob_single_state_rev[]={1.0};
  double *trans_prob_single_state_array;
  double *trans_prob_single_state_rev_array;
  int trans_id_single_state[]={0};
  double c[]={1.0};
  double mue[]={0.0};
  double u[]={1.0};
  sequence_d_t* my_output;

  /* initialise transition array */
  trans_prob_single_state_array=trans_prob_single_state;
  trans_prob_single_state_rev_array=trans_prob_single_state_rev;

  /* initialise this state */
  single_state.pi = 1.0;
  single_state.out_states=1;
  single_state.out_a=&trans_prob_single_state_array;
  single_state.out_id=trans_id_single_state;
  single_state.in_states=1;
  single_state.in_id=trans_id_single_state;
  single_state.in_a=&trans_prob_single_state_rev_array;
  single_state.c=c;  /* weight of distribution */
  single_state.mue=mue; /* mean */
  single_state.u=u; /* variance */
  single_state.fix=0; /* training of output functions */

  /* initialise model */
  my_model.N=1; /* states */
  my_model.M=1; /* density functions per state */
  my_model.cos=1; /* class of states */
  my_model.density=0; /* normal distributions */
  my_model.prior=-1; /* a priori probability */
  my_model.s=&single_state; /* states array*/

#if 0
  /* print model */
  smodel_print(stdout,&my_model);
#endif

  /* generate sequences */
  my_output=smodel_generate_sequences(&my_model,
				      1,  /* random seed */
				      10, /* length of sequences */
				      10, /* sequences */
				      0,  /* label */
				      0  /* maximal sequence length 0: no limit*/
				      );
  /* print out sequences */
  sequence_d_print(stdout,    /* output file */
		   my_output, /* sequence */
		   0          /* do not truncate to integer*/
		   );


  {  
    FILE *my_file;
    char filename_buffer[]="/tmp/achimXXXXXX";
    int descriptor;
    int model_counter;
    smodel **model_array;

    descriptor=mkstemp(filename_buffer);
    /* write this model */
    my_file=fdopen(descriptor,"w+");
    fprintf(stdout,"printing model to file %s\n",filename_buffer);
    smodel_print(my_file,&my_model);
    (void)fseek(my_file, 0L, SEEK_SET);
    fclose(my_file);
    close(descriptor);

    /* read this model */
    model_array=smodel_read(filename_buffer,&model_counter);
  }
#if 0
  sequence_d_gnu_print(stdout,
		       my_output
		       );
#endif

  return 0;
}


int main()
{
  /* Important! initialise rng  */
  gsl_rng_init();

  return single_state_continuous();
}
