/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/two_states_three_symbols.c
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include <ghmm/vector.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/model.h>

/*
  model with two states and three symbols
  state 0 has symbol 0 and 1
  state 1 has symbol 2
*/
int my_model()
{
  model my_model;
  state model_states[2];
  double symbols_0_state[3]={0.5,0.5,0.0};
  double trans_prob_0_state[2]={0.5,0.5};
  int trans_id_0_state[2]={0,1};
  double symbols_1_state[3]={0.0,0.0,1.0};
  double trans_prob_1_state[2]={0.5,0.5};
  int trans_id_1_state[2]={0,1};
  sequence_t* my_output;

  /* initialise state 0 */
  model_states[0].pi = 0.5;
  model_states[0].b=symbols_0_state;
  model_states[0].out_states=2;
  model_states[0].out_a=trans_prob_0_state;
  model_states[0].out_id=trans_id_0_state;
  model_states[0].in_states=0;
  model_states[0].in_a=NULL;
  model_states[0].fix=0;

  /* initialise state 1 */
  model_states[1].pi = 0.5;
  model_states[1].b=symbols_1_state;
  model_states[1].out_states=2;
  model_states[1].out_a=trans_prob_1_state;
  model_states[1].out_id=trans_id_1_state;
  model_states[1].in_states=0;
  model_states[1].in_a=NULL;
  model_states[1].fix=0;

  /* initialise model */
  my_model.N=2; /* number of states */
  my_model.M=3; /* number of symbols */
  my_model.s=model_states;
  my_model.prior=-1;

  /* print out model parameters */
  fprintf(stdout,"two states and three symbols:\n");
  model_print(stdout,&my_model);

  /* generate sequences */
  fprintf(stdout,"generated sequence:\n");
  my_output=model_generate_sequences(&my_model,1,10,10);
  sequence_print(stdout,my_output);

#if 0
  /* randomize  transition vector */
  vector_random_values(trans_prob_0_state,2);
  vector_normalize(trans_prob_0_state,2);

  /* print again */
  fprintf(stdout,"randomized transition matrix:\n");
  model_A_print(stdout,&my_model,""," ","\n");
#else
  /* slight change of probabilities in state 0 */
  symbols_0_state[0] = 0.4;
  symbols_0_state[1] = 0.6;
  symbols_0_state[2] = 0.0;
#endif

  if (model_check(&my_model))
    fprintf(stderr,"model_check failed!\n");
  else
    reestimate_baum_welch(&my_model,my_output);

  model_print(stdout,&my_model);

  return 0;
}

int main()
{

  /* Important! initialise rng  */
  gsl_rng_init();

  return my_model();
}
