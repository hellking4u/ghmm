/*******************************************************************************
  author       : Achim Gaedke
  filename     : ghmm/tests/coin_toss_test.c
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/model.h>

/*
  Simple model with one state and 2 symbols, like a coin toss
*/

int single_state_coin_toss()
{
  state single_state;
  model my_model;
  double symbols_single_state[2]={0.5,0.5};
  double trans_prob_single_state[1]={1.0};
  int trans_id_single_state[1]={0};
  sequence_t* my_output;

  /* initialise this state */
  single_state.pi = 1.0;
  single_state.b=symbols_single_state;
  single_state.out_states=1;
  single_state.out_a=trans_prob_single_state;
  single_state.out_id=trans_id_single_state;
  single_state.in_states=0;
  single_state.in_a=0x0;
  single_state.fix=1;

  /* initialise model */
  my_model.N=1;
  my_model.M=2;
  my_model.s=&single_state;
  my_model.prior=-1;

  fprintf(stdout,"transition matrix:\n");
  model_A_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  model_B_print(stdout,&my_model,""," ","\n");

  my_output=model_generate_sequences(&my_model,1,10,10);
  sequence_print(stdout,my_output);

  return 0;
}


/*
  another coin toss model
  flip between two states and output a different symbol
 */

int two_states_coin_toss()
{
  model my_model;
  state model_states[2];
  double symbols_head_state[2]={1.0,0.0};
  double trans_prob_head_state[2]={0.5,0.5};
  int trans_id_head_state[2]={0,1};
  double symbols_tail_state[2]={0.0,1.0};
  double trans_prob_tail_state[2]={0.5,0.5};
  int trans_id_tail_state[2]={0,1};
  sequence_t* my_output;

  /* initialise head state */
  model_states[0].pi = 0.5;
  model_states[0].b=symbols_head_state;
  model_states[0].out_states=2;
  model_states[0].out_a=trans_prob_head_state;
  model_states[0].out_id=trans_id_head_state;
  model_states[0].in_states=0;
  model_states[0].in_a=NULL;
  model_states[0].fix=1;

  /* initialise tail state */
  model_states[1].pi = 0.5;
  model_states[1].b=symbols_tail_state;
  model_states[1].out_states=2;
  model_states[1].out_a=trans_prob_tail_state;
  model_states[1].out_id=trans_id_tail_state;
  model_states[1].in_states=0;
  model_states[1].in_a=NULL;
  model_states[1].fix=1;

  /* initialise model */
  my_model.N=2;
  my_model.M=2;
  my_model.s=model_states;
  my_model.prior=-1;

  fprintf(stdout,"transition matrix:\n");
  model_A_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  model_B_print(stdout,&my_model,""," ","\n");

  my_output=model_generate_sequences(&my_model,1,10,10);
  sequence_print(stdout,my_output);

  return 0;
}

int main()
{
  /* Important! initialise rng  */
  gsl_rng_init();

  if (single_state_coin_toss() || two_states_coin_toss())
    return 1;
  else
    return 0;
}
