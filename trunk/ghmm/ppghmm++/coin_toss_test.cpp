/*******************************************************************************
  authors      : Achim Gaedke, Peter Pipenbacher
  filename     : ppghmm++/coin_toss_test.cpp
  created      : DATE: 2001-04-25
  $Id$
*******************************************************************************/

#include <stdio.h>
#include "ghmm/rng.h"
#include "ghmm/vector.h"
#include "ghmm/matrix.h"
#include "ppghmm++/GHMM.h"
#include "ppghmm++/GHMM_DiscreteModel.h"
#include "ppghmm++/GHMM_Sequences.h"


/* Simple model with one state and 2 symbols, like a coin toss */

int single_state_coin_toss() {
  double symbols_single_state[2]        ={0.5,0.5};
  double trans_prob_single_state[1]     = {1.0};
  double trans_prob_single_state_rev[1] = {1.0};
  int trans_id_single_state[1]          = {0};
  GHMM_Sequences* my_output = NULL;

  /* initialise model */
  GHMM_DiscreteModel my_model(1,2,-1);

  /* initialise this state */
  my_model.getState(0)->pi         = 1.0;
  my_model.getState(0)->b          = symbols_single_state;
  my_model.getState(0)->out_states = 1;
  my_model.getState(0)->out_a      = trans_prob_single_state;
  my_model.getState(0)->out_id     = trans_id_single_state;
  my_model.getState(0)->in_states  = 1;
  my_model.getState(0)->in_id      = trans_id_single_state;
  my_model.getState(0)->in_a       = trans_prob_single_state_rev;
  my_model.getState(0)->fix        = 1;

  fprintf(stdout,"transition matrix:\n");
  my_model.A_print(stdout,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  my_model.B_print(stdout,""," ","\n");

  my_output = my_model.generate_sequences(1,10,10);
  my_output->print(stdout);

  delete my_output;

  return 0;
}

/*
  another coin toss model
  flip between two states and output a different symbol
 */

int two_states_coin_toss() {
  double symbols_head_state[2]        = {1.0,0.0};
  double trans_prob_head_state[2]     = {0.5,0.5};
  double trans_prob_head_state_rev[2] = {0.5,0.5};
  int trans_id_head_state[2]={0,1};
  double symbols_tail_state[2]={0.0,1.0};
  double trans_prob_tail_state[2]={0.5,0.5};
  double trans_prob_tail_state_rev[2]={0.5,0.5};
  int trans_id_tail_state[2]={0,1};
  GHMM_Sequences* my_output;
  double log_p_viterbi, log_p_forward;
  double forward_scale[10];
  int* viterbi_path;
  double **forward_alpha;

  /* initialise model */
  GHMM_DiscreteModel my_model(2,2,-1);

  /* initialise head state */
  my_model.getState(0)->pi = 0.5;
  my_model.getState(0)->b=symbols_head_state;
  my_model.getState(0)->out_states=2;
  my_model.getState(0)->out_a=trans_prob_head_state;
  my_model.getState(0)->out_id=trans_id_head_state;
  my_model.getState(0)->in_states=2;
  my_model.getState(0)->in_id=trans_id_head_state;
  my_model.getState(0)->in_a=trans_prob_head_state_rev;
  my_model.getState(0)->fix=1;

  /* initialise tail state */
  my_model.getState(1)->pi = 0.5;
  my_model.getState(1)->b=symbols_tail_state;
  my_model.getState(1)->out_states=2;
  my_model.getState(1)->out_id=trans_id_tail_state;
  my_model.getState(1)->out_a=trans_prob_tail_state;
  my_model.getState(1)->in_states=2;
  my_model.getState(1)->in_id=trans_id_tail_state;
  my_model.getState(1)->in_a=trans_prob_tail_state_rev;
  my_model.getState(1)->fix=1;

  fprintf(stdout,"transition matrix:\n");
  my_model.A_print(stdout,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  my_model.B_print(stdout,""," ","\n");

  my_output = my_model.generate_sequences(1,10,10);
  my_output->print(stdout);

  /* try viterbi algorithm in a clear situation */
  viterbi_path = my_model.Viterbi(my_output,0,&log_p_viterbi);
  if (!viterbi_path) {
    fprintf(stderr,"viterbi failed!"); 
    return 1;
  }

  fprintf(stdout,"viterbi path:\n");
  vector_i_print(stdout,viterbi_path,10,""," ","\n");

  fprintf(stdout,
	  "probability of this sequence (viterbi algorithm): %f\n",
  	  log_p_viterbi);

  /* allocate matrix for forward algorithm */
  fprintf(stdout,"applying forward algorithm to the sequence...");
  forward_alpha = matrix_d_alloc(10,2);
  if (!forward_alpha) {
    fprintf(stderr,"\n could not alloc forward_alpha matrix\n");
    return 1;
  }

  /* run foba_forward */
  if (my_model.fobaForward(my_output,0,forward_alpha,forward_scale,&log_p_forward)) {
    fprintf(stderr,"foba_logp failed!");
    matrix_d_free(&forward_alpha,10);
    return 1;
  }
  
  /* alpha matrix */
  fprintf(stdout,"Done.\nalpha matrix from forward algorithm:\n");
  matrix_d_print(stdout,forward_alpha,10,2,""," ","\n");
  fprintf(stdout,"probability of this sequence (forward algorithm): %f\n",log_p_forward);
  
  /* clean up */
  matrix_d_free(&forward_alpha,10);
  delete my_output;

  return 0;
}


int main() {
  /* Important! initialise rng  */
  gsl_rng_init();

  if (single_state_coin_toss() || two_states_coin_toss())
    return 1;
  else
    return 0;
}
