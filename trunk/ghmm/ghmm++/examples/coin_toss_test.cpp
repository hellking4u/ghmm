/*******************************************************************************
  authors      : Achim Gaedke, Peter Pipenbacher
  filename     : ghmm++/examples/coin_toss_test.cpp
  created      : DATE: 2002-03-08
  $Id$

  __copyright__
*******************************************************************************/

#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

/* Simple model with one state and 2 symbols, like a coin toss */
int single_state_coin_toss() {
  GHMM_Sequences* my_output = NULL;

  /* initialise model with 1 state and 2 symbols */
  GHMM_DiscreteModel my_model(1,2);

  /* initialise this state */
  my_model.getState(0)->setInitialProbability(1.0);
  my_model.getState(0)->setOutputProbability(0,0.5);
  my_model.getState(0)->setOutputProbability(1,0.5);

  my_model.setTransition(0,0,1.0);

  if (my_model.check() == -1)
    return 1;

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
  GHMM_Sequences* my_output;
  double log_p_viterbi, log_p_forward;
  GHMM_IntVector* viterbi_path;
  GHMM_DoubleMatrix* forward_alpha = NULL;

  /* initialise model */
  GHMM_DiscreteModel my_model(2,2);

  /* initialise head state */
  my_model.getState(0)->setInitialProbability(0.5);
  my_model.getState(0)->setOutputProbability(0,1.0);

  /* initialise tail state */
  my_model.getState(1)->setInitialProbability(0.5);
  my_model.getState(1)->setOutputProbability(1,1.0);

  my_model.setTransition(0,0,0.5);
  my_model.setTransition(0,1,0.5);
  my_model.setTransition(1,0,0.5);
  my_model.setTransition(1,1,0.5);

  if (my_model.check() == -1)
    return 1;

  fprintf(stdout,"transition matrix:\n");
  my_model.A_print(stdout,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  my_model.B_print(stdout,""," ","\n");

  my_output = my_model.generate_sequences(1,10,100);
  my_output->print(stdout);

  /* try viterbi algorithm in a clear situation */
  viterbi_path = my_model.viterbi(my_output,0,&log_p_viterbi);
  if (!viterbi_path) {
    fprintf(stderr,"viterbi failed!"); 
    return 1;
  }

  fprintf(stdout,"viterbi path:\n");
  viterbi_path->print(stdout,""," ","\n");

  fprintf(stdout,
	  "probability of this sequence (viterbi algorithm): %f\n",
  	  log_p_viterbi);

  /* allocate matrix for forward algorithm */
  fprintf(stdout,"applying forward algorithm to the sequence...");

  /* run foba_forward */
  forward_alpha = my_model.fobaForward(my_output,0,NULL,&log_p_forward);
  
  if (!forward_alpha) {
    fprintf(stderr,"foba_logp failed!");
    return 1;
  }
  
  /* alpha matrix */
  fprintf(stdout,"Done.\nalpha matrix from forward algorithm:\n");
  forward_alpha->print(stdout,""," ","\n");
  fprintf(stdout,"probability of this sequence (forward algorithm): %f\n",log_p_forward);

  /* change initialise probabilities */
  my_model.getState(0)->setInitialProbability(0.1);
  my_model.getState(1)->setInitialProbability(0.9);

  if (my_model.check() == -1)
    return 1;

  my_model.print(stdout);

  my_model.reestimate_baum_welch(my_output);

  if (my_model.check() == -1)
    return 1;

  my_model.print(stdout);
  
  /* clean up */
  delete my_output;
  SAFE_DELETE(viterbi_path);
  SAFE_DELETE(forward_alpha);

  return 0;
}


int main() {
  /* Important! initialise rng  */
  GHMM_Toolkit::gsl_rng_init();

  int result = (single_state_coin_toss() || two_states_coin_toss());

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return result;
}
