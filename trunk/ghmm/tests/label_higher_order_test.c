/*******************************************************************************
*  author       : Janne Grunau
*  filename     : ghmm/tests/label_higher_order_test.c
*  created      : DATE: 2004-05-07
*
*       This file is version $Revision$ 
*                       from $Date$
*             last change by $Author$
*******************************************************************************/

/* test_baumWelch
   generates a model to test c-functions with valgrind
*/

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#include <ghmm/ghmm.h>
#include <ghmm/model.h>
#include <ghmm/viterbi.h>
#include <ghmm/reestimate.h>
#include <ghmm/gradescent.h>
/*#include <ghmm/ghmm_internals.h>*/

#include <ghmm/rng.h>
#include <ghmm/sequence.h>

#define NR_SEQUENCES 50


void generateModel (model *mo, int noStates, unsigned int seed) {

  state * states;
  int h, i, j;
  double rnd;

  /* flags indicating whether a state is silent */
  int *silent_array;

  /* init the random number generator */
  srandom(seed);

  /*allocate memory for states and array of silent flags*/
  if (!(states = malloc (sizeof (state) * noStates)))
    {printf ("malloc failed in line %d", __LINE__); exit(1);}
  if (!(silent_array = calloc (sizeof (int), noStates))) 
    {printf ("malloc failed in line %d", __LINE__); exit(1);}

  mo->N = noStates;
  mo->M = 4;
  mo->maxorder = 1;
  mo->prior = -1;

  /* Model has Higher order Emissions and labeled states*/
  mo->model_type =  kLabeledStates;
  if (mo->maxorder>0)
    mo->model_type += kHigherOrderEmissions;
  /* kHigherOrderEmissions + kHasBackgroundDistributions*/

  /* allocate memory for pow look-up table and fill it */
  if (!(mo->pow_lookup = malloc (sizeof (int) * (mo->maxorder+2)))) 
    {printf ("malloc failed in line %d", __LINE__); exit(1);}
  
  mo->pow_lookup[0] = 1;
  for (i=1; i<mo->maxorder+2; i++)
    mo->pow_lookup[i] =  mo->pow_lookup[i-1] * mo->M;

  /*initialize states*/
  for (i=0; i < mo->N; i++) {
    states[i].pi = (0==i ? 1.0:0.0);
    states[i].fix = 0;
    states[i].label = i%4;
    states[i].order = i%2;
    states[i].out_states = 2;
    states[i].in_states = 2;

    /* allocate memory for the a, the out- and incoming States and b array for higher emmission order states*/
    states[i].b      = malloc(sizeof(double) * pow(mo->M, (states[i].order+1) ));
    states[i].out_id = malloc(sizeof(int)*states[i].out_states);
    states[i].in_id  = malloc(sizeof(int)*states[i].in_states);
    states[i].out_a  = malloc(sizeof(double)*states[i].out_states);
    states[i].in_a   = malloc(sizeof(double)*states[i].in_states);

    for (h=0; h<pow(mo->M,states[i].order); h++) {
      rnd = random()/(double)RAND_MAX;
      for (j=h*mo->M; j<(h*mo->M+mo->M); j++){
	states[i].b[j] = ( (0==((i+j)%mo->M)) ? rnd : (1-rnd) / (mo->M-1));
      }
    }

    if ((mo->N-1)==i) {
      states[i].out_id[0] = 0;
      states[i].out_id[1] = i;
    }
    else {
      states[i].out_id[0] = i;
      states[i].out_id[1] = i+1;
    }

    if (0==i) {
      states[i].in_id[0]  = i;
      states[i].in_id[1]  = mo->N-1;
    }
    else {
      states[i].in_id[1]  = i-1;
      states[i].in_id[0]  = i;
    }

    states[i].out_a[0] = 0.5;
    states[i].out_a[1] = 0.5;
    states[i].in_a[0]  = 0.5;
    states[i].in_a[1]  = 0.5;

#ifdef DEBUG
    printf("State %d goto    : %d, %d\n", i, states[i].out_id[0], states[i].out_id[1]);
    printf("State %d comefrom: %d, %d\n", i, states[i].in_id[0],  states[i].in_id[1]);
    printf("State %d goto    : %g, %g\n", i, states[i].out_a[0], states[i].out_a[1]);
    printf("State %d comefrom: %g, %g\n", i, states[i].in_a[0],  states[i].in_a[1]);
#endif
  }

  mo->s = states;
  mo->silent = silent_array;

#ifdef DEBUG
  for (i = 0; i < mo->N; i++) {
    printf("\n State %d:\n", i);
    for (j = 0; j < pow(mo->M,states[i].order+1); j++){
      printf("%g ",mo->s[i].b[j]);
    }
  }
#endif

  /* model_print(stdout, mo); */

}


void testBaumwelch(int seqlen){

  int  error;
 
  model * mo_gen = NULL;
  model * mo_time = NULL;
  model * mo_mem  = NULL;
  sequence_t * my_output = NULL;
   

  if (!(mo_gen = malloc (sizeof (model))))
    {printf ("malloc failed in line %d", __LINE__); exit(1);}
  if (!(mo_time = malloc (sizeof (model))))
    {printf ("malloc failed in line %d", __LINE__); exit(1);}
  if (!(mo_mem = malloc (sizeof (model))))
    {printf ("malloc failed in line %d", __LINE__); exit(1);}
      
  /* generate a model with variable number of states*/
  generateModel(mo_gen,  7, 92304);
  generateModel(mo_time, 5, 1704);
  generateModel(mo_mem,  5, 1704);

  /*generate a random sequence*/
  my_output = model_label_generate_sequences(mo_gen, 0, seqlen, NR_SEQUENCES, seqlen);
  /* model_add_noise(mo_time, .499, 0); */
  /* randomize the second */
  /* model_add_noise(mo_mem, .499, 0); */

  model_print(stdout, mo_time);
  model_print(stdout, mo_mem);
  printf("Distance between the two models: %g\n\n", model_distance(mo_time, mo_mem));

  /* shifting both models in diffrent directions */
  /* train the first */	 
  /*ghmm_dl_baum_welch(mo_time, my_output);*/
  error = ghmm_d_baum_welch(mo_time, my_output);

  /* train the second and hope they are equal */
  error = ghmm_d_baum_welch(mo_mem, my_output);

  model_print(stdout, mo_time);
  model_print(stdout, mo_mem);
  printf("Distance between the two trained models: %g\n", model_distance(mo_time, mo_mem));

  printf("Log-Likelyhood generating:    %g\n", model_likelihood (mo_gen, my_output));
  printf("Log-Likelyhood fb-Baum-Welch: %g\n", model_likelihood (mo_time, my_output));
  printf("Log-Likelyhood me-Baum-Welch: %g\n", model_likelihood (mo_mem, my_output));


  /* freeing memory */
  model_free(&mo_gen);
  model_free(&mo_time);
  model_free(&mo_mem);
  
  sequence_free(&my_output);
}

int main(){

  ghmm_rng_init();
  testBaumwelch(2700);

  return 0;
}
