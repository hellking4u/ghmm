
#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include <ghmm/matrix.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/model.h>
#include <ghmm/viterbi.h>
#include <ghmm/foba.h>
#include <ghmm/reestimate.h>
#include <ghmm/obsolete.h>

void uniformPsuedoCount(double **pA, double **pB, double *pPi, int N, int M){
        long i, j;
        
        for(i=0;i<N;i++){
            pPi[i] = 1;
            for(j=0;j<N;j++){
                pA[i][j] = 1;
            }
            for(j = 0; j < M; j++){
               pB[i][j] = 1;
            }
        }
    }
int* getList(char* fileName, int count){
 FILE *fp;
    int c;
    int* list = malloc(sizeof(int)*count);
    if (!(fp = fopen(fileName, "rt"))) {
        perror(fileName);
        exit(1);
    }
    int i = 0;
    while ((c = fgetc(fp)) != EOF && i < count ) {
      if(c == 32){
          c = 26;
      }
      else{
          c -= 97;
      }
      //printf("%d ", c);
      list[i] = c;
      i++;
    }
    fclose(fp);
    return list;
}






int main() {
  int result;
  /* Important! initialise rng  */
  ghmm_rng_init();

  int count = 140400;
  int* charList = getList("words.txt", count);

  ghmm_dmodel my_model;
  my_model.name = NULL;
  my_model.model_type = GHMM_kDiscreteHMM;
  my_model.silent = NULL; 
  my_model.maxorder = 0; 
  my_model.emission_history = 0;
  my_model.tied_to = NULL;
  my_model.order = NULL;
  my_model.bp = NULL;
  my_model.background_id = NULL;
  my_model.topo_order =NULL; 
  my_model.topo_order_length = 0; 
  my_model.pow_lookup = NULL; 
  my_model.label = NULL; 
  my_model.label_alphabet = NULL;
  my_model.alphabet = NULL; 


  ghmm_dstate model_states[2];
  double symbols_vowel_state[27]={0.03906, 0.03537, 0.03537, 0.03909, 0.03583, 0.03630, 0.04048, 0.03537,0.03816, 0.03909, 0.03490, 0.03723, 0.03537, 0.03909, 0.03397, 0.03397, 0.03816, 0.03676, 0.04048, 0.03443, 0.03537, 0.03955, 0.03816, 0.03723, 0.03769, 0.03955, 0.03688};
  double trans_prob_vowel_state[2]={0.47, 0.53};
  double trans_prob_vowel_state_rev[2]={0.47, 0.53};
  int trans_id_vowel_state[2]={0,1};
  double symbols_consonant_state[27]={0.03732, 0.03408, 0.03455, 0.03828, 0.03782, 0.03922, 0.03688, 0.03408, 0.03875, 0.04062, 0.03735, 0.03968, 0.03548, 0.03735, 0.04062, 0.03595, 0.03641, 0.03408, 0.04062, 0.03548, 0.03922, 0.04062, 0.03455, 0.03595, 0.03408, 0.03408, 0.03397};
  double trans_prob_consonant_state[2]={0.51,0.49};
  double trans_prob_consonant_state_rev[2]={0.51,0.49};
  int trans_id_consonant_state[2]={0,1};
  ghmm_dseq *my_output;
  double log_p_viterbi, log_p_forward;
  double **forward_alpha;
  double forward_scale[count];
  int *viterbi_path;
  int i, pathlen;
  /* flags indicating whether a state is silent */
  int silent_array[2] =  {0,0}; 

  my_model.model_type = 0;
  /* initialise vowel state */
  model_states[0].pi = 0.49;
  model_states[0].b=symbols_vowel_state;
  model_states[0].out_states=2;
  model_states[0].out_a=trans_prob_vowel_state;
  model_states[0].out_id=trans_id_vowel_state;
  model_states[0].in_states=2;
  model_states[0].in_id=trans_id_vowel_state;
  model_states[0].in_a=trans_prob_vowel_state_rev;
  model_states[0].fix=0;

  /* initialise consonant state */
  model_states[1].pi = 0.51;
  model_states[1].b=symbols_consonant_state;
  model_states[1].out_states=2;
  model_states[1].out_id=trans_id_consonant_state;
  model_states[1].out_a=trans_prob_consonant_state;
  model_states[1].in_states=2;
  model_states[1].in_id=trans_id_consonant_state;
  model_states[1].in_a=trans_prob_consonant_state_rev;
  model_states[1].fix=0;

  /* initialise model */
  my_model.N=2;
  my_model.M=27;
  my_model.s=model_states;
  my_model.prior=-1;
  my_model.silent = silent_array;
  
  fprintf(stdout,"transition matrix:\n");
  ghmm_dmodel_A_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  ghmm_dmodel_B_print(stdout,&my_model,""," ","\n");

  my_output = ghmm_dseq_calloc(1);
  my_output->seq[0] = charList;  
  my_output->seq_len[0] = count;

  //ghmm_dseq_print(my_output, stdout);


//====================tests for fbgibbs==================================================
  printf("fbgibbs \n");
  ghmm_dmodel* mo = ghmm_dmodel_copy(&my_model);
  
  double **pA = ighmm_cmatrix_alloc(my_model.N, my_model.N);
  double **pB = ighmm_cmatrix_alloc(my_model.N, my_model.M);
  double pPi[my_model.N];
  uniformPsuedoCount(pA, pB, pPi, my_model.N, my_model.M);
  int z;
  int iter = 100;
  int* paths[iter];
  for(z = 0; z<iter; z++)
    paths[z] = calloc(my_output->seq_len[0], sizeof(int));
  ghmm_dmodel* new_mos[iter];
  new_mos[0] = ghmm_dmodel_copy(&my_model);
  ghmm_dmodel_fbgibbstep(new_mos[0], 0, my_output->seq[0], my_output->seq_len[0], pA, pB, pPi, paths[0]);
  for(z = 1; z < iter; z++){
    new_mos[z] = ghmm_dmodel_copy(new_mos[z-1]);
    ghmm_dmodel_fbgibbstep( new_mos[z], 0, my_output->seq[0], my_output->seq_len[0],  pA, pB, pPi, paths[z]);
  }
  //for(z = 0; z < my_output->seq_len[0]; z++){
    //printf("%d ", paths[iter-1][z]);
  //}
  printf("viterbi prob mcmc%f \n", ghmm_dmodel_viterbi_logp(new_mos[iter-1], my_output->seq[0], my_output->seq_len[0], paths[iter-1]));
  printf("likelihood mcmc%f \n", ghmm_dmodel_likelihood(new_mos[iter-1], my_output));

  ghmm_dmodel_A_print(stdout,new_mos[iter-1],""," ","\n");
  ghmm_dmodel_B_print(stdout,new_mos[iter-1],""," ","\n");

//=====================end test fbgibbs================================================


  
 
//=====================viterbi/em=====================================================
  printf("Em/viterni\n\n");
  viterbi_path = ghmm_dmodel_viterbi(&my_model, my_output->seq[0],
				my_output->seq_len[0],&pathlen, &log_p_viterbi);
  if (viterbi_path==NULL)
    {fprintf(stderr,"viterbi failed!"); return 1;}
  
  //for(i=0;i<my_output->seq_len[0];i++){
    //printf(" %d, ", viterbi_path[i]);
  //}
  printf("\n");
  fprintf(stdout,
	  "no training (viterbi algorithm): %f\n",
	  log_p_viterbi);

  /* allocate matrix for forward algorithm */
  forward_alpha=ighmm_cmatrix_stat_alloc(count,2);
  if (forward_alpha==NULL)
    {
      fprintf(stderr,"\n could not alloc forward_alpha matrix\n");
      return 1;
    }

  /* run ghmm_dmodel_forward */
  if (ghmm_dmodel_forward(&my_model,
		   my_output->seq[0],
		   my_output->seq_len[0],
		   forward_alpha,
		   forward_scale,
		   &log_p_forward))
    {
      fprintf(stderr,"ghmm_dmodel_logp failed!");
      ighmm_cmatrix_stat_free(&forward_alpha);
      return 1;
    }
  ghmm_dmodel_baum_welch_nstep(&my_model, my_output, 100000, 0.0000001);
  viterbi_path = ghmm_dmodel_viterbi(&my_model, my_output->seq[0],
				my_output->seq_len[0],&pathlen, &log_p_viterbi);
  //print
  ghmm_dmodel_A_print(stdout,&my_model,""," ","\n");
  ghmm_dmodel_B_print(stdout,&my_model,""," ","\n");
  fprintf(stdout,
	  "(viterbi algorithm): %f\n",
	  log_p_viterbi);
  printf("likelihood %f \n", ghmm_dmodel_likelihood(&my_model, my_output));
//==================================================================================


  /* clean up */
  ghmm_dseq_free(&my_output);
  free(viterbi_path);
  ighmm_cmatrix_stat_free(&forward_alpha);
  for(z = 0; z < iter; z++){
    free(paths[z]);
    ghmm_dmodel_free(&new_mos[z]);
  }
  ighmm_cmatrix_free(&pA, my_model.N);
  ighmm_cmatrix_free(&pB, my_model.N);
  return 0;
}
