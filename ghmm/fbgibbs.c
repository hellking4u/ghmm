#include "model.h"
#include "fbgibbs.h"
#include "foba.h"
#include "matrix.h"
#include "ghmm.h"
#include "mes.h"
#include "ghmm_internals.h"
#include "sequence.h"
#include "rng.h"
#include "randvar.h"
#include "obsolete.h"
//===========================================================================================
//=====================            sampleing           ======================================
//===========================================================================================
//returns a sample from a cdf
int sample(int seed, double* dist, int N){
    printf("sample\n\n");

    double total = dist[N-1];
    printf("total = %f\n", total);
    double rn = ighmm_rand_uniform_cont(seed, total, 0.0f);
    printf("rn = %f \n\n", rn);
    int i;
    if(rn <= dist[0])
      return 0;
    for(i = 1; i < N; i++){
       if(dist[i-1] <rn && rn <= dist[i])
         return i;
   }
   return N-1;
}
//currently sampling from the column of ptrs
double* getCDF(double*** pmats, int t, int state, int N){
#define CUR_PROC "getCDF"
    printf("Start getcdf t = %d, state = %d\n", t, state);
    int j;
    double* dist;
    ARRAY_MALLOC(dist, N);
    dist[0] = pmats[t][0][state];
    printf("dist 0: %f\n", dist[0]);
    for(j = 1; j < N; j++){
      dist[j] = dist[j-1] + pmats[t][j][state];
      printf("dist %d: %f\n", j, dist[j]);
    }
    return dist;
STOP:
  printf("Error allocating dist");
#undef CUR_PROC
}

int* sampleStatePath(int seed, ghmm_dmodel *mo, double *alpha, double ***pmats, int T){
#define CUR_PROC "sampleStatePath"
  printf("sampleStatePath\n\n");
  int* states;
  int i, j;
  double temp[mo->N];
  temp[0] = alpha[0];
  for(i = 1; i < mo->N; i++)
    temp[i] = temp[i-1] + alpha[i];
  printf("tmp1 = %f, tmp2 = %f\n", alpha[0], temp[1]);
  ARRAY_CALLOC(states, T);
  states[T-1] = sample(seed, temp, mo->N);
  printf("State T-1 = %d\n", states[T-1]);
  for(i = T-2; i >=0; i--){
    //get cdf to sample from pmats **could be optimized by makeing cdfs in pmats instead**
    double* dist = getCDF(pmats, i+1, states[i+1], mo->N);
    states[i] = sample(seed, dist, mo->N);
  }
  return states;
STOP:
  return NULL;
#undef CUR_PROC
}



//===========================================================================================
//========= Modified forward to get an array of matrices used for sampling ==================
//===========================================================================================
                                 /* ghmm_dmodel_forwardGibbs_init */
int ghmm_dmodel_forwardGibbs_init (ghmm_dmodel * mo, double *alpha_1, int symb, double *scale)
{
# define CUR_PROC "ghmm_dmodel_forwardGibbs_init"
  int i, j, id, in_id;
  double c_0;
  *scale = 0.0;

  /*printf(" *** foba_initforward\n");*/

  /*iterate over non-silent states*/
  /*printf(" *** iterate over non-silent states \n");*/
  for (i = 0; i < mo->N; i++) {
    if (!(mo->model_type & GHMM_kSilentStates) || !(mo->silent[i])) {
      /*no starting in states with order > 0 !!!*/
      if (!(mo->model_type & GHMM_kHigherOrderEmissions) || mo->order[i] == 0) {
        alpha_1[i] = mo->s[i].pi * mo->s[i].b[symb];
        /*printf("\nalpha1[%i]=%f\n",i,alpha_1[i]);*/

        *scale += alpha_1[i];
      }
      else {
        alpha_1[i] = 0;
      }
    }
  }
  /*iterate over silent states*/
  /*printf(" *** iterate over silent states \n");*/
  if (mo->model_type & GHMM_kSilentStates) {
    for (i = 0; i < mo->topo_order_length; i++) {
      id = mo->topo_order[i];
      alpha_1[id] = mo->s[id].pi;

      /*printf("\nsilent_start alpha1[%i]=%f\n",id,alpha_1[id]);*/

      for (j = 0; j < mo->s[id].in_states; j++) {
        in_id = mo->s[id].in_id[j];
        alpha_1[id] += mo->s[id].in_a[j] * alpha_1[in_id];

        /*printf("\n\tsilent_run alpha1[%i]=%f\n",id,alpha_1[id]);*/

      }
      *scale += alpha_1[id];
    }
  }

  /* printf("\nwo label: scale[0] = %f\n",scale[0]); */

  if (*scale >= GHMM_EPS_PREC) {
    c_0 = 1 / *scale;
    for (i = 0; i < mo->N; i++)
      alpha_1[i] *= c_0;
  }
  return (0);                   /* attention: scale might be 0 */
# undef CUR_PROC
}                               /* ghmm_dmodel_forwardGibbs_init */

/*----------------------------------------------------------------------------*/

double ghmm_dmodel_forwardGibbs_step (ghmm_dstate * s, double *alpha_t, const double b_symb, double*** pmats, int t, int j)
{
  int i, id;
  double value = 0.0;

  //if (b_symb < GHMM_EPS_PREC)
    //return 0.0;

  /*printf(" *** fobagibbs_stepforward\n");*/
  printf("%d\n", s->in_states);
  for (i = 0; i < s->in_states; i++) { 
    id = s->in_id[i];
    pmats[t][i][j] = s->in_a[i] * alpha_t[id];
    value += pmats[t][i][j];
    //pmats[t][i][j] *= b_symb;//pavel didnt do this ask why?.. it was in scott paper
    printf("    state %d, value %f, p_symb %f, pmats %f\n",id, value, b_symb, pmats[t][i][j]); 
  }
  value *= b_symb;
  return (value);

}                               /* ghmm_dmodel_forwardGibbs_step */

/*============================================================================*/

int ghmm_dmodel_forwardGibbs (ghmm_dmodel * mo, const int *O, int len, double **alpha, double*** pmats)
{
# define CUR_PROC "ghmm_dmodel_forwardGibbs"
  char * str;
  int i, t, id;
  int e_index;
  double c_t;
  double scale;


  if (mo->model_type & GHMM_kSilentStates)
    ghmm_dmodel_order_topological(mo);

  ghmm_dmodel_forwardGibbs_init (mo, alpha[0], O[0], &scale);

  if (scale < GHMM_EPS_PREC) {
    /* means: first symbol can't be generated by hmm */
    //printf("\nscale kleiner als eps (line_no: 123)\n");
    return -1;
  }
  else {
    for (t = 1; t < len; t++) {

      scale = 0.0;
      update_emission_history (mo, O[t - 1]);

      //printf("\n\nStep t=%i mit len=%i, O[i]=%i\n",t,len,O[t]);
      //printf("iterate over non-silent state\n"); 
      /* iterate over non-silent states */
      for (i = 0; i < mo->N; i++) {
        if (!(mo->model_type & GHMM_kSilentStates) || !(mo->silent[i])) {
          e_index = get_emission_index (mo, i, O[t], t);
          if (e_index != -1) {
            printf("pmat %d, %d \n", t, i);
            alpha[t][i] =
              ghmm_dmodel_forwardGibbs_step (&mo->s[i], alpha[t - 1], mo->s[i].b[e_index], pmats, t, i);
            scale += alpha[t][i];
          }
          else {
            alpha[t][i] = 0;
          }
        }
      }
      /* iterate over silent states */
      //printf("iterate over silent state\n"); 
      if (mo->model_type & GHMM_kSilentStates) {
        for (i = 0; i < mo->topo_order_length; i++) {
          /*printf("\nget id\n");*/
          id = mo->topo_order[i];
          /*printf("  akt_ state %d\n",id);*/
          /*printf("\nin stepforward\n");*/
          alpha[t][id] = ghmm_dmodel_forwardGibbs_step (&mo->s[id], alpha[t], 1, pmats, t, i);
          /*printf("\nnach stepforward\n");*/
          scale += alpha[t][id];
        }
      }

      if (scale < GHMM_EPS_PREC) {
        /* O-string  can't be generated by hmm */
        str = ighmm_mprintf(NULL, 0, "scale smaller than epsilon (%g < %g) in position %d. Can't generate symbol %d\n", scale, GHMM_EPS_PREC, t, O[t]);
	GHMM_LOG(LCONVERTED, str);
	m_free (str);
        return -1;
      }
      c_t = 1 / scale;
      for (i = 0; i < mo->N; i++) {
        alpha[t][i] *= c_t;
      }
    }
  }
  return 0;
# undef CUR_PROC
}                               

			/* ghmm_dmodel_forwardGibbs */




//======================== update parameters ===========================================
//given states, prior matrices, and model calculates new A,B,Pi
void update(int seed, int T, int *states, int* O, double **priorA, double **priorB, double *priorPi, ghmm_dmodel *mo, double **A, double **B, double *Pi){
  double transition[mo->N][mo->N];
  double obsinstate[mo->N]; 
  double obsinstatealpha[mo->N][mo->M]; 
    int i,k;
    //initialize
    for(i=0;i<mo->N;i++){
        obsinstate[i]  = 0;
        for(k=0;k<mo->N;k++){
            transition[i][k] = 0;
        }
        for(k=0;k<mo->M;k++){
            obsinstatealpha[i][k] = 0;
        }
    }
    //count states and transitions
    for(i=0;i<T;i++){        
        obsinstate[states[i]] += 1;
        obsinstatealpha[states[i]][O[i]] += 1;
    }
    for(i=0;i<T-1;i++){
        transition[states[i]][states[i+1]] += 1;
    }

    for(i=0;i<mo->N;i++){
        for(k=0;k<mo->M;k++){
            obsinstatealpha[i][k] += priorB[i][k];
        }
        for(k=0;k<mo->N;k++){
            transition[i][k] += priorA[i][k];
        }
        ighmm_rand_dirichlet(seed, mo->M, obsinstatealpha[i], B[i]);
        ighmm_rand_dirichlet(seed, mo->N, transition[i], A[i]);
    }

    for(k=0; k<mo->N; k++){
        obsinstate[k] += priorPi[k];
    }
    ighmm_rand_dirichlet(seed, mo->N, obsinstate, Pi);
}
//===========================fbgibbstep==================================================

void fbgibbstep (int seed, ghmm_dmodel * mo, int *O, int len, double **A, double **B, double *Pi, double **priorA, double **priorB, double *priorPi, int steps){
  //sample states from model
  //update parameters store in A, B, Pi 
  //return new parameters
  printf("fbgibbsStep \n\n");
  GHMM_RNG_SET (RNG, seed);
  //initilizations
  double **alpha = ighmm_dmatrix_alloc(len, mo->N);
  double ***pmats = ighmm_dmatrix_3d_alloc(len, mo->N, mo->N);
  int i,j,k;
  for(i = 0; i < len; i++){
    for(j = 0; j < mo->N; j++){
      for(k = 0; k < mo->N; k++){
        pmats[i][j][k] = 0;
      }
    }
  }
  ghmm_dmodel_forwardGibbs(mo, O, len, alpha, pmats);
  int *Q = sampleStatePath(seed, mo, alpha[len-1], pmats, len);
  printf("done samplestatepath\n\n");
  for(i = 0; i < len; i++){
    printf("%d ", Q[i]);
  }
  printf("\n");
  update(seed, len, Q, O, priorA, priorB, priorPi, mo, A, B, Pi);
  ighmm_cmatrix_print(stdout,A,mo->N,mo->N,""," ","\n\n\n");
  ighmm_cmatrix_print(stdout,B,mo->N,mo->N,""," ","\n");
}
