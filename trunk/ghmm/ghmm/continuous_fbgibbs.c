#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <math.h>
#include <float.h>
#include "ghmm.h"
#include "mprintf.h"
#include "mes.h"
#include "matrix.h"
#include "randvar.h"
#include "ghmm_internals.h"

#include "continuous_fbgibbs.h"
/*----------------------------------------------------------------------------*/
/*                     modified forwards to get cdf pmats                     */
/*----------------------------------------------------------------------------*/

static double sfoba_stepforward_gibbs (ghmm_cstate * s, double *alpha_t, int osc,
                                 double b_omega, double* pmats, int N)
{
  int i, id, prv;
  double value = 0.0;
  prv = s->in_id[0];
  for (i = 0; i < s->in_states; i++) {
    id = s->in_id[i];
    pmats[id] = s->in_a[osc][i] * alpha_t[id];
    value += pmats[id];
    
    //fill in values of pmats previous id to current id
    for(; prv < id; prv++){
        pmats[prv+1] += pmats[prv];
    }
    prv = id;
  }
  for(prv+=1;prv<N;prv++){
      pmats[prv] += pmats[prv-1];
  }
  value *= b_omega;             /* b_omega outside the sum */
  return (value);
}                               /* sfoba_stepforward */


/*============================================================================*/
int ghmm_cmodel_forwardgibbs (ghmm_cmodel * smo, double *O, int T, double ***b,
                   double **alpha, double *scale, double *log_p, double***pmats)
{
# define CUR_PROC "ghmm_cmodel_forward"
  int res = -1;
  int i, t = 0, osc = 0;
  double c_t;
  int pos;

  /* T is length of sequence; divide by dimension to represent the number of time points */
  T /= smo->dim;
  /* calculate alpha and scale for t = 0 */
  if (b == NULL)
    sfoba_initforward(smo, alpha[0], O, scale, NULL);
  else
    sfoba_initforward(smo, alpha[0], O, scale, b[0]);
  if (scale[0] <= DBL_MIN) {
    /* means f(O[0], mue, u) << 0, first symbol very unlikely */
    /* GHMM_LOG(LCONVERTED, "scale[0] == 0.0!\n"); */
    goto STOP;
  }
  else {
    *log_p = -log (1 / scale[0]);

    if (smo->cos == 1) {
      osc = 0;
    }
    else {
      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        return (-1);
      }
      /* printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t); */
      osc = smo->class_change->get_class (smo, O, smo->class_change->k, t);
      if (osc >= smo->cos){
        printf("ERROR: get_class returned index %d but model has only %d classes !\n",osc,smo->cos);
        goto STOP;
      }

    }


    for (t = 1; t < T; t++) {
      scale[t] = 0.0;
      pos = t * smo->dim;
      /* b not calculated yet */
      if (b == NULL) {
        for (i = 0; i < smo->N; i++) {
          alpha[t][i] = sfoba_stepforward_gibbs(smo->s+i, alpha[t-1], osc,
                                          ghmm_cmodel_calc_b(smo->s+i, O+pos),
                                          pmats[t][i],smo->N);
          scale[t] += alpha[t][i];
        }
      }
      /* b precalculated */
      else {
        for (i = 0; i < smo->N; i++) {
          alpha[t][i] = sfoba_stepforward_gibbs (smo->s+i, alpha[t - 1], osc,
                                           b[t][i][smo->M], pmats[t][i], smo->N);
          scale[t] += alpha[t][i];
        }
      }
      if (scale[t] <= DBL_MIN) {        /* seq. can't be build */
        goto STOP;
        break;
      }
      c_t = 1 / scale[t];
      /* scale alpha */
      for (i = 0; i < smo->N; i++)
        alpha[t][i] *= c_t;
      /* summation of log(c[t]) for calculation of log( P(O|lambda) ) */
      *log_p -= log (c_t);

      if (smo->cos == 1) {
        osc = 0;
      }
      else {
        if (!smo->class_change->get_class) {
          printf ("ERROR: get_class not initialized\n");
          return (-1);
        }
        /* printf("1: cos = %d, k = %d, t = %d\n",smo->cos,smo->class_change->k,t); */
        osc = smo->class_change->get_class (smo, O, smo->class_change->k, t);
        if (osc >= smo->cos){
          printf("ERROR: get_class returned index %d but model has only %d classes !\n",osc,smo->cos);
          goto STOP;
        }		
      }

    }
  }
  /* log_p should not be smaller than value used for seqs. that 
     can't be build ???
     if (*log_p < (double)PENALTY_LOGP)
     *log_p = (double)PENALTY_LOGP;
   */
  return 0;
STOP:
  *log_p = (double) -DBL_MAX;
  return (res);
# undef CUR_PROC
}                               /* ghmm_cmodel_forward */
//====================================================================================
//==================================end forwards======================================
//====================================================================================

//normal data holds information to update a normal hmm
typedef struct normal_update {
    double **transitions;
    double *in_state;
    double *sample_mean;
    double *sample_var;
}normal_update;


void alloc_normal_data(ghmm_cmodel* mo, normal_update **update){
    #define CUR_PROC "allocNormalData"
    *update = (normal_update*)malloc(sizeof(normal_update));
    (*update)->transitions = ighmm_cmatrix_alloc(mo->N,mo->N);
    ARRAY_CALLOC((*update)->in_state, mo->N);
    ARRAY_CALLOC((*update)->sample_mean, mo->N);
    ARRAY_CALLOC((*update)->sample_var, mo->N);
STOP:
    return;//XXX error handle
#undef CUR_PROC
}

void clear_normal_data(ghmm_cmodel*mo, normal_update* update){
    int i,j;
    for(i=0;i<mo->N;i++){
        update->in_state[i] = 0;
        update->sample_mean[i] = 0;
        update->sample_var[i] = 0;
        for(j=0;j<mo->N;j++){
            update->transitions[i][j] = 0;
        }
    }
}


void free_normal_data(ghmm_cmodel* mo, normal_update *update){
#define CUR_PROC "free_normal_data"
    m_free(update->sample_var);
    m_free(update->sample_mean);
    m_free(update->in_state);
    ighmm_cmatrix_free(&update->transitions, mo->N);
#undef CUR_PROC
}

void get_normal_data(ghmm_cmodel* mo, int *states, double* O, int T, normal_update *update){
    int i;
    double tmp;
    for(i=0;i<T;i++){
        update->sample_mean[states[i]] += O[i];
        update->in_state[states[i]]++;
    }
    for(i=0;i<mo->N;i++){
        if(update->in_state[i] > 0){
            update->sample_mean[i]/=update->in_state[i];
        }
    }

    for(i=0;i<T-1;i++)
        update->transitions[states[i]][states[i+1]]++;

    for(i=0;i<T;i++){
        tmp = O[i] - update->sample_mean[states[i]];
        update->sample_var[states[i]] += tmp*tmp;
    }
}

void update_normal(ghmm_cmodel *mo, double **pA, normal_hyper **hyperparams, double *pPi, normal_update *update){
    int i,k;
    double mue, nu, a, b;
    double tmp, sum;
    normal_hyper *p;
    //B
    for(i = 0; i < mo->N; i++){
        p = hyperparams[i];
        tmp = update->sample_mean[i] - p->mue;
        sum = p->mue + p->nu;

        mue = ( p->mue*p->nu + update->in_state[i]*update->sample_mean[i] ) / sum; 
        nu = sum;
        a = p->a + update->in_state[i]/2;
        b = p->b + .5*update->sample_var[i] + ( p->mue * p->nu * tmp * tmp) / (2*sum);
        mo->s[i].e->mean.val = ighmm_rand_normal(mue, nu, 0);
        mo->s[i].e->variance.val = ighmm_rand_gamma(a, b, 0);  
        printf("a %f, b %f, var %f, mue %f\n", a, b, mo->s[i].e->variance.val, mo->s[i].e->mean.val);
    } 

    //A
    for(i = 0; i < mo->N; i++){
        update->in_state[i] += pPi[i];
        for(k = 0; k < mo->N; k++){
            update->transitions[i][k] += pA[i][k];
        }
    }
    double tmp_n[mo->N];
    for(i=0;i<mo->N;i++){
        ighmm_rand_dirichlet(0, mo->N, update->transitions[i], tmp_n);
        for(k = 0; k < mo->N; k++){
            ghmm_cmodel_set_transition(mo, i, k, 0,tmp_n[k]);
        }        
    }

    //Pi
    ighmm_rand_dirichlet(0, mo->N, update->in_state, tmp_n);
    for(k=0;k<mo->N;k++){
        mo->s[k].pi = tmp_n[k];
    }
}

void ghmm_cmodel_fbgibbstep (ghmm_cmodel * mo, double *O, int len,int *Q, double** alpha, double***pmats){
    int i,j,k;
    for(i = 0; i < len; i++){
        for(j = 0; j < mo->N; j++){
            alpha[i][j] = 0;
            for(k = 0; k < mo->N; k++){
               pmats[i][j][k] = 0;
            }
        }
    }

    double scale[len];
    double logP;
    ghmm_cmodel_forwardgibbs(mo, O, len, NULL, alpha, scale, &logP, pmats);
    sampleStatePath(mo->N, alpha[len-1], pmats, len, Q);
}

int* ghmm_cmodel_fbgibbs (ghmm_cmodel * mo, ghmm_cseq* seq,
        double **pA, normal_hyper **pB, double *pPi, int burnIn){
#define CUR_PROC "ghmm_cmodel_fbgibbs"
    double **alpha = ighmm_cmatrix_alloc(seq->seq_len[0],mo->N);
    double ***pmats = ighmm_cmatrix_3d_alloc(seq->seq_len[0], mo->N, mo->N);
    int *Q;  
    ARRAY_CALLOC(Q, seq->seq_len[0]);
    normal_update *update;
    alloc_normal_data(mo, &update);
    for(; burnIn > 0; burnIn--){
        ghmm_cmodel_fbgibbstep(mo,seq->seq[0],seq->seq_len[0], Q, alpha, pmats);//XXX only using seq 1
        get_normal_data(mo, Q, seq->seq[0], seq->seq_len[0], update);
        update_normal(mo, pA, pB, pPi, update);
        clear_normal_data(mo, update);
    }
    free_normal_data(mo, update);
    ighmm_cmatrix_free(&alpha, seq->seq_len[0]);
    ighmm_cmatrix_3d_free(&pmats, seq->seq_len[0],mo->N);
    return Q;
STOP:
    return NULL; //XXX error handle
#undef CUR_PROC
}
