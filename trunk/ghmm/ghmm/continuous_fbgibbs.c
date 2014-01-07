/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/model.c
*       Authors:  Benhard Knab, Bernd Wichern, Benjamin Georgi, Alexander Schliep
*
*       Copyright (C) 1998-2004 Alexander Schliep
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*       Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik,
*                               Berlin
*
*       Contact: schliep@ghmm.org
*
*       This library is free software; you can redistribute it and/or
*       modify it under the terms of the GNU Library General Public
*       License as published by the Free Software Foundation; either
*       version 2 of the License, or (at your option) any later version.
*
*       This library is distributed in the hope that it will be useful,
*       but WITHOUT ANY WARRANTY; without even the implied warranty of
*       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*       Library General Public License for more details.
*
*       You should have received a copy of the GNU Library General Public
*       License along with this library; if not, write to the Free
*       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
*
*
*       This file is version $Revision: 2304 $
*                       from $Date: 2013-05-31 13:48:13 -0400 (Fri, 31 May 2013) $
*             last change by $Author: ejb177 $.
*
*******************************************************************************/

#ifndef GHMM_CONTINUOUS_FBGIBBS
#define GHMM_CONTINUOUS_FBGIBBS
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
#include "rng.h"
#include "sfoba.h"
#include "continuous_fbgibbs.h"
#include "smodel.h"

//data for each state for posterior
typedef struct sample_emission_data{
    int emitted;
    union{
        double val;
        double *vec;
    }mean;
    union{
        double val;
        double *vec;
        double **mat;
    }variance;
    double a;
    double b;
} sample_emission_data;

//data for model posterior
typedef struct ghmm_sample_data{
    double **transition;
    sample_emission_data **state_data;  //[state][mixture] rename to emission_data
}ghmm_sample_data;
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
#undef CUR_PROC
}                               /* ghmm_cmodel_forward */
//====================================================================================
//==================================end forwards======================================
//====================================================================================



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



int ghmm_alloc_sample_data(ghmm_bayes_hmm *mo, ghmm_sample_data *data){
#define CUR_PROC "ghmm_alloc_sample_data"
//XXX must do alloc matrices for dim >1
    int i;
    data->transition = ighmm_cmatrix_alloc(mo->N, mo->N);
    ARRAY_MALLOC(data->state_data, mo->N);
    for(i = 0; i < mo->N; i++){
        ARRAY_MALLOC(data->state_data[i], mo->M[i]);
        /*for(i = 0; i < mo->M[i]; i++){//only needed for dim >1
            ghmm_alloc_emission_data(data->state_data[i][j], ghmm_bayes_hmm->params[i][j])
        }*/
    }
    return 0;
STOP:
    return -1;
#undef CUR_PROC
}
 
void ghmm_clear_emission_data(sample_emission_data *data){
    data->emitted = 0;
    data->mean.val = 0;
    data->variance.val = 0;
    data->a = 0;
    data->b = 0;
}          

void ghmm_clear_sample_data(ghmm_sample_data * data, ghmm_bayes_hmm *bayes){
    int i, j;
    for(i = 0; i < bayes->N; i++){
        for(j = 0; j < bayes->M[i]; j++){
            ghmm_clear_emission_data(&data->state_data[i][j]);
        }
        for(j=0; j<bayes->N; j++){
            data->transition[i][j] = 0;
        }
    }
}
       
void ghmm_get_emission_data_first_pass(sample_emission_data *data, ghmm_density_t type,
        double *observation){
    switch(type){
        case(normal):
            data->mean.val += *observation;
            data->emitted++;
        default://not supported
            return;
    }
}

void ghmm_get_emission_data_second_pass(sample_emission_data *data, ghmm_density_t type,
        double *observation){
    double tmp;
    switch(type){
        case(normal)://divide by emitted before 2 pass
            tmp = *observation - data->mean.val;
            data->variance.val += tmp*tmp;
        default://not supported
            return;
    }
}

void ghmm_get_sample_data(ghmm_sample_data *data, ghmm_bayes_hmm *bayes,int *Q, double *O, int T){
    int i;
    for(i=0; i<T-1; i++){
        data->transition[Q[i]][Q[i+1]]++;
        ghmm_get_emission_data_first_pass(&(data->state_data[Q[i]][0]),
                bayes->params[Q[i]][0].type, O+i);
    }
    //when adding other types might need to make this a function
    ghmm_get_emission_data_first_pass(&data->state_data[Q[T-1]][0],
            bayes->params[Q[i]][0].type, &O[i]);
    for(i=0; i< bayes->N; i++){
        if(data->state_data[i][0].emitted>0)
            data->state_data[i][0].mean.val /= data->state_data[i][0].emitted;
    }
    for(i=0;i<T;i++){
        ghmm_get_emission_data_second_pass(&data->state_data[Q[i]][0],
                bayes->params[Q[i]][0].type, &O[i]);
    }
}

   
/* using data colected in sample_emission_data sample from posterior distribution*/
void ghmm_update_emission(sample_emission_data *data, ghmm_hyperparameters *params,
        ghmm_c_emission *emission){
    switch(params->type){
        case(normal):
            {
                double mean, var, a, b;
                double tmp;
                //var
                var = params->emission[0].variance.val + data->emitted;
                //mean
                mean = params->emission[0].variance.val * params->emission[0].mean.val;
                mean += data->emitted*data->mean.val;
                mean /= var;
                //a
                a = params->emission[1].min + data->emitted/2;
                
                //b
                tmp = data->mean.val - params->emission[0].mean.val;
                b = params->emission[1].max + data->variance.val;
                b += (data->emitted*params->emission[0].variance.val*tmp*tmp)/(2*var);

                // sample from posterior hyperparameters
                emission->mean.val = ighmm_rand_normal(mean, (b/var) /
                        (params->emission[0].variance.val+data->emitted),0);
                tmp = 1/ighmm_rand_gamma(a, 1/b, 0);
                emission->variance.val = tmp;
            }
        default:
            return;
    }
}

void ghmm_update_model(ghmm_cmodel *mo, ghmm_bayes_hmm *bayes, ghmm_sample_data *data){
    int i, k;
    double tmp_n[mo->N];
    double tmp2_n[mo->N];
    //emission
    for(i=0; i<bayes->N; i++){
        for(k=0; k<bayes->M[k]; k++){
            ghmm_update_emission(&data->state_data[i][k], &bayes->params[i][k],&mo->s[i].e[k]);
        }
    }
    //add prior for A, Pi
    for(i = 0; i < bayes->N; i++){
        tmp_n[i] = 0;
        for(k=0;k<bayes->M[i];k++){
            tmp_n[i] += data->state_data[i][k].emitted; 
        }
        for(k = 0; k < bayes->N; k++){
            data->transition[i][k] += bayes->A[i][k];
            //printf("data %f\n", data->transition[i][k]);
        }
    }

    //Pi
    ighmm_rand_dirichlet(0, mo->N, tmp_n, tmp2_n);
    for(k=0;k<mo->N;k++){
        mo->s[k].pi = tmp2_n[k];
    }

    //A
    for(i=0;i<mo->N;i++){
        ighmm_rand_dirichlet(0, mo->N, data->transition[i], tmp_n);
        for(k = 0; k < mo->N; k++){
            ghmm_cmodel_set_transition(mo, i, k, 0, tmp_n[k]);
        }
    }
}


//only uses first sequence
int* ghmm_bayes_hmm_fbgibbs(ghmm_bayes_hmm *bayes, ghmm_cmodel *mo, ghmm_cseq* seq,
         int burnIn, int seed){
#define CUR_PROC "ghmm_cmodel_fbgibbs"
    //XXX seed
    GHMM_RNG_SET (RNG, seed);
    double **alpha = ighmm_cmatrix_alloc(seq->seq_len[0],mo->N);
    double ***pmats = ighmm_cmatrix_3d_alloc(seq->seq_len[0], mo->N, mo->N);
    int *Q; 
    ARRAY_CALLOC(Q, seq->seq_len[0]);
    ghmm_sample_data data;
    ghmm_alloc_sample_data(bayes, &data);
    ghmm_clear_sample_data(&data, bayes);//XXX swap parameter 
    for(; burnIn > 0; burnIn--){
        //XXX only using seq 0
        ghmm_cmodel_fbgibbstep(mo,seq->seq[0],seq->seq_len[0], Q, alpha, pmats);
        ghmm_get_sample_data(&data, bayes, Q, seq->seq[0], seq->seq_len[0]); 
        ghmm_update_model(mo, bayes, &data);
        ghmm_clear_sample_data(&data, bayes);
    }
    ighmm_cmatrix_free(&alpha, seq->seq_len[0]);
    ighmm_cmatrix_3d_free(&pmats, seq->seq_len[0],mo->N);
    return Q;
STOP:
    return NULL; //XXX error handle
#undef CUR_PROC
}
#endif
