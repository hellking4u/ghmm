/* author       : Wasinee Rungsarityotin and Benjamin Georgi
 *  filename    : ghmmwrapper/gql.c
 *  created      : DATE: September, 2003
 *
 * $Id$
 */

/*
  __copyright__
*/

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/mes.h>
#include <ghmm/matrix.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/smodel.h>
#include <ghmm++/GHMM_convertXMLtoC.h>
#include <float.h>
#include <assert.h>


static double *log_p_pt;
static int seq_rank(const void *a1, const void *a2);
  
int seq_rank(const void *a1, const void *a2) {
  
  int i, j;

  i = *((int*)a1);
  j = *((int*)a2);

  if ( log_p_pt[i] == log_p_pt[j] ) return 0;
  else
    if ( log_p_pt[i] < log_p_pt[j] ) return -1;
    else
      return 1;
}

static int smodel_state_alloc(sstate *state,
			      int M,
			      int in_states,
			      int out_states,
			      int cos) {
# define CUR_PROC "smodel_state_alloc"
  int res = -1;
  if(!m_calloc(state->c, M)) {mes_proc(); goto STOP;}
  if(!m_calloc(state->mue, M)) {mes_proc(); goto STOP;}
  if(!m_calloc(state->u, M)) {mes_proc(); goto STOP;}
  if (out_states > 0) {
    if(!m_calloc(state->out_id, out_states)) {mes_proc(); goto STOP;}
    state->out_a = matrix_d_alloc(cos, out_states);
    if(!state->out_a) {mes_proc(); goto STOP;}
  }
  if (in_states > 0) {
    if(!m_calloc(state->in_id, in_states)) {mes_proc(); goto STOP;}
    state->in_a = matrix_d_alloc(cos, in_states);
    if(!state->in_a) {mes_proc(); goto STOP;}
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
}

smodel *smodel_alloc_fill(int N, int M, int cos, double prior, int density) {
  int i;
  smodel *smo=NULL;
  if ((smo = malloc(sizeof(smodel))) == NULL) {
    fprintf(stderr, "smodel_alloc_fill(1): out of memory\n");
    return NULL;
  }
  smo->M   = M;
  smo->N   = N;
  smo->cos = cos;
  smo->prior = prior;
  smo->density = density;  /* 0 = normal, 1 = normal_pos */
  if (!m_calloc(smo->s, smo->N)) {
    fprintf(stderr, "smodel_alloc_fill(1): out of memory\n");
    return NULL;
  }
  for(i=0; i < smo->N; i++) {
    smodel_state_alloc(&smo->s[i], smo->M, smo->N, smo->N, cos);
  }
  return smo;
}

void smodel_set_pivector(smodel *smo, int i, double piv) {
  if (smo->s != NULL) {
    smo->s[i].pi = piv;
  }
}

void smodel_set_fixvector(smodel *smo, int i, double fixv) {
  if (smo->s != NULL) {
    smo->s[i].fix = fixv;
  }
}

void smodel_set_transition(smodel *smo, int i, int j, int cos, double prob) {
  int in, out;
  if (cos >= smo->cos) {
    fprintf(stderr, "smodel_set_transition(cos): cos > state->cos\n");
    exit(-1);	
  }
  if (smo->s && smo->s[i].out_a && smo->s[j].in_a) {
    for(out=0; out < smo->s[i].out_states; out++) {
      if ( smo->s[i].out_id[out] == j ) {
	smo->s[i].out_a[cos][out] = prob;
	fprintf(stderr, "smodel_set_transition(0):State %d, %d, = %f\n", i, j, prob);
	break;
      }
    }

    for(in=0; in < smo->s[j].in_states; in++) {
      if ( smo->s[j].in_id[in] == i ) {
	smo->s[j].in_a[cos][in] = prob;
	break;
      }
    }
  }
}


double smodel_get_transition(smodel *smo, int i, int j, int cos) {
  int in, out;
  if (cos >= smo->cos) {
    fprintf(stderr, "smodel_get_transition(0): cos > state->cos\n");
    exit(-1);	
  }
  if (smo->s && smo->s[i].out_a && smo->s[j].in_a) {
    for(out=0; out < smo->s[i].out_states; out++) {
      if ( smo->s[i].out_id[out] == j ) {
	return 	smo->s[i].out_a[cos][out];
      }
    }
  }
  fprintf(stderr, "smodel_get_transition(1): data structure not initialized\n");
  return 0.0;
}

void smodel_set_mean(smodel *smo, int i, double *mu) {
  int m;
  if (smo->s != NULL) {
    for(m = 0; m < smo->M; m++)
      smo->s[i].mue[m] = mu[m];
  }
}

void smodel_set_variance(smodel *smo, int i, double *variance) {
  int m;
  if (smo->s != NULL) {
    for(m = 0; m < smo->M; m++) {
      smo->s[i].u[m] = variance[m];
      assert( smo->s[i].u[m] > 0.0 );
    }
  }
}


void call_smodel_print(char *filename, smodel *smo) {
  FILE *fp=fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);    
  } else {
    smodel_print(fp, smo);
    fclose(fp);
  }
}

int smodel_sorted_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps, int *seq_rank) {
  int matched, res;
  double log_p_i;
  int i;
  res=-1;
  matched=0;
  for ( i = 0; i < sqd->seq_number; i++) {
    seq_rank[i] = i;
    if (sfoba_logp(smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) { 
      log_ps[i] = log_p_i;
      matched++;
    }
    else  {
      /* Test: very small log score for sequence cannot be produced */
      log_ps[i] = -DBL_MAX;
      /*      mes(MES_WIN, "sequence[%d] can't be build.\n", i); */
    }
  }

  res=matched;
  if (matched == 0) { 
    fprintf(stderr, "smodel_likelihood: NO sequence can be build.\n"); 
  } else  qsort( seq_rank, sqd->seq_number, sizeof(int), (int (*)(const void *, const void *))seq_rank);

  /* return number of "matched" sequences */
  return res;
}


int smodel_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps) {
  int matched, res;
  double log_p_i;
  int i;
  res=-1;
  matched=0;
  for ( i = 0; i < sqd->seq_number; i++) {
    if (sfoba_logp(smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) { 
      log_ps[i] = log_p_i;
      matched++;
    }
    else  {
      /* Test: very small log score for sequence cannot be produced. */
      log_ps[i] = -DBL_MAX;
      /* fprintf(stderr,"sequence[%d] cannot be build.\n", i); */
    }
  }

  smodel_print( stderr, smo );

  res=matched;
  if (matched == 0) { 
    fprintf(stderr, "smodel_likelihood: NO sequence can be build.\n"); 
  }

  /* return number of "matched" sequences */
  return res;
}
