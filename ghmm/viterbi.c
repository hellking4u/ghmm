/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/viterbi.c
  created      : TIME: 09:07:57     DATE: Fri 19. December 1997
  $Id$

__copyright__

*******************************************************************************/


#include "mprintf.h"
#include "mes.h"
#include "viterbi.h"
#include "matrix.h"
#include "model.h"
#include <float.h>
#include <math.h>

/* matrices are allocated using stat_matrix_d_alloc */
typedef struct local_store_t {
  double **log_in_a;
  double **log_b;
  double *phi;
  double *phi_new;
  int **psi;
} local_store_t;

static local_store_t *viterbi_alloc(model *mo, int len);
static int viterbi_free(local_store_t **v, int n, int len);

/*----------------------------------------------------------------------------*/
static local_store_t *viterbi_alloc(model *mo, int len) {
#define CUR_PROC "viterbi_alloc"
  local_store_t* v = NULL;
  int j;
  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}

  /* Allocate the log_in_a's -> individal lenghts */
  if (!m_calloc(v->log_in_a, mo->N)) {mes_proc(); goto STOP;}
  for (j = 0; j < mo->N; j++)
    if (!m_calloc(v->log_in_a[j], mo->s[j].in_states)) {mes_proc(); goto STOP;}
  v->log_b = stat_matrix_d_alloc(mo->N, len);
  if (!(v->log_b)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi, mo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi_new, mo->N)) {mes_proc(); goto STOP;}
  v->psi = matrix_i_alloc(len, mo->N);
  if (!(v->psi)) {mes_proc(); goto STOP;}
  return(v);
STOP:
  viterbi_free(&v, mo->N, len);
  return(NULL);
#undef CUR_PROC
} /* viterbi_alloc */


/*----------------------------------------------------------------------------*/
static int viterbi_free(local_store_t **v, int n, int len) {
#define CUR_PROC "viterbi_free"
  int j;
  mes_check_ptr(v, return(-1));
  if( !*v ) return(0);
  for (j = 0; j < n; j++)
    m_free((*v)->log_in_a[j]);
  m_free((*v)->log_in_a);
  stat_matrix_d_free( &((*v)->log_b));
  m_free((*v)->phi);
  m_free((*v)->phi_new);
  matrix_i_free( &((*v)->psi), len );
  m_free(*v);
  return(0);
#undef CUR_PROC
} /* viterbi_free */


/*----------------------------------------------------------------------------*/
static void viterbi_precompute(model *mo, int *o, int len, local_store_t *v) {
#define CUR_PROC "viterbi_precompute"
  int i, j, t;
  /* Precomputing the log(a_ij) */
  for (j = 0; j < mo->N; j++)
    for (i = 0; i < mo->s[j].in_states; i++)
      if ( mo->s[j].in_a[i] == 0.0 )   /* DBL_EPSILON ? */
	v->log_in_a[j][i] = +1; /* Not used any further in the calculations */
      else
	v->log_in_a[j][i] = log( mo->s[j].in_a[i] );
  /* Precomputing the log(bj(ot)) */
  for (j = 0; j < mo->N; j++)
    for (t = 0; t < len; t++)
      if ( mo->s[j].b[o[t]] == 0.0 )   /* DBL_EPSILON ? */ 
	v->log_b[j][t] = +1; 
      else
	v->log_b[j][t] = log( mo->s[j].b[o[t]] );
#undef CUR_PROC
} /* viterbi_precompute */
  

/*============================================================================*/
int *viterbi(model *mo, int *o, int len, double *log_p) {
#define CUR_PROC "viterbi"
  int *state_seq = NULL;
  int t, j, i;
  double value, max_value;
  local_store_t *v;

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  v = viterbi_alloc(mo, len);
  if (!v) {mes_proc(); goto STOP;}
  if (!m_calloc(state_seq, len)) { mes_proc(); goto STOP; }
  /* Precomputing the log(a_ij) and log(bj(ot)) */
  viterbi_precompute(mo, o, len, v);

  /* Initialization, that is t = 0 */
  for (j = 0; j < mo->N; j++) {
    if ( mo->s[j].pi == 0.0 || v->log_b[j][0] == +1 ) /* instead of 0, DBL_EPS.? */
      v->phi[j] = +1;
    else
      v->phi[j] = log(mo->s[j].pi) + v->log_b[j][0];
  }
  /* psi[0][i] = 0, also unneccessary here */

  /* Recursion step */
  for (t = 1; t < len; t++) {
    for (j = 0; j < mo->N; j++) {
      /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      max_value = -DBL_MAX;
      v->psi[t][j] = -1;
      for (i = 0; i < mo->s[j].in_states; i++) {
	if ( v->phi[ mo->s[j].in_id[i] ] != +1 && v->log_in_a[j][i] != +1) {
	  value = v->phi[ mo->s[j].in_id[i] ] + v->log_in_a[j][i];
	  if (value > max_value) {
	    max_value = value;
	    v->psi[t][j] = mo->s[j].in_id[i];
	  }
	}
      }
      /* No maximum found (that is, state never reached)
         or the output O[t] = 0.0: */
      if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	  v->log_b[j][t] == +1 ) {
	v->phi_new[j] = +1;
      }
      else
	v->phi_new[j] = max_value + v->log_b[j][t];
    }

    /* First now replace the old phi with the new phi */
    for (j = 0; j < mo->N; j++) 
      v->phi[j] = v->phi_new[j];
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[len-1] = -1;
  for (j = 0; j < mo->N; j++)
    if (v->phi[j] != +1 && v->phi[j] > max_value) { 
      max_value = v->phi[j];
      state_seq[len-1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* Sequence can't be generated from the model! */
    *log_p = +1;
    /* Backtracing doesn't work, because state_seq[*] allocated with -1 */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = -1;    
  }
  else {
    *log_p = max_value;
    /* Backtracing */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = v->psi[t+1][state_seq[t+1]];
  }

  viterbi_free(&v, mo->N, len);
  return(state_seq);
STOP:
  /* Free the memory space */
  viterbi_free(&v, mo->N, len);
  m_free(state_seq);
  return NULL;
#undef CUR_PROC
} /* viterbi */


/*============================================================================*/
double viterbi_logp(model *mo, int *o, int len, int *state_seq) {
#define CUR_PROC "viterbi_logp"
  int t, i, j, id;
  double log_p = 0.0;

  /* t == 0 */
  i = state_seq[0];
  if ( mo->s[i].pi == 0.0 ||  mo->s[i].b[o[0]] == 0.0 ) /* 0 instead of DBL_EPS.? */
    return 0.0;
  else
    log_p = log(mo->s[i].pi) + log(mo->s[i].b[o[0]]);

  /* t == 1 .. len-1 */
  for (t = 1; t < len; t++) {
    i = state_seq[t-1];
    j = state_seq[t];
    /* search for a[i][j] */
    for (id = 0; id < mo->s[j].in_states; id++) {
      if (mo->s[j].in_id[id] == i)
	break;
    }
    if (id == mo->s[j].in_states || 
	mo->s[j].in_a[id] == 0.0 ||
	mo->s[j].b[o[t]] == 0.0)
      return 0.0;
    else 
      log_p += log(mo->s[j].in_a[id]) + log(mo->s[j].b[o[t]]);
  }

  return log_p;

#undef CUR_PROC
} /* viterbi_logp */
