/*******************************************************************************
  author       : Bernhard Knab
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/sviterbi.c
  created      : TIME: 17:08:46     DATE: Mon 22. February 1999
  last-modified: TIME: 15:14:46     DATE: Wed 23. May 2001
*******************************************************************************/
/* $Id$ */

#include "config.h"

#if __EXPERIMENTAL__ == 3

#include <float.h>
#include <math.h>
#include "cviterbi.h"
#include "matrix.h"
#include "cmodel.h"

typedef struct local_store_t {
  double **log_in_a;
  double **log_b;
  double *phi;
  double *phi_new;
  int **psi;
} local_store_t;

static local_store_t *cviterbi_alloc(cmodel *cmo, int T);
static int cviterbi_free(local_store_t **v, int n, int T);

/*----------------------------------------------------------------------------*/
static local_store_t *cviterbi_alloc(cmodel *cmo, int T) {
#define CUR_PROC "cviterbi_alloc"
  local_store_t* v = NULL;
  int j;
  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}
  /* Allozieren der log_in_a's -> individuelle Laenge */
  if (!m_calloc(v->log_in_a, cmo->N)) {mes_proc(); goto STOP;}
  for (j = 0; j < cmo->N; j++)
    if (!m_calloc(v->log_in_a[j], cmo->s[j].in_states)) {mes_proc(); goto STOP;}
  v->log_b = matrix_d_alloc(cmo->N, T);
  if (!(v->log_b)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi, cmo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi_new, cmo->N)) {mes_proc(); goto STOP;}
  v->psi = matrix_i_alloc(T, cmo->N);
  if (!(v->psi)) {mes_proc(); goto STOP;}
  return(v);
STOP:
  cviterbi_free(&v, cmo->N, T);
  return(NULL);
#undef CUR_PROC
} /* cviterbi_alloc */


/*----------------------------------------------------------------------------*/
static int cviterbi_free(local_store_t **v, int n, int T) {
#define CUR_PROC "cviterbi_free"
  int j;
  mes_check_ptr(v, return(-1));
  if( !*v ) return(0);
  for (j = 0; j < n; j++)
    m_free((*v)->log_in_a[j]);
  m_free((*v)->log_in_a);
  matrix_d_free( &((*v)->log_b), n );
  m_free((*v)->phi);
  m_free((*v)->phi_new);
  matrix_i_free( &((*v)->psi), T );
  m_free(*v);
  return(0);
#undef CUR_PROC
} /* cviterbi_free */


/*----------------------------------------------------------------------------*/
static void cviterbi_precompute(cmodel *cmo, double *o, int T, 
				local_store_t *v) {
#define CUR_PROC "cviterbi_precompute"
  int i, j, t;
  double b;
  /* Precomputing der log(a_ij) */
  for (j = 0; j < cmo->N; j++)
    for (i = 0; i < cmo->s[j].in_states; i++)
      if ( cmo->s[j].in_a[i] == 0.0 )   /* DBL_EPSILON ? */
	v->log_in_a[j][i] = -DBL_MAX; 
      else
	v->log_in_a[j][i] = log( cmo->s[j].in_a[i] );
  /* Precomputing der log(bj(ot)) */
  for (j = 0; j < cmo->N; j++) {
    for (t = 0; t < T; t++) {
      b = cmodel_calc_b(cmo, j, o[t]);
      if ( b == 0.0 )   /* DBL_EPSILON ? */ 
	v->log_b[j][t] = -DBL_MAX; 
      else
	v->log_b[j][t] = log(b);
    }
  }
#undef CUR_PROC
} /* cviterbi_precompute */
  

/*============================================================================*/
int *cviterbi(cmodel *cmo, double *o, int T, double *log_p) {
#define CUR_PROC "cviterbi"
  int *state_seq = NULL;
  int t, j, i;
  double value, max_value;
  local_store_t *v;

  /* Allozieren der Matritzen log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  v = cviterbi_alloc(cmo, T);
  if (!v) {mes_proc(); goto STOP;}
  if (!m_calloc(state_seq, T)) { mes_proc(); goto STOP; }
  /* Precomputing der log(a_ij) und log(bj(ot)) */
  cviterbi_precompute(cmo, o, T, v);

  /* Initialisierung, d.h. t = 0 */
  for (j = 0; j < cmo->N; j++) {
    if ( cmo->s[j].pi == 0.0 || v->log_b[j][0] == -DBL_MAX ) /*0 oder DBL_EPS?*/
      v->phi[j] = -DBL_MAX;
    else
      v->phi[j] = log(cmo->s[j].pi) + v->log_b[j][0];
  }
  /* psi[0][i] = 0, also ueberfluessig hier */

  /* Rekursionsschritt */
  for (t = 1; t < T; t++) {
    for (j = 0; j < cmo->N; j++) {
      /* Maximum bestimmen */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      max_value = -DBL_MAX;
      v->psi[t][j] = -1;
      for (i = 0; i < cmo->s[j].in_states; i++) {
	if ( v->phi[ cmo->s[j].in_id[i] ] > -DBL_MAX && 
	     v->log_in_a[j][i] > -DBL_MAX) {
	  value = v->phi[ cmo->s[j].in_id[i] ] + v->log_in_a[j][i];
	  if (value > max_value) {
	    max_value = value;
	    v->psi[t][j] = cmo->s[j].in_id[i];
	  }
	}
      }
      /* Kein Max. gefunden (d.h. Zust. wird nie erreicht)
         oder Ausgabe O[t] = 0.0: */
      if (max_value == -DBL_MAX || /* und damit auch: (v->psi[t][j] == -1) */
	  v->log_b[j][t] == -DBL_MAX ) {
	v->phi_new[j] = -DBL_MAX;
      }
      else
	v->phi_new[j] = max_value + v->log_b[j][t];
    }

    /* jetzt erst altes phi ersetzen durch neues */
    for (j = 0; j < cmo->N; j++) 
      v->phi[j] = v->phi_new[j];
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[T-1] = -1;
  for (j = 0; j < cmo->N; j++)
    if (v->phi[j] != -DBL_MAX && v->phi[j] > max_value) { 
      max_value = v->phi[j];
      state_seq[T-1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* Sequenz kann vom Modell nicht gebildet werden! */
    *log_p = -DBL_MAX;
    /* Backtracking fkt. nicht, da state_seq[*] mit -1 belegt wurde */
    for (t = T - 2; t >= 0; t--)
      state_seq[t] = -1;    
  }
  else {
    *log_p = max_value;
    /* Backtracking */
    for (t = T - 2; t >= 0; t--)
      state_seq[t] = v->psi[t+1][state_seq[t+1]];
  }
  /* cviterbi_free(&v, cmo->N, T); ??? Hier auch notwendig? */
  return(state_seq);
STOP:
  /* Freiraeumen der Speicherplaetze ... */
  cviterbi_free(&v, cmo->N, T);
  m_free(state_seq);
  return NULL;
#undef CUR_PROC
} /* cviterbi */


/*============================================================================*/
/*  Kontrollrechnung: logP zu geg. Zustandspfad, Sequenz und Modell..*/
int cviterbi_logp(cmodel *cmo, double *o, int T, int *state_seq, double *log_p){
#define CUR_PROC "cviterbi_logp"
  int t, i, j, id;
  double b;
  *log_p = -DBL_MAX;
  /* t == 0 */
  i = state_seq[0];
  b = cmodel_calc_b(cmo, i, o[0]);
  if ( cmo->s[i].pi == 0.0 || b == 0.0 ) /* statt 0 DBL_EPS.? */
    return(-1);
  else
    *log_p = log(cmo->s[i].pi) + log(b);
  /* t == 1 .. T-1 */
  for (t = 1; t < T; t++) {
    i = state_seq[t-1];
    j = state_seq[t];
    /* a[i][j] suchen */
    for (id = 0; id < cmo->s[j].in_states; id++) {
      if (cmo->s[j].in_id[id] == i)
	break;
    }
    b = cmodel_calc_b(cmo, j, o[t]);
    if (id == cmo->s[j].in_states || cmo->s[j].in_a[id] == 0.0 || b == 0.0) {
      *log_p = -DBL_MAX;
      return(-1);
    }
    else 
      *log_p += log(cmo->s[j].in_a[id]) + log(b);
  }
  return(0);
#undef CUR_PROC
} /* cviterbi_logp */

#endif /* __EXPERIMENTAL__ == 3 */
