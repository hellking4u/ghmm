
/*******************************************************************************
  author       : Wasinee Rungsarityotin
  filename     : ghmm/ghmm/sdviterbi.c
  created      : TIME: 09:07:57     DATE: April, 2003
  $Id$

__copyright__

*******************************************************************************/

#include "mprintf.h"
#include "mes.h"
#include <float.h>
#include <math.h>
#include <assert.h>

#include <ghmm/ghmm.h>
#include <ghmm/matrix.h>
#include <ghmm/sdmodel.h>

typedef enum DFSFLAG
  { DONE, NOTVISITED, VISITED } DFSFLAG;


typedef struct local_store_t {
  double ***log_in_a;
  double **log_b;
  double *phi;
  double *phi_new;
  int    **psi;

  DFSFLAG *colors;
  int    *stateind;
  int    *doneTime;
  int    *topo_order;
} local_store_t;

static local_store_t *viterbi_alloc(sdmodel *mo, int len);
static int viterbi_free(local_store_t **v, int n, int cos, int len);

/*----------------------------------------------------------------------------*/
static local_store_t *viterbi_alloc(sdmodel *mo, int len) {
#define CUR_PROC "viterbi_alloc"
  local_store_t* v = NULL;
  int j;
  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}

  /* Allocate the log_in_a's -> individal lenghts */

  if ( v->log_in_a = (double***) malloc(sizeof(double**)*mo->N )) {
    for (j = 0; j < mo->N; j++)
      v->log_in_a[j] = matrix_d_alloc(mo->cos,  mo->s[j].in_states);
  } else {mes_proc(); goto STOP;}

  v->log_b = matrix_d_alloc(mo->N, len);
  if (!(v->log_b)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi, mo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi_new, mo->N)) {mes_proc(); goto STOP;}
  v->psi = matrix_i_alloc(len, mo->N);
  if (!(v->psi)) {mes_proc(); goto STOP;}
  return(v);
STOP:
  viterbi_free(&v, mo->N, mo->cos, len);
  return(NULL);
#undef CUR_PROC
} /* viterbi_alloc */


/*----------------------------------------------------------------------------*/
static int viterbi_free(local_store_t **v, int n, int cos, int len) {
#define CUR_PROC "viterbi_free"
  int j;
  mes_check_ptr(v, return(-1));
  if( !*v ) return(0);
  for (j = 0; j < n; j++)
    matrix_d_free(&((*v)->log_in_a[j]), cos);
  m_free((*v)->log_in_a);
  matrix_d_free( &((*v)->log_b), n );
  m_free((*v)->phi);
  m_free((*v)->phi_new);
  matrix_i_free( &((*v)->psi), len );
  m_free(*v);
  return(0);
#undef CUR_PROC
} /* viterbi_free */


static void Viterbi_precompute( sdmodel *mo, int *o, int len, local_store_t *v)
{
#define CUR_PROC "viterbi_precompute"
  int i, j, k, t, osc;
  double log_p = +1;
  
/* Precomputing the log(a_ij) */
  
  //for (j = 0; j < mo->N; j++)
  // for (i = 0; i < mo->s[j].in_states; i++)
  //  if ( mo->s[j].in_a[i] == 0.0 )   /* DBL_EPSILON ? */
  //log_in_a[j][i] = +1; /* Not used any further in the calculations */
  //  else
  //log_in_a[j][i] = log( mo->s[j].in_a[i] );

  /* Collecting emitting states */
  for(j = 0; j < mo->N; j++) 
    {
      // if (!mo->silent[j]) stateind.push_back(j);
    }
  
  for(j = 0; j < mo->N; j++) 
  {
    for(k = 0; k < mo->cos; k++)
      for(i = 0; i < mo->s[j].in_states; i++)	
	if ( mo->s[j].in_a[k][i] == 0.0 )   /* DBL_EPSILON ? */
	  v->log_in_a[j][k][i] = +1; /* Not used any further in the calculations */
	else 
	  v->log_in_a[j][k][i] = log( mo->s[j].in_a[k][i] );	
  }


  /* Precomputing the log(bj(ot)) */
  for (j = 0; j < mo->N; j++)
    for (t = 0; t < len; t++)
      if ( mo->s[j].b[o[t]] == 0.0 )   /* DBL_EPSILON ? */ 
	v->log_b[j][t] = +1; 
      else
	v->log_b[j][t] = log( mo->s[j].b[o[t]] );
#undef CUR_PROC
}/* viterbi_precompute */


/** Return the log score of the sequence */
int *sdviterbi( sdmodel *mo, int *o, int len, double *log_p)
{
#define CUR_PROC "sdviterbi"

  int *state_seq = NULL;
  int t, j, i;
  int last_osc = -1;
  double value, max_value;  
  double sum, osum = 0.0;
  double dummy = 0.0;
  local_store_t *v;

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  v = viterbi_alloc(mo, len);
  if (!v)                        {mes_proc(); goto STOP;}
  if (!m_calloc(state_seq, len)) { mes_proc(); goto STOP; }


  /* Precomputing the log(a_ij) and log(bj(ot)) */
  Viterbi_precompute(mo, o, len, v);

  /* Initialization, that is t = 0 */
  for (j = 0; j < mo->N; j++) {
    if ( mo->s[j].pi == 0.0 || v->log_b[j][0] == +1 ) /* instead of 0, DBL_EPS.? */
      v->phi[j] = +1;
    else
      v->phi[j] = log(mo->s[j].pi) + v->log_b[j][0];

    /*printf("phi[%d],%f\n", j, phi[j]);*/
  }
  /* psi[0][i] = 0, also unneccessary here */

  /* Recursion step */
  for (t = 1; t < len; t++) {

    int osc = mo->get_class(&dummy,t,&osum);

    for (j = 0; j < mo->N; j++) {
      /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      max_value = -DBL_MAX;
      v->psi[t][j] = -1;
      for (i = 0; i < mo->s[j].in_states; i++) {

	if ( last_osc != -1 && 
	     last_osc != osc ) { /* just switch class */
	  
	  if (   v->phi[ mo->s[j].in_id[i] ] != +1 ) {
	    if ( v->phi[ mo->s[j].in_id[i] ] > max_value ) {
	      max_value = v->phi[ mo->s[j].in_id[i] ];
	      v->psi[t][j] = mo->s[j].in_id[i];			       
	    }
	  }

	} else {
	  if ( v->phi[ mo->s[j].in_id[i] ] != +1 &&
	       v->log_in_a[j][osc][i]    != +1) {
	    value = v->phi[ mo->s[j].in_id[i] ] + v->log_in_a[j][osc][i];
	    if (value > max_value) {
	      max_value = value;
	      v->psi[t][j] = mo->s[j].in_id[i];
	    }
	  }
	  else
	    fprintf(stderr, " %d --> %d = %f, \n", i,j,v->log_in_a[j][osc][i]);
	}
      }

      last_osc = osc; /* save last transition class */

      /* No maximum found (that is, state never reached)
         or the output O[t] = 0.0: */
      if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	  v->log_b[j][t] == +1 ) {
	v->phi_new[j] = +1;
	fprintf(stderr, "b[%d][%d]=0.0, \n", j, t);
      }
      else
	v->phi_new[j] = max_value + v->log_b[j][t];
    } /* complete time step */

    /* First now replace the old phi with the new phi */
    for (j = 0; j < mo->N; j++) 
      {      
	v->phi[j] = v->phi_new[j];
	printf("\npsi[%d],%d, phi, %f\n", t, v->psi[t][j], v->phi[j]);
      }
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[len-1] = -1;
  for (j = 0; j < mo->N; j++)
    if ( v->phi[j] != +1 && v->phi[j] > max_value) { 
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

  return (state_seq);
STOP:
  /* Free the memory space */
  viterbi_free(&v, mo->N, mo->cos, len);
  m_free(state_seq);
  return NULL;
#undef CUR_PROC
} /* viterbi */

int *sdviterbi_silent(sdmodel *mo, int *o, int len, double *log_p)
{
#define CUR_PROC "sdviterbi_silent"


#undef CUR_PROC
}


int sdfoba_forward(sdmodel *mo, const int *O, int length, double **alpha, 
		   double *scale, double *log_p)
{
#define CUR_PROC "sdfoba_forward"


STOP:
  /* Free the memory space */
  return 0;
#undef CUR_PROC

}

