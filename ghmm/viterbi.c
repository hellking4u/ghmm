/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/viterbi.c
*       Authors:  Wasinee Rungsarityotin, Benjamin Georgi
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
*       This file is version $Revision$
*                       from $Date$
*             last change by $Author$.
*
*******************************************************************************/


#include "mprintf.h"
#include "mes.h"
#include <float.h>
#include <math.h>
#include <assert.h>

#include <ghmm/ghmm.h>
#include <ghmm/matrix.h>
#include <ghmm/sdmodel.h>

typedef struct local_store_t {
  double **log_in_a;
  double **log_b;
  double *phi;
  double *phi_new;
  int    **psi;

  int    *topo_order;
  int    topo_order_length;
} local_store_t;

static local_store_t *viterbi_alloc(model *mo, int len);
static int viterbi_free(local_store_t **v, int n,  int len);

/*----------------------------------------------------------------------------*/
static local_store_t *viterbi_alloc(model *mo, int len) {
#define CUR_PROC "sdviterbi_alloc"
  local_store_t* v = NULL;
  int j;
  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}

  /* Allocate the log_in_a's -> individal lenghts */

  if (!m_calloc(v->log_in_a, mo->N)) {mes_proc(); goto STOP;}
  for (j = 0; j < mo->N; j++){
    if (!m_calloc(v->log_in_a[j], mo->s[j].in_states)) {mes_proc(); goto STOP;}
  }

  v->log_b = matrix_d_alloc(mo->N, len);
  if (!(v->log_b)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi, mo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi_new, mo->N)) {mes_proc(); goto STOP;}
  v->psi = stat_matrix_i_alloc(len, mo->N);
  if (!(v->psi)) {mes_proc(); goto STOP;}

  v->topo_order_length = 0;
  if (!m_calloc(v->topo_order, mo->N)) {mes_proc(); goto STOP;}

  return(v);
STOP:
  viterbi_free(&v, mo->N,  len);
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
  matrix_d_free( &((*v)->log_b), n );
  m_free((*v)->phi);
  m_free((*v)->phi_new);
  //matrix_i_free( &((*v)->psi), len );
  stat_matrix_i_free( &((*v)->psi));
  m_free((*v)->topo_order);
  m_free(*v);
  return(0);
#undef CUR_PROC
} /* viterbi_free */


static void Viterbi_precompute( model *mo, int *o, int len, local_store_t *v)
{
#define CUR_PROC "viterbi_precompute"
  int i, j, k, t;
  double log_p = +1;
  
  /* Precomputing the log(a_ij) */
  
  for (j = 0; j < mo->N; j++)
   for (i = 0; i < mo->s[j].in_states; i++)
    if ( mo->s[j].in_a[i] == 0.0 )   /* DBL_EPSILON ? */
 	  v->log_in_a[j][i] = +1; /* Not used any further in the calculations */
    else
  	  v->log_in_a[j][i] = log( mo->s[j].in_a[i] );


  /* Precomputing the log(bj(ot)) */
  for (j = 0; j < mo->N; j++)
    for (t = 0; t < len; t++) {
      if (o[t] != mo->M){
	if ( mo->s[j].b[o[t]] == 0.0 )   /* DBL_EPSILON ? */ 
	  v->log_b[j][t] = +1; 
	else
	  v->log_b[j][t] = log( mo->s[j].b[o[t]] );
      }
      else{
	v->log_b[j][t] = 0.0; 
      }
    }

#undef CUR_PROC
}/* viterbi_precompute */

/** */
static void __viterbi_silent( model *mo, int t, local_store_t *v )
{
  int topocount;
  int i, k;
  double max_value, value;

  for ( topocount = 0; topocount < mo->topo_order_length; topocount++) {
      k = mo->topo_order[topocount];
      if ( mo->silent[k] ) { /* Silent states */	
	  /* Determine the maximum */
	  /* max_phi = phi[i] + log_in_a[j][i] ... */
	  max_value = -DBL_MAX;
	  v->psi[t][k] = -1;
	  for (i = 0; i < mo->s[k].in_states; i++) {
	    
	  if ( v->phi[ mo->s[k].in_id[i] ] != +1 &&
		 v->log_in_a[k][i]      != +1) {
	      value = v->phi[ mo->s[k].in_id[i] ] + v->log_in_a[k][i];
	      if (value > max_value) {
		max_value = value;
		v->psi[t][k] = mo->s[k].in_id[i];
	      }
	    }
	  }
	  
	  /* No maximum found (that is, state never reached)
	     or the output O[t] = 0.0: */
	  if (max_value    == -DBL_MAX) {
	    v->phi[k] = +1;
	  } else {
	    v->phi[k] = max_value;
	  }
	  
	}
  }
}

/** Return the log score of the sequence */
int *viterbi( model *mo, int *o, int len, double *log_p)
{
#define CUR_PROC "viterbi"

  int *state_seq = NULL;
  int t, j, i, k, St;
  int topocount = 0;
  double value, max_value;  
  double sum, osum = 0.0;
  double dummy = 0.0;
  local_store_t *v;
  int len_path  = mo->N*len;
  int lastemState;
  
  // printf("---- viterbi -----\n");

  if (mo->model_type == kSilentStates && 
      mo->silent != NULL && 
      mo->topo_order == NULL) {
    model_topo_ordering( mo ); /* Should we call it here ???? */
  }

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  v = viterbi_alloc(mo, len);
  if (!v)                        { mes_proc(); goto STOP; }
  if (!m_calloc(state_seq, len_path)) { mes_proc(); goto STOP; }
  for(i=0; i < len_path ; i++) {
    state_seq[i] = -1;
  }
  
  /* Precomputing the log(a_ij) and log(bj(ot)) */
  Viterbi_precompute(mo, o, len, v);

  /* Initialization, that is t = 0 */
  for (j = 0; j < mo->N; j++) {
    if ( mo->s[j].pi == 0.0 || v->log_b[j][0] == +1 ) /* instead of 0, DBL_EPS.? */
      v->phi[j] = +1;
    else
      v->phi[j] = log(mo->s[j].pi) + v->log_b[j][0];
  }
  if ( mo->model_type == kSilentStates ) { /* could go into silent state at t=0 */
    int osc;
    __viterbi_silent( mo, t=0, v);
  }
  /*for (j = 0; j < mo->N; j++)
    {
      printf("\npsi[%d],in:%d, phi=%f\n", t, v->psi[t][j], v->phi[j]);
    }

  for( i = 0; i < mo->N; i++){
    printf("%d\t", former_matchcount[i]);
  }
  for (i = 0; i < mo->N; i++){
    printf("%d\t", recent_matchcount[i]);
  }*/
  
  /* t > 0 */
  for (t = 1; t < len; t++) {

    //int osc = mo->get_class(mo->N,t);

    for (j = 0; j < mo->N; j++) /** initialization of phi, psi **/
    {
      v->phi_new[j] = +1;
      v->psi[t][j]  = -1;
    }

    for (k = 0; k < mo->N; k++) {

	  /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      if ( mo->model_type != kSilentStates || !mo->silent[k] ) {
		St        = k;
		max_value = -DBL_MAX;
		v->psi[t][St] = -1;
		for (i = 0; i < mo->s[St].in_states; i++) {

		  if ( v->phi[ mo->s[St].in_id[i] ] != +1 && v->log_in_a[St][i] != +1) {
		    value = v->phi[ mo->s[St].in_id[i] ] + v->log_in_a[St][i];
	   	    if (value > max_value) {
	      	  max_value = value;
	      	  v->psi[t][St] = mo->s[St].in_id[i];
	        }
	    }
	    else
	      {;} // fprintf(stderr, " %d --> %d = %f, \n", i,St,v->log_in_a[St][i]);
	  }

	/* No maximum found (that is, state never reached)
	   or the output O[t] = 0.0: */
	if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	    v->log_b[St][t] == +1 ) {
	    v->phi_new[St] = +1;
	}
	else
	  v->phi_new[St] = max_value + v->log_b[St][t];

      }
    } /* complete time step for emitting states */

    /* First now replace the old phi with the new phi */
    for (j = 0; j < mo->N; j++) {      
		v->phi[j] = v->phi_new[j];
		//printf("\npsi[%d],%d, phi, %f\n", t, v->psi[t][j], v->phi[j]); 
    }
    
	//last_osc = osc; /* save last transition class */

    if ( mo->model_type == kSilentStates ) { 
      __viterbi_silent( mo, t, v );
    } /* complete time step for silent states */
    
    /**************
    for (j = 0; j < mo->N; j++) 
      {      
	printf("\npsi[%d],in:%d, phi=%f\n", t, v->psi[t][j], v->phi[j]);
       }
      
    for (i = 0; i < mo->N; i++){
      printf("%d\t", former_matchcount[i]);
    }

    for (i = 0; i < mo->N; i++){
      printf("%d\t", recent_matchcount[i]);
    }
    ****************/

  } /* Next observation , increment time-step */

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[len_path-1] = -1;
  for (j = 0; j < mo->N; j++)
    if ( v->phi[j] != +1 && v->phi[j] > max_value) { 
      max_value = v->phi[j];
      state_seq[len_path-1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* Sequence can't be generated from the model! */
    *log_p = +1;
    /* Backtracing doesn't work, because state_seq[*] allocated with -1 */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = -1;    
  }
  else {
    /* Backtracing, should put DEL path nicely */
    *log_p = max_value;
    lastemState = state_seq[len_path-1];
    for(t = len - 2, i=len_path-2; t >= 0; t--) {
      if ( mo->model_type == kSilentStates &&
	   mo->silent[ v->psi[t+1][ lastemState ] ]) {

	St = v->psi[t+1][ lastemState ];
	/* fprintf(stderr, "t=%d:  DEL St=%d\n", t+1, St ); */
	while( St != -1  && mo->silent[ St ] ) { /* trace back up to the last emitting state */
	  /* fprintf(stderr, "t=%d:  DEL St=%d\n", t, St ); */
	  state_seq[i--] = St;
	  St = v->psi[t][ St ];
	}
	state_seq[i--] = St;
	lastemState = St;
      } else {
	state_seq[i--] = v->psi[t+1][ lastemState ];
	lastemState    = v->psi[t+1][ lastemState ];
      }
    }
     
  }

  /* PRINT PATH */
  /* 
  fprintf(stderr, "Viterbi path: " );
  for(t=0; t < len_path; t++)
    if (state_seq[t] >= 0) fprintf(stderr, " %d ",  state_seq[t]);
  fprintf(stderr, "\n Freeing ... \n"); */ 

  /* Free the memory space */
  viterbi_free(&v, mo->N, len);
  return (state_seq);
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
  int* vpath;
  
  vpath = viterbi(mo,o,len, &log_p);
  

  return log_p;

#undef CUR_PROC
} /* viterbi_logp */

/*============================================================================*/
