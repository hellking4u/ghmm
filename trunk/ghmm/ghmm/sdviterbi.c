
/*******************************************************************************
  author       : Wasinee Rungsarityotin, 
  an extension from viterbi() in model.c with ordering of silent states
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

  int    *topo_order;
  int    topo_order_length;
} local_store_t;

static local_store_t *sdviterbi_alloc(sdmodel *mo, int len);
static int sdviterbi_free(local_store_t **v, int n, int cos, int len);

/*----------------------------------------------------------------------------*/
static local_store_t *sdviterbi_alloc(sdmodel *mo, int len) {
#define CUR_PROC "sdviterbi_alloc"
  local_store_t* v = NULL;
  int j;
  if (!m_calloc(v, 1)) {mes_proc(); goto STOP;}

  /* Allocate the log_in_a's -> individal lenghts */

  if ( v->log_in_a = (double***) malloc(sizeof(double**)*mo->N )) {
    for (j = 0; j < mo->N; j++)
      v->log_in_a[j] = stat_matrix_d_alloc(mo->cos,  mo->s[j].in_states);
  } else {mes_proc(); goto STOP;}

  v->log_b = stat_matrix_d_alloc(mo->N, len);
  if (!(v->log_b)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi, mo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(v->phi_new, mo->N)) {mes_proc(); goto STOP;}
  v->psi = matrix_i_alloc(len, mo->N);
  if (!(v->psi)) {mes_proc(); goto STOP;}

  v->topo_order_length = 0;
  if (!m_calloc(v->topo_order, mo->N)) {mes_proc(); goto STOP;}

  return(v);
STOP:
  sdviterbi_free(&v, mo->N, mo->cos, len);
  return(NULL);
#undef CUR_PROC
} /* viterbi_alloc */


/*----------------------------------------------------------------------------*/
static int sdviterbi_free(local_store_t **v, int n, int cos, int len) {
#define CUR_PROC "sdviterbi_free"
  int j;
  mes_check_ptr(v, return(-1));
  if( !*v ) return(0);
  for (j = 0; j < n; j++)
    stat_matrix_d_free(&((*v)->log_in_a[j]));
  m_free((*v)->log_in_a);
  stat_matrix_d_free( &((*v)->log_b));
  m_free((*v)->phi);
  m_free((*v)->phi_new);
  matrix_i_free( &((*v)->psi), len );
  m_free((*v)->topo_order);
  m_free(*v);
  return(0);
#undef CUR_PROC
} /* viterbi_free */

<<<<<<< sdviterbi.c
=======
/*----------------------------------------------------------------------------*/
static void __VisitNext(sdmodel *mo, int j, int *counter, local_store_t *v)
{
  int i, nextState, ins;

  v->colors[j] = VISITED; 
  v->topo_order[(*counter)++] = j;

  for(i = 0; i < mo->s[j].out_states; i++)	
    {
      if (v->colors[mo->s[j].out_id[i]] == NOTVISITED &&
	  mo->silent[ mo->s[j].out_id[i] ]) { /* looping back taken care of */
	nextState = mo->s[j].out_id[i];
	/* Check if all in-coming silent states has been visited */
	for( ins=0; ins < mo->s[nextState].in_states; ins++) {	
	  if ( nextState != mo->s[nextState].in_id[ins] &&
	       mo->silent[ mo->s[nextState].in_id[ins] ] ) {
	    if ( v->colors[ mo->s[nextState].in_id[ins] ] == NOTVISITED ) {
	      fprintf(stderr, "%d, %d to %d\n",j, ins, nextState);
	      goto find_next_silent;
	    }
	  }
	}
	
	v->colors[nextState] = VISITED; 
	v->topo_order[(*counter)++] = nextState; /* All in-coming silent states
						  * has been visited,
						  * and so we have the ordering
						  */

      }
 find_next_silent:;
    }
}



/*----------------------------------------------------------------------------*/
static void __sdmodel_topo_ordering(sdmodel *mo, local_store_t *v) 
{
  int i,j,k;

  assert(mo->model_type == kSilentStates); /* otherwise, why are you here? */

  v->colors   =   (DFSFLAG*)malloc( sizeof(DFSFLAG)*mo->N );
  v->topo_order = (int*)malloc( sizeof(int)*mo->N );
  v->topo_order_length = 0;

  for(i=0; i < mo->N; i++)
    {
      v->colors[i] = NOTVISITED;
    }
    
  for(i=0; i < mo->N; i++) {
    if ( mo->silent[i] && 
	 v->colors[i] == NOTVISITED ) {
      for(j = 0; j < mo->s[i].in_states; j++) {
	if ( i != mo->s[i].in_id[j] &&
	     mo->silent[ mo->s[i].in_id[j] ] ) {
	  goto find_start_points;
	}
      }
      fprintf(stderr, "Starting at %d\n", i);
      v->colors[i] = VISITED;
      v->topo_order[v->topo_order_length++] = i;
      for(j = 0; j < mo->s[i].out_states; j++)	
	{
	  if ( mo->silent[ mo->s[i].out_id[j] ]) {  /* HACK!!! for a starting point */
	    v->colors[ mo->s[i].out_id[j] ] = VISITED;
	    v->topo_order[v->topo_order_length++] = mo->s[i].out_id[j];
	  } 
	}
    }
find_start_points:;
  }

  for(i=0; i < mo->N; i++) {
    if ( mo->silent[i] && 
	 v->colors[i] == NOTVISITED ) {
      for(j = 0; j < mo->s[i].in_states; j++) {
	if ( i != mo->s[i].in_id[j] &&
	     mo->silent[ mo->s[i].in_id[j] ] ) {
	  /*
	   * If an in-coming transition is from a silent state, 
	   * it must be visited before.
	   */
	  if ( v->colors[ mo->s[i].in_id[j] ] == NOTVISITED ) {
	    /* fprintf(stderr, "%d to %d\n", mo->s[i].in_id[j], i); */
	    goto find_start; 
	  }
	}
      }
      __VisitNext(mo, i, &v->topo_order_length, v);	
    }
find_start:;
  }

}


/*----------------------------------------------------------------------------
void sdmodel_topo_ordering(sdmodel *mo) 
{
#define CUR_PROC "sdmodel_topo_ordering"
  int i;
  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi 
  local_store_t *v;

  v = viterbi_alloc(mo, 1);
  if (!v) { mes_proc(); goto STOP; }

  __sdmodel_topo_ordering( mo, v);
  
  mo->topo_order_length = v->topo_order_length;
  if (!m_calloc(mo->topo_order, mo->topo_order_length)) {mes_proc(); goto STOP;}

  for(i=0; i < v->topo_order_length; i++) {
    mo->topo_order[i] = v->topo_order[i];
  }
  fprintf(stderr,"Ordering silent states....\n\t");
  for(i=0; i < mo->topo_order_length; i++) {
    fprintf(stderr, "%d, ", mo->topo_order[i]);
  }
  /* viterbi_free(&v, mo->N, mo->cos, 1); Memory problem !!!! 
 STOP:
#undef CUR_PROC
} */

>>>>>>> 1.3

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
    for (t = 0; t < len; t++) {
      if ( mo->s[j].b[o[t]] == 0.0 )   /* DBL_EPSILON ? */ 
	v->log_b[j][t] = +1; 
      else
	v->log_b[j][t] = log( mo->s[j].b[o[t]] );
    }

#undef CUR_PROC
}/* viterbi_precompute */

/** */
static void __viterbi_silent( sdmodel *mo, int t, local_store_t *v, int *recent_matchcount, int *countstates, int nr_of_countstates )
{
  int topocount;
  int i, k, osc;
  double max_value, value;

  osc = 0;
  for ( topocount = 0; topocount < mo->topo_order_length; topocount++)
    {
      k = mo->topo_order[topocount];
      if ( mo->silent[k] ) /* Silent states */ 
	{
	  /* Determine the maximum */
	  /* max_phi = phi[i] + log_in_a[j][i] ... */
	  max_value = -DBL_MAX;
	  v->psi[t][k] = -1;
	  for (i = 0; i < mo->s[k].in_states; i++) {
	    // printf("\nBerrechnung von transclass von Zustand %d", mo->s[k].in_id[i]);
	    if (mo->cos != 1){
	      osc = mo->get_class(mo->N, recent_matchcount[mo->s[k].in_id[i]]);
	    }
	    if ( v->phi[ mo->s[k].in_id[i] ] != +1 &&
		 v->log_in_a[k][osc][i]      != +1) {
	      value = v->phi[ mo->s[k].in_id[i] ] + v->log_in_a[k][osc][i];
	      if (value > max_value) {
		max_value = value;
		v->psi[t][k] = mo->s[k].in_id[i];
	      }
	    }
	  }
	  /*find out, if we are in a delete state unless this state isn't reached anyway*/
	  if (v->psi[t][k] != -1){
	  for (i = 0; i < nr_of_countstates; i++){
	    if (k == countstates[i]){
	      recent_matchcount[k] = 1;
	      break;
	    }
	  }
	  recent_matchcount[k] += recent_matchcount[v->psi[t][k]];

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
}

/** Return the log score of the sequence */
int *sdviterbi( sdmodel *mo, int *o, int len, double *log_p)
{
#define CUR_PROC "sdviterbi_silent"

  int *state_seq = NULL;
  int t, j, i, k, St, osc;
  int topocount = 0;
  int last_osc = -1;
  double value, max_value;  
  double sum, osum = 0.0;
  double dummy = 0.0;
  local_store_t *v;
  int len_path  = mo->N*len;
  /*lists to remember how long we have been staying in the circular part of the model*/
  int *former_matchcount = NULL;
  int *recent_matchcount = NULL;
  int *tmp_matchcount = NULL;
  int *countstates = NULL;
  int nr_of_countstates = 2*((mo->N - 5)/3);	// # of matchstates + deletestates

  osc = 0;

  if (mo->model_type == kSilentStates && mo->topo_order == NULL) {
    /* sdmodel_topo_ordering( mo ); */
    /* Should we call it here ???? */
    fprintf(stderr, "Viterbi Error: Contain silent states, but no topological ordering\n");
    goto STOP;
  }

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  v = sdviterbi_alloc(mo, len);
  if (!v)                        { mes_proc(); goto STOP; }
  if (!m_calloc(state_seq, len_path)) { mes_proc(); goto STOP; }
  for(i=0; i < len_path ; i++) {
    state_seq[i] = -1;
  }
  
  if (!m_calloc(former_matchcount, mo->N)) { mes_proc(); goto STOP; }
  if (!m_calloc(recent_matchcount, mo->N)) { mes_proc(); goto STOP; }
  /*We always start outside of the circle, no way to get in in t=0 => matchcounts = 0*/
  for(i=0; i < mo->N; i++){
    former_matchcount[i] = 0;
    recent_matchcount[i] = 0;
  }
  
  if (!m_calloc(countstates, nr_of_countstates)) { mes_proc(); goto STOP; }
  for(i = 0; i < nr_of_countstates/2; i++){
    /* 5th state is first matchstate, then come all other matchstates and afterwards all deletestates
    so we have a list of states that inkrement the counter of how long we have been in the circle
    */
    countstates[i] = i+5;
    countstates[i+nr_of_countstates/2] = i+nr_of_countstates+4;
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
    __viterbi_silent( mo, t=0, v, recent_matchcount, countstates, nr_of_countstates );
  }
  for (j = 0; j < mo->N; j++) 
    {
      printf("\npsi[%d],in:%d, phi=%f\n", t, v->psi[t][j], v->phi[j]);
    }

  for( i = 0; i < mo->N; i++){
    printf("%d\t", former_matchcount[i]);
  }
  for (i = 0; i < mo->N; i++){
    printf("%d\t", recent_matchcount[i]);
  }
  /* t > 0 */
  for (t = 1; t < len; t++) {

<<<<<<< sdviterbi.c
    //int osc = mo->get_class(mo->N,t);
=======
    int osc = mo->get_class(*o, t);
>>>>>>> 1.3

    for (j = 0; j < mo->N; j++) /** initialization of phi, psi **/
    {
      v->phi_new[j] = +1;
      v->psi[t][j]  = -1;
    }

    for (k = 0; k < mo->N; k++) {

      recent_matchcount[k] = 0;
      /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      if ( mo->model_type != kSilentStates || !mo->silent[k] ) {
	St        = k;
	max_value = -DBL_MAX;
	v->psi[t][St] = -1;
	for (i = 0; i < mo->s[St].in_states; i++) {
	  // get_class of in state
	  // printf("\nBerechnung von transclass fuer Zustand %d", mo->s[St].in_id[i]);
	  if (mo->cos != 1){
	    osc = mo->get_class(mo->N, former_matchcount[mo->s[St].in_id[i]]);
	  }
	  if ( v->phi[ mo->s[St].in_id[i] ] != +1 &&
	       v->log_in_a[St][osc][i]    != +1) {
	    value = v->phi[ mo->s[St].in_id[i] ] + v->log_in_a[St][osc][i];
	    if (value > max_value) {
	      max_value = value;
	      v->psi[t][St] = mo->s[St].in_id[i];
	    }
	  }
	  else
	    fprintf(stderr, " %d --> %d = %f, \n", i,St,v->log_in_a[St][osc][i]);
	}

	/* No maximum found (that is, state never reached)
	   or the output O[t] = 0.0: */
	if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	    v->log_b[St][t] == +1 ) {
	    v->phi_new[St] = +1;
	}
	else
	  v->phi_new[St] = max_value + v->log_b[St][t];

	/*find out how long we have been staying in the circle unless we didn't reach this state*/
	if (v->psi[t][St] != -1){
	for(i = 0; i < nr_of_countstates; i++){
	  if (countstates[i] == St){
	    recent_matchcount[St] = 1;
	    break;
	  }
	}
	recent_matchcount[St] += former_matchcount[v->psi[t][St]];
	}

      }
    } /* complete time step for emitting states */

    /* First now replace the old phi with the new phi */
    for (j = 0; j < mo->N; j++) 
      {      
	v->phi[j] = v->phi_new[j];
	/* printf("\npsi[%d],%d, phi, %f\n", t, v->psi[t][j], v->phi[j]); */
      }
    last_osc = osc; /* save last transition class */

    if ( mo->model_type == kSilentStates ) { 
      __viterbi_silent( mo, t, v, recent_matchcount, countstates, nr_of_countstates );
    } /* complete time step for silent states */
    
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
    // swap the lists of matchcounts
    tmp_matchcount = former_matchcount;
    former_matchcount = recent_matchcount;
    recent_matchcount = tmp_matchcount;

  } /* Next observation , increment time-step */

  /* Termination */
  int lastemState;
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
  fprintf(stderr, "Viterbi path: " );
  for(t=0; t < len_path; t++)
    if (state_seq[t] >= 0) fprintf(stderr, " %d ",  state_seq[t]);
  fprintf(stderr, "\n");

  return (state_seq);
STOP:
  /* Free the memory space */
  sdviterbi_free(&v, mo->N, mo->cos, len);
  m_free(state_seq);
  m_free(former_matchcount);
  m_free(recent_matchcount);
  m_free(tmp_matchcount);
  m_free(countstates);
  return NULL;
#undef CUR_PROC
} /* viterbi */


int *sdviterbi_silent(sdmodel *mo, int *o, int len, double *log_p)
{
#define CUR_PROC "sdviterbi_silent"


#undef CUR_PROC
}



