/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/model.c
  created      : TIME: 10:47:27     DATE: Fri 19. December 1997
  $Id$

__copyright__

*******************************************************************************/


#include <float.h>
#include <math.h>
#include "model.h"
#include "matrix.h"
#include "sequence.h"
#include "const.h"
#include "rng.h"
#include "foba.h"
#include "mes.h"
#include "mprintf.h"
#include "string.h"
#include "ghmm.h"
#include "modelutil.h"

#define  __EPS 10e-6

typedef enum DFSFLAG
  { DONE, NOTVISITED, VISITED } DFSFLAG;


/*typedef struct local_store_t {
  DFSFLAG *colors;
  int    *topo_order;
  int    topo_order_length;
} local_store_t;

static local_store_t *topo_alloc(model *mo, int len);
static int topo_free(local_store_t **v, int n, int cos, int len); */


/*----------------------------------------------------------------------------*/
static int model_state_alloc(state *state, int M, int in_states,
			     int out_states) {
# define CUR_PROC "model_state_alloc"
  int res = -1;
  if(!m_calloc(state->b, M)) {mes_proc(); goto STOP;}
  if (out_states > 0) {
    if(!m_calloc(state->out_id, out_states)) {mes_proc(); goto STOP;}
    if(!m_calloc(state->out_a, out_states)) {mes_proc(); goto STOP;}
  }
  if (in_states > 0) {
    if(!m_calloc(state->in_id, in_states)) {mes_proc(); goto STOP;}
    if(!m_calloc(state->in_a, in_states)) {mes_proc(); goto STOP;}
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* model_state_alloc */

/*----------------------------------------------------------------------------*/

static int model_copy_vectors(model *mo, int index, double **a_matrix, 
			      double **b_matrix, double *pi, int *fix) {
#define CUR_PROC "model_copy_vectors"
  int i, cnt_out = 0, cnt_in = 0;
  mo->s[index].pi = pi[index];
  mo->s[index].fix = fix[index];
  for (i = 0; i < mo->M; i++)
    mo->s[index].b[i] = b_matrix[index][i];
  for (i = 0; i < mo->N; i++) {
    if (a_matrix[index][i]) { /* Transitions to a following state possible */
      if (cnt_out >= mo->s[index].out_states) {mes_proc(); return(-1);}
      mo->s[index].out_id[cnt_out] = i;
      mo->s[index].out_a[cnt_out] = a_matrix[index][i];
      cnt_out++;
    }
    if (a_matrix[i][index]) { /* Transitions to a previous state possible */
      if (cnt_in >= mo->s[index].in_states) {mes_proc(); return(-1);}
      mo->s[index].in_id[cnt_in] = i;
      mo->s[index].in_a[cnt_in] = a_matrix[i][index];
      cnt_in++;
    }
  }
  return(0);
#undef CUR_PROC
} /* model_copy_vectors */



/*============================================================================*/

/* Old prototype:

model **model_read(char *filename, int *mo_number, int **seq,
			 const int *seq_len, int seq_number) { */

model **model_read(char *filename, int *mo_number) {
#define CUR_PROC "model_read"
  int j;
  long new_models = 0;
  scanner_t *s = NULL;
  model **mo = NULL;
  *mo_number = 0;
  s = scanner_alloc(filename);  if(!s) {mes_proc(); goto STOP;}
  while(!s->err && !s->eof) {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP; 
    if (!strcmp(s->id, "HMM") || !strcmp(s->id, "hmm")) {
      (*mo_number)++;
      /* more mem */	 
      if (m_realloc(mo, *mo_number)) { mes_proc(); goto STOP; }
      mo[*mo_number - 1] = model_direct_read(s, (int *) &new_models); 
      if (!mo[*mo_number - 1]) { mes_proc(); goto STOP; }
      /* Copies the model, that has already been read. */
      if (new_models > 1) { 
	/* "-1" because mo_number++ has already been done. */
	if (m_realloc(mo, *mo_number - 1 + new_models)) { 
	  mes_proc(); goto STOP; 
	}
	for (j = 1; j < new_models; j++) {
	  mo[*mo_number] = model_copy(mo[*mo_number - 1]);
	  if (!mo[*mo_number]) { mes_proc(); goto STOP; }
	  (*mo_number)++;
	}
      }
    }  
    else if (!strcmp(s->id, "HMM_SEQ")) {    
      model **tmp_mo = NULL;
      tmp_mo = model_from_sequence_ascii(s, &new_models);
      if (m_realloc(mo, (*mo_number + new_models))) { mes_proc(); goto STOP; }
      for (j = 0; j < new_models; j++) {
	if (!tmp_mo[j]) { mes_proc(); goto STOP; }
	mo[*mo_number] = tmp_mo[j];
	(*mo_number)++;
      }	 
    }	 
    else {
      scanner_error(s, "unknown identifier");
      goto STOP;
    }
    scanner_consume(s, ';'); if(s->err) goto STOP; 
  } /* while(!s->err && !s->eof) */
  return mo;
STOP:
  return NULL;
#undef CUR_PROC
} /* model_read */



/*============================================================================*/

model *model_direct_read(scanner_t *s, int *multip){
#define CUR_PROC "model_direct_read"
  int i, m_read, n_read, a_read, b_read, pi_read, len, prior_read, fix_read;
  model *mo      = NULL;
  double **a_matrix = NULL, **b_matrix = NULL;
  double *pi_vector = NULL;
  int *fix_vector = NULL;
  m_read = n_read = a_read = b_read = pi_read = prior_read = fix_read = 0;
  *multip = 1; /* default */
  if (!(m_calloc(mo, 1))) { mes_proc(); goto STOP; }
  scanner_consume( s, '{' ); if(s->err) goto STOP; 
  while(!s->err && !s->eof && s->c - '}') {    
    scanner_get_name(s);
    if (strcmp(s->id, "M") && strcmp(s->id, "N") && strcmp(s->id, "Pi")
	&& strcmp(s->id, "A") && strcmp(s->id, "B") && 
	strcmp(s->id, "multip") && strcmp(s->id, "prior") &&
	strcmp(s->id, "fix_state")) {
      scanner_error(s, "unknown identifier"); goto STOP;
    }    
    scanner_consume(s, '='); if(s->err) goto STOP;
    if (!strcmp(s->id, "multip")) {
      *multip = scanner_get_int(s);
      if (*multip < 1) {/* Doesn't make any sense! */
	*multip = 1;
	mes_prot("Multip < 1 ignored\n");
      }
    }
    else if (!strcmp(s->id, "M")) {/*Number of output values*/
      if (m_read) {scanner_error(s, "identifier M twice"); goto STOP;}
      mo->M = scanner_get_int(s);
      m_read = 1;
    }
    else if (!strcmp(s->id, "N")) {/*Number of states*/
      if (n_read) {scanner_error(s, "identifier N twice"); goto STOP;}
      mo->N = scanner_get_int(s);
      if (!m_calloc(mo->s, mo->N)) {mes_proc(); goto STOP;}
      n_read = 1;
    }   
    else if (!strcmp(s->id, "A")) {/*Transition probability*/
      if (!n_read) {scanner_error(s, "need N as a range for A"); goto STOP;}	    
      if (a_read) {scanner_error(s, "identifier A twice"); goto STOP;}
      if (!m_calloc(a_matrix, mo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if (matrix_d_read(s, a_matrix, mo->N, mo->N)){
	  scanner_error(s, "unable to read matrix A"); goto STOP;
	}
	a_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else if (!strcmp(s->id, "B")) {/*Output probability*/
      if ((!n_read)||(!m_read)){
	scanner_error(s, "need M and N as a range for B"); goto STOP;
      }
      if (b_read) {scanner_error(s, "identifier B twice"); goto STOP;}
      if (!m_calloc(b_matrix, mo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if(matrix_d_read(s, b_matrix, mo->N, mo->M)) {
	  scanner_error(s, "unable to read matrix B"); goto STOP;
	}
	b_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else if (!strcmp(s->id, "prior")) {/*A prior model*/
      if (prior_read) {scanner_error(s,"identifier prior twice");goto STOP;}
      mo->prior = scanner_get_edouble(s);
      if ((mo->prior < 0 || mo->prior > 1) && mo->prior != -1)
	{ scanner_error(s, "invalid model prior"); goto STOP; }
      prior_read = 1;
    }   
    else if (!strcmp(s->id, "Pi")) {/*Initial state probabilty*/
      if (!n_read) {scanner_error(s, "need N as a range for Pi"); goto STOP;}
      if (pi_read) {scanner_error(s, "identifier Pi twice"); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "vector")) {
	scanner_consume( s, '{' ); if(s->err) goto STOP; 
	pi_vector = scanner_get_double_earray(s, &len);
	if (len != mo->N) {scanner_error(s, "wrong number of elements in PI"); goto STOP;}
	scanner_consume( s, ';' ); if(s->err) goto STOP;
	scanner_consume( s, '}' ); if(s->err) goto STOP; 
	pi_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else if (!strcmp(s->id, "fix_state")) {
      if (!n_read) {scanner_error(s, "need N as a range for fix_state"); goto STOP;}
      if (fix_read) {scanner_error(s, "identifier fix_state twice"); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "vector")) {
	scanner_consume( s, '{' ); if(s->err) goto STOP; 
	fix_vector = scanner_get_int_array(s, &len);
	if (len != mo->N) 
	  {scanner_error(s, "wrong number of elements in fix_state"); goto STOP;}
	scanner_consume( s, ';' ); if(s->err) goto STOP;
	scanner_consume( s, '}' ); if(s->err) goto STOP; 
	fix_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else {
      scanner_error(s, "unknown identifier");
      goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  } /* while(..s->c-'}') */
  scanner_consume( s, '}' ); if(s->err) goto STOP;  

  /* No prior read --> give it the value -1 */
  if (prior_read == 0)
    mo->prior = -1;
  /* Allocate memory for the model */
  for (i = 0; i < mo->N; i++) {
    mo->s[i].out_states = matrix_d_notzero_columns(a_matrix, i, mo->N);
    mo->s[i].in_states = matrix_d_notzero_rows(a_matrix, i, mo->N);
    if (model_state_alloc(mo->s + i, mo->M, mo->s[i].in_states,
			  mo->s[i].out_states)) { mes_proc(); goto STOP; }

    /* Assign the parameters to the model */
    if (!a_matrix) {
      fprintf(stderr,"no A matrix specified in file!\n");
      exit(1);
    }
    if (!b_matrix) {
      fprintf(stderr,"no B matrix specified in file!\n");
      exit(1);
    }
    if (!fix_vector) {
      fprintf(stderr,"no fix_state vector specified in file!\n");
      exit(1);
    }
    if (!pi_vector) {
      fprintf(stderr,"no Pi vector specified in file!\n");
      exit(1);
    }

    if(model_copy_vectors(mo, i, a_matrix, b_matrix, pi_vector, fix_vector)) {
      mes_proc(); goto STOP;
    }
  }
  matrix_d_free(&a_matrix, mo->N);
  matrix_d_free(&b_matrix, mo->N );
  m_free(pi_vector);
  return(mo);
STOP:
  matrix_d_free(&a_matrix, mo->N);
  matrix_d_free(&b_matrix, mo->N);
  m_free(pi_vector);
  model_free(&mo);
  return NULL;
#undef CUR_PROC
} /* model_direct_read */

/*============================================================================*/
/* Produces models from given sequences */
model **model_from_sequence(sequence_t *sq, long *mo_number){
#define CUR_PROC "model_from_sequence"
  long i;
  int max_symb;
  model **mo = NULL;
  if (!m_calloc(mo, sq->seq_number)) { mes_proc(); goto STOP; }
  max_symb = sequence_max_symbol(sq);
  for (i = 0; i < sq->seq_number; i++) 
    mo[i] = model_generate_from_sequence(sq->seq[i], sq->seq_len[i],
					 max_symb + 1);
  *mo_number = sq->seq_number;
  return mo;
 STOP:
  for (i = 0; i < *mo_number; i++)
    model_free(&(mo[i]));
  return NULL;
#undef CUR_PROC
} /* model_from_sequence */

/*============================================================================*/
/* Produces models form given sequences */
model **model_from_sequence_ascii(scanner_t *s, long *mo_number){
#define CUR_PROC "model_from_sequence_ascii"
  long i;
  int max_symb;
  model **mo = NULL;
  sequence_t *sq = NULL;

  scanner_consume( s, '{' ); if(s->err) goto STOP; 
  while(!s->err && !s->eof && s->c - '}') {    
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP; 
    /* Reads sequences on normal format */
    if (!strcmp(s->id, "SEQ")) {
      sq = sequence_read_alloc(s);
      if (!sq) { mes_proc(); goto STOP; }
    }
    else {
      scanner_error(s, "unknown identifier");
      goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  } /* while(..s->c-'}') */
  scanner_consume( s, '}' ); if(s->err) goto STOP;  

  if (!m_calloc(mo, sq->seq_number)) { mes_proc(); goto STOP; }
  /* The biggest symbol that occurs */
  max_symb = sequence_max_symbol(sq);
  for (i = 0; i < sq->seq_number; i++)
    mo[i] = model_generate_from_sequence(sq->seq[i], sq->seq_len[i],
					 max_symb + 1);

  *mo_number = sq->seq_number;
  sequence_free(&sq);
  return mo;
 STOP:
  sequence_free(&sq);
  for (i = 0; i < *mo_number; i++)
    model_free(&(mo[i]));
  return NULL;
#undef CUR_PROC
} /* model_from_sequence_ascii */

/*============================================================================*/

int model_free(model **mo) {
#define CUR_PROC "model_free"
  int i;
  mes_check_ptr(mo, return(-1));
  if( !*mo ) return(0);
  for (i = 0; i < (*mo)->N; i++)
    state_clean(&(*mo)->s[i]);
  
  if ( (*mo)->s)
    m_free((*mo)->s);
  if ((*mo) ->  silent)
    m_free((*mo)->silent);

  if ((*mo) -> tied_to)
    m_free((*mo)->tied_to);
  
  /*if ((*mo) -> emission_order) Moved to state
    m_free((*mo)->emission_order);*/
  
  if ((*mo) -> topo_order)
    m_free((*mo)->topo_order);
              
  m_free(*mo);
  return(0);
#undef CUR_PROC
} /* model_free */   


/*===========================================================================

int model_free(model **mo) {
#define CUR_PROC "sdmodel_free"
  state *my_state;
  int i;
  mes_check_ptr(mo, return(-1));
  if( !*mo ) return(0);
  for (i = 0; i < (*mo)->N; i++) {
    my_state = &((*mo)->s[i]);
    if (my_state->b)
      m_free(my_state->b);
    if (my_state->out_id)
      m_free(my_state->out_id);
    if (my_state->in_id)
      m_free(my_state->in_id);
    
    if (my_state->out_a)
      m_free(my_state->out_a);
    if (my_state->in_a)
      m_free(my_state->in_a);
    
	/*if (my_state->out_a)
	  matrix_d_free(&((*mo)->s[i].out_a), (*mo)->cos);
    if (my_state->in_a)
      matrix_d_free(&((*mo)->s[i].in_a), (*mo)->cos); 
    
    printf("Free every thing set it to NULL\n");

    my_state->pi         = 0;
    my_state->b          = NULL;
    my_state->out_id     = NULL;  
    my_state->in_id      = NULL;
    my_state->out_a      = NULL;
    my_state->in_a       = NULL;
    my_state->out_states = 0;
    my_state->in_states  = 0;
    my_state->fix        = 0;
  }
  m_free((*mo)->s);
  m_free(*mo);
  return(0);
#undef CUR_PROC
} /* model_free */   

/*============================================================================*/
model *model_copy(const model *mo) {
# define CUR_PROC "model_copy"
  int i, j, nachf, vorg, m;
  model *m2 = NULL;
  if(!m_calloc(m2, 1)) {mes_proc(); goto STOP;}
  if (!m_calloc(m2->s, mo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo->N; i++) {
    nachf = mo->s[i].out_states;
    vorg = mo->s[i].in_states;
    if(!m_calloc(m2->s[i].out_id, nachf)) {mes_proc(); goto STOP;}
    if(!m_calloc(m2->s[i].out_a, nachf)) {mes_proc(); goto STOP;}
    if(!m_calloc(m2->s[i].in_id, vorg)) {mes_proc(); goto STOP;}
    if(!m_calloc(m2->s[i].in_a, vorg)) {mes_proc(); goto STOP;}
    if(!m_calloc(m2->s[i].b, mo->M)) {mes_proc(); goto STOP;}
    /* Copy the values */
    for (j = 0; j < nachf; j++) {
      m2->s[i].out_a[j] = mo->s[i].out_a[j];
      m2->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      m2->s[i].in_a[j] = mo->s[i].in_a[j];
      m2->s[i].in_id[j] = mo->s[i].in_id[j];
    }
    for (m = 0; m < mo->M; m++)
      m2->s[i].b[m] = mo->s[i].b[m];
    m2->s[i].pi = mo->s[i].pi;
    m2->s[i].out_states = nachf;
    m2->s[i].in_states = vorg;
  }
  m2->N = mo->N;
  m2->M = mo->M;
  m2->prior = mo->prior;
  return(m2);
STOP:
  model_free(&m2);
  return(NULL);
# undef CUR_PROC
} /* model_copy */


/*============================================================================*/
int model_check(const model* mo) {
# define CUR_PROC "model_check"
  int res = -1;
  double sum;
  int i,j;
  char *str;
  /* The sum of the Pi[i]'s is 1 */

  sum = 0.0;
  for (i = 0; i < mo->N; i++) {
    sum += mo->s[i].pi;
  }

  if ( fabs(sum - 1.0) >= EPS_PREC )
    { mes_prot("sum Pi[i] != 1.0\n"); goto STOP; }

  /* check each state */
  for (i = 0; i < mo->N; i++) {
    sum = 0.0;
    if (mo->s[i].out_states == 0) {
      char *str = 
	mprintf(NULL,0,"out_states = 0 (state %d -> final state!)\n",i); 
      mes_prot(str);
    }
    /* Sum the a[i][j]'s : normalized out transitions */
    for (j = 0; j < mo->s[i].out_states; j++) {
      sum += mo->s[i].out_a[j];
      /* printf("    out_a[%d][%d] = %8.5f\n", i,j, mo->s[i].out_a[j]); */
    }
    if ( fabs(sum - 1.0) >= EPS_PREC ) { 
      char *str = mprintf(NULL, 0, "sum out_a[j] = %.2f != 1.0 (state %d)\n", 
			  sum, i); 
      mes_prot(str);
      m_free(str);
      /* goto STOP; */
    }
    /* Sum the b[j]'s: normalized emission probs */
    sum = 0.0;
    for (j = 0; j < mo->M; j++)
      sum += mo->s[i].b[j];
    if ( fabs(sum - 1.0) >= EPS_PREC ) { 
      res = -2; /* So we can ignore a silent state */
      str = mprintf(NULL, 0, "sum b[j] = %.2f != 1.0 (state %d)\n",
		    sum, i); 
      mes_prot(str);
      m_free(str);
      goto STOP;
    } /* i over all states */
  }
  
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* model_check */

/*============================================================================*/

int model_check_compatibility(model **mo, int model_number) {
#define CUR_PROC "model_check_compatibility"
  int i, j;
  for (i = 0; i < model_number; i++) 
    for (j = i+1; j < model_number; j++) {
      if (mo[i]->N != mo[j]->N) {
	char *str = 
	  mprintf(NULL, 0, "ERROR: different number of states in model %d (%d) and model %d (%d)", 
		  i, mo[i]->N, j, mo[j]->N); 
	mes_prot(str);
	m_free(str);
	return (-1);
      }
      if (mo[i]->M != mo[j]->M) {
	char *str = 
	  mprintf(NULL, 0, "ERROR: different number of possible outputs in model  %d (%d) and model %d (%d)", 
		  i, mo[i]->M, j, mo[j]->M); 
	mes_prot(str);
	m_free(str);
	return (-1);
      }
    }
  return 0;
#undef CUR_PROC    
} /* model_check_compatibility */

/*============================================================================*/

model *model_generate_from_sequence(const int *seq, int seq_len, int anz_symb) {
#define CUR_PROC "model_generate_from_sequence"
  int i;
  model *mo = NULL;
  state *s = NULL;
  if(!m_calloc(mo, 1)) {mes_proc(); goto STOP;}
  mo->N = seq_len; 
  mo->M = anz_symb;

  /* Allocate memory for all vectors */
  if (!m_calloc(mo->s, mo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo->N; i++) {
    if (i == 0) {               /* Initial state */
      if (model_state_alloc(mo->s, mo->M, 0, 1)) { mes_proc(); goto STOP; }
    }          
    else if (i == mo->N - 1) {  /* End state */
      if (model_state_alloc(mo->s + i, mo->M, 1, 0)) { mes_proc(); goto STOP; }
    }
    else {                      /* others */ 
      if (model_state_alloc(mo->s + i, mo->M, 1, 1)) { mes_proc(); goto STOP; }
    }
  }

  /* Allocate states with the right values, the initial state and the end 
     state extra */
  for (i = 1; i < mo->N - 1; i++) {
    s = mo->s + i;
    s->pi = 0.0;
    s->out_states = 1; s->in_states = 1;
    s->b[seq[i]] = 1.0; /* others stay 0 */    
    *(s->out_id) = i + 1;
    *(s->in_id) = i - 1;
    *(s->out_a) = *(s->in_a) = 1.0;
  }
    
  /* Initial state */
  s = mo->s;
  s->pi = 1.0; 
  s->out_states = 1; s->in_states = 0;
  s->b[seq[0]] = 1.0;
  *(s->out_id) = 1; 
  *(s->out_a) = 1.0;
  /* No in_id and in_a */

  /* End state */
  s = mo->s + mo->N - 1;
  s->pi = 0.0; 
  s->out_states = 0; s->in_states = 1;
  s->b[seq[mo->N-1]] = 1.0; /* All other b's stay zero */
  *(s->in_id) = mo->N - 2;
  *(s->in_a) = 1.0;
  /* No out_id and out_a */

  if (model_check(mo)) { mes_proc(); goto STOP; }
  return mo;
STOP:
  model_free(&mo);
  return NULL;
#undef CUR_PROC
} /* model_generate_from_sequence */


/*============================================================================*/
sequence_t *model_generate_sequences(model* mo, int seed, int global_len,
				     long seq_number, int Tmax) {
# define CUR_PROC "model_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  sequence_t *sq = NULL;
  int i, j, m,temp;
  double p, sum;
  int len = global_len;
  //int silent_len = 0;
  int n=0;
  int incomplete_seq = 0;
  int state=0;
  
  sq = sequence_calloc(seq_number);
  if (!sq) { mes_proc(); goto STOP; }
  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int)MAX_SEQ_LEN;
  
  if (seed > 0) {
	gsl_rng_set(RNG,seed);
  }
   
while (n < seq_number) {
	//printf("sequenz n = %d\n",n);
    //printf("incomplete_seq: %d\n",incomplete_seq);
	
    if (incomplete_seq == 0) { 
	    if(!m_calloc(sq->seq[n], len)) {mes_proc(); goto STOP;}
   		state = 0;
    }
   
	/* Get a random initial state i */
    p = gsl_rng_uniform(RNG);
	sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
	   break;
    }
    
  	if ( mo->model_type == kSilentStates ){
	  if (!mo->silent[i]) { /* first state emits */
	    //printf("first state %d not silent\n",i);
		  
		/* Get a random initial output m */
        p = gsl_rng_uniform(RNG);
	    sum = 0.0;   
        for (m = 0; m < mo->M; m++) {
          sum += mo->s[i].b[m];
          if (sum >= p)
	        break;
        }
          sq->seq[n][state] = m;
          state = state+1;
	  }
	  else {  /* silent state: we do nothing, no output */
		//printf("first state %d silent\n",i);
		//state = 0; 
		//silent_len = silent_len + 1; 
	  }
   }
	else {
	  /* Get a random initial output m */
	  p = gsl_rng_uniform(RNG);
	  sum = 0.0;   
	  for (m = 0; m < mo->M; m++) {
	    sum += mo->s[i].b[m];
		 if (sum >= p)
	      break;
	  }
	  sq->seq[n][state] = m;	
	  state = state + 1;
    }
	/* check whether sequence was completed by inital state*/ 
	if (state >= len && incomplete_seq == 1){
 	    
		//printf("assinging length %d to sequence %d\n",state,n);
		//printf("sequence complete...\n");
		sq->seq_len[n] = state;    
		incomplete_seq = 0;
		n++;
		continue;
	}
    while (state < len) {		
	  /* Get a random state i */
     p = gsl_rng_uniform(RNG);
	  sum = 0.0;   
     for (j = 0; j < mo->s[i].out_states; j++) {
	    sum += mo->s[i].out_a[j];   
	    if (sum >= p)
	      break;
     }

	//printf("state %d selected (i: %d, j: %d) at position %d\n",mo->s[i].out_id[j],i,j,state);
	 
	if (sum == 0.0) {
	    if (state < len) {
			incomplete_seq = 1;			
			//printf("final state reached - aborting\n");
			break;
		}
		else {
			n++;
			break;
		}
	 }	
     i = mo->s[i].out_id[j];
     if (mo->model_type == kSilentStates && mo->silent[i]) { /* Get a silent state i */
	    //printf("silent state \n");
        //silent_len += 1;
		/*if (silent_len >= Tmax) {
		   printf("%d silent states reached -> silent circle - aborting...\n",silent_len);
		   incomplete_seq = 0;
		   sq->seq_len[n] = state; 
		   n++; 
		   break;
		}*/
	  }  
	  else {
		/* Get a random output m from state i */
      	p = gsl_rng_uniform(RNG);
        sum = 0.0;   
        for (m = 0; m < mo->M; m++) {
	    	sum += mo->s[i].b[m];
	      	if (sum >= p)
	        	break;
        }
        sq->seq[n][state] = m;
        state++;
      }	
	  if (state == len ){
	     incomplete_seq = 0;
  	 	 sq->seq_len[n] = state;    
		 n++;
	  }  
	}  /* while (state < len) */  
 } /* while( n < seq_number )*/

return(sq);
STOP:
  sequence_free(&sq);
  return(NULL);
# undef CUR_PROC
} /* data */
  
/*============================================================================*/


double model_likelihood(model *mo, sequence_t *sq) {
# define CUR_PROC "model_likelihood"
  double log_p_i, log_p;
  int found, i,j;
	
  //printf("***  model_likelihood:\n");
  
  found = 0;
  log_p = 0.0;
  for (i = 0; i < sq->seq_number; i++) {
    
	//printf("sequence:\n");
	//for (j=0;j < sq->seq_len[i];j++) { 
	//	printf("%d, ",sq->seq[i][j]);
	//}
	//printf("\n"); 
	   
	  
	if (foba_logp(mo, sq->seq[i], sq->seq_len[i], &log_p_i)  == -1) {
	  mes_proc();
	  goto STOP;
	}
	
	//printf("\nlog_p_i = %f\n", log_p_i);
    
	if (log_p_i != +1) {
      log_p += log_p_i;
      found = 1;
    }
    else {
      char *str = mprintf(NULL, 0, "sequence[%d] can't be build.\n", i);
      mes_prot(str);
    }
  }
  if (!found)
    log_p = +1.0;
  return(log_p);
STOP:
  return -1;
# undef CUR_PROC
} /* model_likelihood */



void model_set_transition(model *mo, int i, int j, double prob) {
  # define CUR_PROC "model_set_transition"
  int in, out;
  
  if (mo->s && mo->s[i].out_a && mo->s[j].in_a) {
    for(out=0; out < mo->s[i].out_states; out++) {
      if ( mo->s[i].out_id[out] == j ) {
	mo->s[i].out_a[out] = prob;
	fprintf(stderr, "model_set_transition(0):State %d, %d, = %f\n", i, j, prob);
	break;
      }
    }

    for(in=0; in < mo->s[j].in_states; in++) {
      if ( mo->s[j].in_id[in] == i ) {
	mo->s[j].in_a[in] = prob;
	break;
      }
    }
  }
  # undef CUR_PROC
}
/* model_set_transition */




/*============================================================================*/
/* Some outputs */
/*============================================================================*/

void model_states_print(FILE *file, model *mo) {
  int i, j;
  fprintf(file, "Modelparameters: \n M = %d \t N = %d\n", 
	 mo->M, mo->N);
  for (i = 0; i < mo->N; i++) {
    fprintf(file, "\nState %d \n PI = %.3f \n out_states = %d \n in_states = %d \n",
	   i, mo->s[i].pi, mo->s[i].out_states, mo->s[i].in_states);
    fprintf(file, " Output probability:\t");
    for (j = 0; j < mo->M; j++)
      fprintf(file, "%.3f \t", mo->s[i].b[j]);
    fprintf(file, "\n Transition probability \n");
    fprintf(file, "  Out states (Id, a):\t");
    for (j = 0; j < mo->s[i].out_states; j++)
      fprintf(file, "(%d, %.3f) \t", mo->s[i].out_id[j], mo->s[i].out_a[j]);
    fprintf(file, "\n");
    fprintf(file, "  In states (Id, a):\t");
    for (j = 0; j < mo->s[i].in_states; j++)
      fprintf(file, "(%d, %.3f) \t", mo->s[i].in_id[j], mo->s[i].in_a[j]);
    fprintf(file, "\n");
  }
} /* model_states_print */

/*============================================================================*/

void model_A_print(FILE *file, model *mo, char *tab, char *separator, 
		   char *ending) {
  int i, j, out_state;
  for (i = 0; i < mo->N; i++) {
    out_state = 0;
    fprintf(file, "%s", tab);
    if (mo->s[i].out_states > 0 && 
	mo->s[i].out_id[out_state] == 0) {
	fprintf(file, "%.2f", mo->s[i].out_a[out_state]);
	out_state++;
    }
    else fprintf(file, "0.00");
    for (j = 1; j < mo->N; j++) {
      if (mo->s[i].out_states > out_state && 
	  mo->s[i].out_id[out_state] == j) {
	fprintf(file, "%s %.2f", separator, mo->s[i].out_a[out_state]);
	out_state++;
      }
      else fprintf(file, "%s 0.00", separator);
    }
    fprintf(file, "%s\n", ending);
  }
} /* model_A_print */

/*============================================================================*/

void model_B_print(FILE *file, model *mo, char *tab, char *separator, 
		   char *ending) {
  int i, j;
  for (i = 0; i < mo->N; i++) {
    fprintf(file, "%s", tab);
    fprintf(file, "%.2f", mo->s[i].b[0]);
    for (j = 1; j < mo->M; j++)
      fprintf(file, "%s %.2f", separator, mo->s[i].b[j]);
    fprintf(file, "%s\n", ending);
  }
} /* model_B_print */

/*============================================================================*/

void model_Pi_print(FILE *file, model *mo, char *tab, char *separator, 
		    char *ending) {
  int i;
  fprintf(file, "%s%.2f", tab, mo->s[0].pi);
  for (i = 1; i < mo->N; i++)
    fprintf(file, "%s %.2f", separator, mo->s[i].pi);
  fprintf(file, "%s\n", ending);
} /* model_Pi_print */

void model_fix_print(FILE *file, model *mo, char *tab, char *separator, 
		     char *ending) {
  int i;
  fprintf(file, "%s%d", tab, mo->s[0].fix);
  for (i = 1; i < mo->N; i++)
    fprintf(file, "%s %d", separator, mo->s[i].fix);
  fprintf(file, "%s\n", ending);
} /* model_Pi_print */

/*============================================================================*/

void model_A_print_transp(FILE *file, model *mo, char *tab, char *separator, 
			  char *ending) {
# define CUR_PROC "model_A_print_transp"
  int i, j;
  int *out_state;

  if(!m_calloc(out_state, mo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo->N; i++)
    out_state[i] = 0;

  for (j = 0; j < mo->N; j++) {
    fprintf(file, "%s", tab);
    if (mo->s[0].out_states != 0 && 
	mo->s[0].out_id[out_state[0]] == j) {
      fprintf(file, "%.2f", mo->s[0].out_a[out_state[0]]);
      (out_state[0])++;
    }
    else 
      fprintf(file, "0.00");
    for (i = 1; i < mo->N; i++) {
      if (mo->s[i].out_states != 0 && 
	  mo->s[i].out_id[out_state[i]] == j) {
	fprintf(file, "%s %.2f", separator, mo->s[i].out_a[out_state[i]]);
	(out_state[i])++;
      }
      else 
	fprintf(file, "%s 0.00", separator);
    }
    fprintf(file, "%s\n", ending);
  }
STOP:
  m_free(out_state);
  return;
# undef CUR_PROC
} /* model_A_print_transp */

/*============================================================================*/

void model_B_print_transp(FILE *file, model *mo, char *tab, char *separator, 
			  char *ending) {
  int i, j;
  for (j = 0; j < mo->M; j++) {
    fprintf(file, "%s", tab);
    fprintf(file, "%.2f", mo->s[0].b[j]);
    for (i = 1; i < mo->N; i++)
      fprintf(file, "%s %.2f", separator, mo->s[i].b[j]);
    fprintf(file, "%s\n", ending);
  }
} /* model_B_print_transp */

/*============================================================================*/

void model_Pi_print_transp(FILE *file, model *mo, char *tab, char *ending) {
  int i;
  for (i = 0; i < mo->N; i++)
    fprintf(file, "%s%.2f%s\n", tab, mo->s[i].pi, ending);
} /* model_Pi_print_transp */

/*============================================================================*/

void model_print(FILE *file, model *mo) {
  fprintf(file, "HMM = {\n\tM = %d;\n\tN = %d;\n", mo->M, mo->N);
  fprintf(file, "\tprior = %.3f;\n", mo->prior);
  fprintf(file, "\tA = matrix {\n");
  model_A_print(file, mo, "\t", ",", ";");
  fprintf(file, "\t};\n\tB = matrix {\n");
  model_B_print(file, mo,"\t", ",", ";");
  fprintf(file, "\t};\n\tPi = vector {\n");
  model_Pi_print(file, mo, "\t", ",", ";");
  fprintf(file, "\t};\n\tfix_state = vector {\n");
  model_fix_print(file, mo, "\t", ",", ";");
  fprintf(file, "\t};\n};\n\n");
} /* model_print */

/*============================================================================*/

void model_direct_print(FILE *file, model_direct *mo_d, int multip) {
  int i, j;
  for (i = 0; i < multip; i++) {
    fprintf(file, "HMM = {\n\tM = %d;\n\tN = %d;\n", mo_d->M, mo_d->N);
    fprintf(file, "\tprior = %.3f;\n", mo_d->prior);
    fprintf(file, "\tA = matrix {\n");
    matrix_d_print(file, mo_d->A, mo_d->N, mo_d->N, "\t", ",", ";");
    fprintf(file, "\t};\n\tB = matrix {\n");
    matrix_d_print(file, mo_d->B, mo_d->N, mo_d->M, "\t", ",", ";");
    fprintf(file, "\t};\n\tPi = vector {\n");
    fprintf(file, "\t%.4f", mo_d->Pi[0]);
    for (j = 1; j < mo_d->N; j++)
      fprintf(file, ", %.4f", mo_d->Pi[j]);
    fprintf(file, ";\n\t};\n");
    fprintf(file, "\tfix_state = vector {\n");
    fprintf(file, "\t%d", mo_d->fix_state[0]);
    for (j = 1; j < mo_d->N; j++)
      fprintf(file, ", %d", mo_d->fix_state[j]);
    fprintf(file, ";\n\t};\n");    
    fprintf(file, "};\n\n");
  }
} /* model_direct_print */

/*============================================================================*/

void model_direct_clean(model_direct *mo_d, hmm_check_t *check) {
#define CUR_PROC "model_direct_clean"
  int i;  
  if (!mo_d) return;
  mo_d->M = mo_d->N = 0;
  mo_d->prior = -1;
  if (mo_d->A) {
    for (i = 0; i < check->r_a; i++)
      m_free(mo_d->A[i]);
    m_free(mo_d->A);
  }
  if (mo_d->B) {
    for (i = 0; i < check->r_b; i++)
      m_free(mo_d->B[i]);
    m_free(mo_d->B);
  }
  if (mo_d->Pi) 
    m_free(mo_d->Pi);
  if (mo_d->fix_state)
    m_free(mo_d->fix_state);
      
  mo_d->A = mo_d->B = NULL;
  mo_d->Pi = NULL;
  mo_d->fix_state = NULL;
#undef CUR_PROC
} /* model_direct_clean */

/*============================================================================*/

int model_direct_check_data(model_direct *mo_d, hmm_check_t *check) {
#define CUR_PROC "model_direct_check_data"
  char *str;
  if (check->r_a != mo_d->N || check->c_a != mo_d->N) {
    str = mprintf(NULL, 0, "Incompatible dim. A (%d X %d) and N (%d)\n", 
		  check->r_a, check->c_a, mo_d->N); 
    mes_prot(str);
    m_free(str);
    return (-1);
  }
  if (check->r_b != mo_d->N || check->c_b != mo_d->M) {
    str = mprintf(NULL,0,"Incompatible dim. B (%d X %d) and N X M (%d X %d)\n",
		  check->r_b, check->c_b, mo_d->N, mo_d->M);
    mes_prot(str);
    m_free(str);
    return (-1);
  }
  if (check->len_pi != mo_d->N) {
    str = mprintf(NULL, 0, "Incompatible dim. Pi (%d) and N (%d)\n", 
		  check->len_pi, mo_d->N); 
    mes_prot(str);
    m_free(str);
    return (-1);
  }
  if (check->len_fix != mo_d->N) {
    str = mprintf(NULL, 0, "Incompatible dim. fix_state (%d) and N (%d)\n", 
		  check->len_fix, mo_d->N); 
    mes_prot(str);
    m_free(str);
    return (-1);
  }
  
  return 0;
#undef CUR_PROC
} /* model_direct_check_data */



/*============================================================================*/
/* XXX symmetric not implemented yet */
double model_prob_distance(model *m0, model *m, int maxT, int symmetric, 
			   int verbose)
{
#define CUR_PROC "model_prob_distance"

#define STEPS 40

  double p0, p;
  double d = 0.0;
  double *d1;
  sequence_t *seq0 = NULL;
  sequence_t *tmp = NULL;
  model *mo1, *mo2;
  int i, t, a, k;
  int true_len;
  int true_number;
  int left_to_right = 0;
  long total, index;
  int step_width = 0;
  int steps = 1;

  //printf("***  model_prob_distance:\n");
   
  if (verbose) { /* If we are doing it verbosely we want to have 40 steps */ 
    step_width = maxT / 40;
    steps = STEPS;
  }
  else         /* else just one */ 
    step_width = maxT;
  
  if( !m_calloc(d1, steps) ) {mes_proc();goto STOP;} 

  mo1 = m0;
  mo2 = m;
 
  for (k = 0; k < 2; k++) { /* Two passes for the symmetric case */
    
    /* seed = 0 -> no reseeding. Call  gsl_rng_timeseed(RNG) externally */
    seq0 = model_generate_sequences(mo1, 0, maxT+1, 1,maxT+1);
    
	
	
	if (seq0 == NULL) { mes_prot(" generate_sequences failed !");goto STOP;} 
		
    if (seq0->seq_len[0] < maxT) { /* There is an absorbing state */
      
      /* NOTA BENE: Assumpting the model delivers an explicit end state, 
		 the condition of a fix initial state is removed. */ 
      
      /* For now check that Pi puts all weight on state */
      /*
		t = 0;
		for (i = 0; i < mo1->N; i++) {
		if (mo1->s[i].pi > 0.001)
		t++;
		}    
		if (t > 1) {
		mes_prot("ERROR: No proper left-to-right model. Multiple start states");
		goto STOP;
		} */
      
      left_to_right = 1;
      total = seq0->seq_len[0];

      while (total <= maxT) {
	
		/* create a additional sequences at once */
		a = (maxT - total) / (total / seq0->seq_number) + 1;
		/* printf("total=%d generating %d", total, a); */
		tmp = model_generate_sequences(mo1, 0, 0, a,a);
    	if (tmp == NULL) { mes_prot(" generate_sequences failed !");goto STOP;} 
		sequence_free(&tmp);  
		sequence_add(seq0,tmp);
	
		total = 0;
		for (i = 0; i < seq0->seq_number; i++)
		  total += seq0->seq_len[i];      
   	   }
   	 }
    
   	 if (left_to_right) {
      
   	   for (t=step_width, i=0; t <= maxT; t+= step_width, i++) {

		index = 0;
		total = seq0->seq_len[0];
	
		/* Determine how many sequences we need to get a total of t
		   and adjust length of last sequence to obtain total of 
		   exactly t */

		while(total < t) {
		  index++;
		  total += seq0->seq_len[index];      
		}

		true_len = seq0->seq_len[index];
		true_number = seq0->seq_number;
	
		if ((total - t) > 0)
		  seq0->seq_len[index] = total - t;
		seq0->seq_number = index;
	
		p0 = model_likelihood(mo1, seq0);
		if (p0 == +1 || p0 == -1) { /* error! */
		  mes_prot("problem: model_likelihood failed !");
		  goto STOP;
		}	  
		p = model_likelihood(mo2, seq0);
		if (p == +1 || p == -1) { /* what shall we do now? */
		  mes_prot("problem: model_likelihood failed !");
		  goto STOP;
		}	  

		d = 1.0 / t * (p0 -  p);

		if (symmetric) {
		  if (k == 0)
		    /* save d */
	  	  d1[i] = d;
		  else {
		    /* calculate d */
		    d = 0.5 * (d1[i] + d);
		  }  
		}

		if (verbose && (!symmetric || k == 1))
		  printf("%d\t%f\t%f\t%f\n", t, p0, p, d);
	
		seq0->seq_len[index] = true_len;
		seq0->seq_number = true_number;
      }
    } 

    else {
      
      true_len = seq0->seq_len[0];

      for (t=step_width, i=0; t <= maxT; t+= step_width, i++) {
		seq0->seq_len[0] = t;
	
	p0 = model_likelihood(mo1, seq0);
	//printf("   P(O|m1) = %f\n",p0);
	if (p0 == +1) { /* error! */
	 	mes_prot("seq0 can't be build from mo1!");
	  goto STOP;
	}	  
	p = model_likelihood(mo2, seq0);
	//printf("   P(O|m2) = %f\n",p);
	if (p == +1) { /* what shall we do now? */
	    mes_prot("problem: seq0 can't be build from mo2!");
	  goto STOP;
	}	  

	d = (1.0 / t) * (p0 -  p);

	if (symmetric) {
	  if (k == 0)
	    /* save d */
	    d1[i] = d;
	  else {
	    /* calculate d */
	    d = 0.5 * (d1[i] + d);
	  }  
	}

	if (verbose && (!symmetric || k == 1))
	  printf("%d\t%f\t%f\t%f\n", t, p0, p, d);
	
      }
      seq0->seq_len[0] = true_len;
    }

    if (symmetric) {
      sequence_free(&seq0);
      mo1 = m ;
      mo2 = m0;
    }
    else
      break;
    
  } /* k = 1,2 */ 
  
  sequence_free(&seq0);
  free(d1);    
  return d;

STOP:
  sequence_free(&seq0);
  free(d1);
  return(0.0);
#undef CUR_PROC
}


/*============================================================================*/

void state_clean(state *my_state) {
#define CUR_PROC "state_clean"
  if (!my_state) return;
  
  if (my_state->b)
    m_free(my_state->b);
  
  if (my_state->out_id)
    m_free(my_state->out_id);
  
  if (my_state->in_id)
    m_free(my_state->in_id);
  
  if (my_state->out_a)
    m_free(my_state->out_a);
  
  if (my_state->in_a)
    m_free(my_state->in_a);

  my_state->pi         = 0;
  my_state->b          = NULL;
  my_state->out_id     = NULL;  
  my_state->in_id      = NULL;
  my_state->out_a      = NULL;
  my_state->in_a       = NULL;
  my_state->out_states = 0;
  my_state->in_states  = 0;
  my_state->fix        = 0;  

#undef CUR_PROC
} /* state_clean */

/*============================================================================*/

/*state* state_copy(state *my_state) {
  state* new_state = (state*) malloc(sizeof(state));

  state_copy_to(my_state,new_state);

  return new_state;

  }*/ /* state_copy */

/*============================================================================*/

/*void state_copy_to(state *source, state* dest) {
  dest->pi         = source->pi;
  dest->out_states = source->out_states;
  dest->in_states  = source->in_states;
  dest->fix        = source->fix;

  dest->b          = malloc(xxx);
  memcpy(dest->b,source->b,xxx);

  dest->out_id     = malloc(xxx);
  memcpy(dest->out_id,source->out_id,xxx);

  dest->in_id      = malloc(xxx);
  memcpy(dest->in_id,source->in_id,xxx);

  dest->out_a      = malloc(xxx);
  memcpy(dest->out_a,source->out_a,xxx);

  dest->in_a       = malloc(xxx);
  memcpy(dest->in_a,source->in_a,xxx);
  }*/ /* state_copy_to */


 /*==========================Labeled HMMS ================================*/         

           
sequence_t *model_label_generate_sequences(model* mo, int seed, int global_len,
				     long seq_number, int Tmax) {
# define CUR_PROC "model_label_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  sequence_t *sq = NULL;
  int i, j, m,temp;
  double p, sum;
  int len = global_len;
  //int silent_len = 0;
  int n=0;
  int incomplete_seq = 0;
  int state=0;
  int label_index = 0;
  
  sq = sequence_calloc(seq_number);

  if (!sq) { mes_proc(); goto STOP; }

  /* allocating additional fields for the labels in the sequence_t struct */
  if(!m_calloc(sq->state_labels, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->state_labels_len, seq_number)) {mes_proc(); goto STOP;}

  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int)MAX_SEQ_LEN;
  
  if (seed > 0) {
	gsl_rng_set(RNG,seed);
  }
   

while (n < seq_number) {
	//printf("sequenz n = %d\n",n);
    //printf("incomplete_seq: %d\n",incomplete_seq);
	
    if (incomplete_seq == 0) { 
	    if(!m_calloc(sq->seq[n], len)) {mes_proc(); goto STOP;}
   		if (mo->model_type == kSilentStates){
          
          printf("Model has silent states. \n");  
          /* for silent models we have to allocate for the maximal possible number of lables*/
          if(!m_calloc(sq->state_labels[n], len * mo->N)) {mes_proc(); goto STOP;}
        }    
        else {
            printf("Model has no silent states. \n");
            if(!m_calloc(sq->state_labels[n], len)) {mes_proc(); goto STOP;}
        }
        label_index = 0;
        state = 0;
    }
   
	/* Get a random initial state i */
    p = gsl_rng_uniform(RNG);
	sum = 0.0;
    for (i = 0; i < mo->N; i++) {
      sum += mo->s[i].pi;
      if (sum >= p)
	   break;
    }
    
    /* add label of fist state to the label list */
    sq->state_labels[n][label_index] = mo->s[i].label;
    label_index++;
    
  	if ( mo->model_type == kSilentStates ){
	  if (!mo->silent[i]) { /* first state emits */
	    //printf("first state %d not silent\n",i);
		  
		/* Get a random initial output m */
        p = gsl_rng_uniform(RNG);
	    sum = 0.0;   
        for (m = 0; m < mo->M; m++) {
          sum += mo->s[i].b[m];
          if (sum >= p)
	        break;
        }
          sq->seq[n][state] = m;
          state = state+1;
	  }
	  else {  /* silent state: we do nothing, no output */
		//printf("first state %d silent\n",i);
		//state = 0; 
		//silent_len = silent_len + 1; 
	  }
   }
	else {
	  /* Get a random initial output m */
	  p = gsl_rng_uniform(RNG);
	  sum = 0.0;   
	  for (m = 0; m < mo->M; m++) {
	    sum += mo->s[i].b[m];
		 if (sum >= p)
	      break;
	  }
	  sq->seq[n][state] = m;	
	  state = state + 1;
    }
	/* check whether sequence was completed by inital state*/ 
	if (state >= len && incomplete_seq == 1){
 	    
		//printf("assinging length %d to sequence %d\n",state,n);
		//printf("sequence complete...\n");
		
        sq->seq_len[n] = state;    
	
        sq->state_labels_len[n] = label_index;
		
        printf("1: seq %d -> %d labels\n",n,sq->state_labels_len[n]);
        
        if (mo->model_type == kSilentStates){
          printf("reallocating\n");
          if (m_realloc(sq->state_labels[n], sq->state_labels_len[n])){mes_proc(); goto STOP;}
        }  
        
  		incomplete_seq = 0;
		n++;
        continue;
	}
    while (state < len) {		
	  /* Get a random state i */
     p = gsl_rng_uniform(RNG);
	  sum = 0.0;   
     for (j = 0; j < mo->s[i].out_states; j++) {
	    sum += mo->s[i].out_a[j];   
	    if (sum >= p)
	      break;
     }

     i = mo->s[i].out_id[j];
    
    /* add label of state to the label list */
    sq->state_labels[n][label_index] = mo->s[i].label;
    label_index++;	
     
    //printf("state %d selected (i: %d, j: %d) at position %d\n",mo->s[i].out_id[j],i,j,state);
	 
	if (sum == 0.0) {
	    if (state < len) {
			incomplete_seq = 1;			
			//printf("final state reached - aborting\n");
			break;
		}
		else {
			n++;
			break;
		}
	 }	
     
     if (mo->model_type == kSilentStates && mo->silent[i]) { /* Got a silent state i */
	    //printf("silent state \n");
        //silent_len += 1;
		/*if (silent_len >= Tmax) {
		   printf("%d silent states reached -> silent circle - aborting...\n",silent_len);
		   incomplete_seq = 0;
		   sq->seq_len[n] = state; 
		   n++; 
		   break;
		}*/
	  }  
	  else {
		/* Get a random output m from state i */
      	p = gsl_rng_uniform(RNG);
        sum = 0.0;   
        for (m = 0; m < mo->M; m++) {
	    	sum += mo->s[i].b[m];
	      	if (sum >= p)
	        	break;
        }
        sq->seq[n][state] = m;
        state++;
      }	
	  if (state == len ){
	     incomplete_seq = 0;
  	 	 
         sq->state_labels_len[n] = label_index;
		
         printf("2: seq %d -> %d labels\n",n,sq->state_labels_len[n]);
        
         if (mo->model_type == kSilentStates){
           printf("reallocating\n");
           if (m_realloc(sq->state_labels[n], sq->state_labels_len[n])){mes_proc(); goto STOP;}
         } 
         
         sq->seq_len[n] = state;    
		 n++;
	  }  
	}  /* while (state < len) */  
 } /* while( n < seq_number )*/

return(sq);
STOP:
  sequence_free(&sq);
  return(NULL);
# undef CUR_PROC
} /* data */
  
/*============================================================================*/
          
 /*===================== E n d   o f  f i l e  "model.c"       ===============*/
