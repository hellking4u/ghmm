/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/model.c
  created      : TIME: 10:47:27     DATE: Fri 19. December 1997
  $Id$

__copyright__

*******************************************************************************/


#include <float.h>
#include <math.h>
#include "sdmodel.h"
#include "matrix.h"
#include "sequence.h"
#include "const.h"
#include "rng.h"
#include "foba.h"
#include "mes.h"
#include "mprintf.h"
#include "string.h"

/*----------------------------------------------------------------------------*/
static int sdmodel_state_alloc(sdstate *state, int M, int in_states,
			       int out_states, int cos) {
# define CUR_PROC "sdmodel_state_alloc"
  int res = -1;
  if(!m_calloc(state->b, M)) {mes_proc(); goto STOP;}

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
} /* model_state_alloc */

/*----------------------------------------------------------------------------*/

static int sdmodel_copy_vectors(sdmodel *mo, int index, double ***a_matrix, 
			      double **b_matrix, double *pi, int *fix) {
#define CUR_PROC "sdmodel_alloc_vectors"
  int i, j, c, cnt_out = 0, cnt_in = 0;

  mo->s[index].pi = pi[index];
  mo->s[index].fix = fix[index];
  for (i = 0; i < mo->M; i++)
    mo->s[index].b[i] = b_matrix[index][i];
  for (c = 0; c < mo->cos; c++)
    for (i = 0; i < mo->N; i++) {
      if (a_matrix[c][index][i]) { /* Transitions to a following state possible */
	if (cnt_out >= mo->s[index].out_states) {mes_proc(); return(-1);}
	mo->s[index].out_id[cnt_out] = i;
	mo->s[index].out_a[c][cnt_out] = a_matrix[c][index][i];
	cnt_out++;
      }
      if (a_matrix[i][index]) { /* Transitions to a previous state possible */
	if (cnt_in >= mo->s[index].in_states) {mes_proc(); return(-1);}
	mo->s[index].in_id[cnt_in] = i;
	mo->s[index].in_a[c][cnt_in] = a_matrix[c][i][index];
	cnt_in++;
      }
    }
  return(0);
#undef CUR_PROC
} /* model_alloc_vectors */


/*============================================================================*/

/* Old prototyp:

model **model_read(char *filename, int *mo_number, int **seq,
			 const int *seq_len, int seq_number) { */



/*============================================================================*/

int sdmodel_free(sdmodel **mo) {
#define CUR_PROC "sdmodel_free"
  sdstate *my_state;
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
    /*
      if (my_state->out_a)
      m_free(my_state->out_a);
      if (my_state->in_a)
      m_free(my_state->in_a);*/
    if (my_state->out_a)
      matrix_d_free(&((*mo)->s[i].out_a), (*mo)->cos);
    if (my_state->in_a)
      matrix_d_free(&((*mo)->s[i].in_a), (*mo)->cos);
    
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
sdmodel *sdmodel_copy(const sdmodel *mo) {
# define CUR_PROC "sdmodel_copy"
  int i, j, k, nachf, vorg, m;
  sdmodel *m2 = NULL;
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
      for (k = 0; k < mo->cos; k++)
	m2->s[i].out_a[k][j] = mo->s[i].out_a[k][j];
      m2->s[i].out_id[j] = mo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      for (k = 0; k < mo->cos; k++)
	m2->s[i].in_a[k][j] = mo->s[i].in_a[k][j];
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
  sdmodel_free(&m2);
  return(NULL);
# undef CUR_PROC
} /* model_copy */



/*============================================================================*/
sequence_t *sdmodel_generate_sequences(sdmodel* mo, int seed, int global_len,
				     long seq_number, int Tmax) {
# define CUR_PROC "sdmodel_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  sequence_t *sq = NULL;
  int state, n, i, j, m,  reject_os, reject_tmax, badseq, class;
  double p, sum, osum = 0.0;
  int len = global_len, up = 0, stillbadseq = 0, reject_os_tmp = 0;
  double dummy = 0.0;

  sq = sequence_calloc(seq_number);
  if (!sq) { mes_proc(); goto STOP; }
  if (len <= 0)
    /* A specific length of the sequences isn't given. As a model should have
       an end state, the konstant MAX_SEQ_LEN is used. */
    len = (int)MAX_SEQ_LEN;
  
  if (seed > 0)
    gsl_rng_set(RNG,seed);

  // for (n = 0; n < seq_number; n++) {

  n = 0;
  reject_os = reject_tmax = 0;

  while (n < seq_number) 
    {
      /* Test: A new seed for each sequence */
      /*   gsl_rng_timeseed(RNG); */
      stillbadseq = badseq = 0;
      if(!m_calloc(sq->seq[n], len)) {mes_proc(); goto STOP;}

      /* Get a random initial state i */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;
      for (i = 0; i < mo->N; i++) {
	sum += mo->s[i].pi;
	if (sum >= p)
	  break;
      }

      /* Get a random initial output m */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (m = 0; m < mo->M; m++) {
	sum += mo->s[i].b[m];
	if (sum >= p)
	  break;
      }
      sq->seq[n][0] = m;
      state = 1;

      /* The first symbol chooses the start class */
      class = sequence_d_class(&dummy, 0, &osum); /*  dummy function */
      while (state < len) {
      
	/* Get a new state */
	p = gsl_rng_uniform(RNG);
	sum = 0.0;   
	for (j = 0; j < mo->s[i].out_states; j++) {
	  sum += mo->s[i].out_a[class][j];   
	  if (sum >= p)
	    break;
	}
      
	if (sum == 0.0) {
	  if (mo->s[i].out_states > 0) {
	    /* Repudiate the sequence, if all smo->s[i].out_a[class][.] == 0,
	       that is, class "class" isn't used in the original data:
	       go out of the while-loop, n should not be counted. */
	    /* printf("Zustand %d, class %d, len %d out_states %d \n", i, class,
	       state, smo->s[i].out_states); */
	    badseq = 1;
	    /* break; */
	  
	    /* Try: If the class is "empty", try the neighbour class;
	       first, sweep down to zero; if still no success, sweep up to
	       COS - 1. If still no success --> Repudiate the sequence. */
	    if (class > 0 && up == 0) {
	      class--;
	      continue;
	    }
	    else if (class < mo->cos - 1) {
	      class++;
	      up = 1;
	      continue;
	    }
	    else {
	      stillbadseq = 1;
	      break;
	    }
	  } else

	    /* An end state is reached, get out of the while-loop */
	    break;
	}
	i = mo->s[i].out_id[j];

	/* Get a random output m from state i */
	p = gsl_rng_uniform(RNG);
	sum = 0.0;   
	for (m = 0; m < mo->M; m++) {
	  sum += mo->s[i].b[m];
	  if (sum >= p)
	    break;
	}
      
	sq->seq[n][state] = m;

	/* Decide the class for the next step */
	class = sequence_d_class(&dummy, state, &osum); /* dummy */
	up = 0;
	state++;
      } /* while (state < len) , global_len depends on the data */  
      
      if (badseq) {
	reject_os_tmp++;
      }
      
      if (stillbadseq) {
	reject_os++;
	m_free(sq->seq[n]);
	/*      printf("cl %d, s %d, %d\n", class, i, n); */
      }
      else if (state > Tmax) {
	reject_tmax++;
	m_free(sq->seq[n]);
      }
      else {
	if (state < len)
	  if(m_realloc(sq->seq[n], state)) {mes_proc(); goto STOP;}
	sq->seq_len[n] = state;
	/* sq->seq_label[n] = label; */
	/* vector_d_print(stdout, sq->seq[n], sq->seq_len[n]," "," ",""); */
	n++;
      }
      /*    printf("reject_os %d, reject_tmax %d\n", reject_os, reject_tmax); */
      if (reject_os > 10000) {
	mes_prot("Reached max. no. of rejections\n");
	break;
      }
      if (!(n%1000)) printf("%d Seqs. generated\n", n);
    } /* n-loop */
  

  if (reject_os > 0) printf("%d sequences rejected (os)!\n", reject_os);
  if (reject_os_tmp > 0) printf("%d sequences changed class\n", 
				reject_os_tmp - reject_os);	
  if (reject_tmax > 0) printf("%d sequences rejected (Tmax, %d)!\n", 
			      reject_tmax, Tmax);

  return(sq);
STOP:
  sequence_free(&sq);
  return(NULL);
# undef CUR_PROC
} /* data */
  

/*============================================================================*/
/* Some outputs */
/*============================================================================*/

void sdmodel_states_print(FILE *file, sdmodel *mo) {
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

void sdmodel_Ak_print(FILE *file, sdmodel *mo, int k, char *tab, char *separator, 
		   char *ending) {
  int i, j, out_state;
  for (i = 0; i < mo->N; i++) {
    out_state = 0;
    fprintf(file, "%s", tab);
    if (mo->s[i].out_states > 0 && 
	mo->s[i].out_id[out_state] == 0) {
	fprintf(file, "%.2f", mo->s[i].out_a[k][out_state]);
	out_state++;
    }
    else fprintf(file, "0.00");
    for (j = 1; j < mo->N; j++) {
      if (mo->s[i].out_states > out_state && 
	  mo->s[i].out_id[out_state] == j) {
	fprintf(file, "%s %.2f", separator, mo->s[i].out_a[k][out_state]);
	out_state++;
      }
      else fprintf(file, "%s 0.00", separator);
    }
    fprintf(file, "%s\n", ending);
  }
} /* model_A_print */

/*============================================================================*/

void sdmodel_B_print(FILE *file, sdmodel *mo, char *tab, char *separator, 
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

void sdmodel_Pi_print(FILE *file, sdmodel *mo, char *tab, char *separator, 
		    char *ending) {
  int i;
  fprintf(file, "%s%.2f", tab, mo->s[0].pi);
  for (i = 1; i < mo->N; i++)
    fprintf(file, "%s %.2f", separator, mo->s[i].pi);
  fprintf(file, "%s\n", ending);
} /* model_Pi_print */


/*============================================================================*/

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

/*===================== E n d   o f  f i l e  "model.c"       ===============*/
