/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/reestimate.c
  created      : TIME: 11:54:33     DATE: Mon 16. February 1998
  $Id$

__copyright__

*******************************************************************************/


#include "mes.h"
#include "mprintf.h"
#include "reestimate.h"
#include "matrix.h"
#include "model.h"
#include "foba.h"
#include "float.h"
#include "const.h"
#include "math.h"

typedef struct local_store_t {
  double *pi_num;
  double pi_denom;
  double **a_num;
  double *a_denom;
  double **b_num;
  double *b_denom;
} local_store_t;

static local_store_t *reestimate_alloc(const model *mo);
static int reestimate_free(local_store_t **r, int N);
static int reestimate_init(local_store_t *r, const model *mo);
static int reestimate_alloc_matvek(double ***alpha, double ***beta, 
				   double **scale, int T, int N);
static int reestimate_free_matvek(double **alpha, double **beta, 
				  double *scale, int T);
static int reestimate_setlambda(local_store_t *r, model *mo);
static int reestimate_one_step(model *mo, local_store_t *r, 
			       int seq_number, int *seq_length, int **O,
			       double *log_p, double *seq_w);

/*----------------------------------------------------------------------------*/
static local_store_t *reestimate_alloc(const model *mo) {
# define CUR_PROC "reestimate_alloc"
  int i;
  local_store_t* r = NULL;
  if (!m_calloc(r, 1)) {mes_proc(); goto STOP;}
  if (!m_calloc(r->pi_num, mo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(r->a_num, mo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo->N; i++)
    if (!m_calloc(r->a_num[i], mo->s[i].out_states)){mes_proc(); goto STOP;}
  if (!m_calloc(r->a_denom, mo->N)) {mes_proc(); goto STOP;}
  r->b_num = matrix_d_alloc(mo->N, mo->M);
  if (!(r->b_num)) {mes_proc(); goto STOP;}
  if (!m_calloc(r->b_denom, mo->N)) {mes_proc(); goto STOP;}
  return(r);
STOP:
  reestimate_free(&r, mo->N);
  return(NULL);
# undef CUR_PROC
} /* reestimate_alloc */

/*----------------------------------------------------------------------------*/
static int reestimate_free(local_store_t **r, int N) {
# define CUR_PROC "reestimate_free"
  int i;
  mes_check_ptr(r, return(-1));
  if( !*r ) return(0);
  m_free((*r)->pi_num);  
  for (i = 0; i < N; i++)
    m_free((*r)->a_num[i]);
  m_free((*r)->a_num);
  m_free((*r)->a_denom);
  matrix_d_free( &((*r)->b_num), N );
  m_free((*r)->b_denom);
  m_free(*r);
  return(0);
# undef CUR_PROC
} /* reestimate_free */

/*----------------------------------------------------------------------------*/
static int reestimate_init(local_store_t *r, const model *mo) {
# define CUR_PROC "reestimate_init"
  int i, j, m;
  r->pi_denom = 0.0;
  for (i = 0; i < mo->N; i++) {
    r->pi_num[i] = 0.0;
    r->a_denom[i] = 0.0;
    r->b_denom[i] = 0.0;
    for (j = 0; j < mo->s[i].out_states; j++)
      r->a_num[i][j] = 0.0;
    for (m = 0; m < mo->M; m++)
      r->b_num[i][m] = 0.0;
  }
  return(0);
# undef CUR_PROC
} /* reestimate_init */

/*----------------------------------------------------------------------------*/
static int reestimate_alloc_matvek(double ***alpha, double ***beta, 
				   double **scale, int T, int N) {
# define CUR_PROC "reestimate_alloc_matvek"
  int res = -1;
  *alpha = matrix_d_alloc(T, N);
  if (!(*alpha)) {mes_proc(); goto STOP;}
  *beta = matrix_d_alloc(T, N);
  if (!(*beta)) {mes_proc(); goto STOP;}
  if (!m_calloc(*scale, T)) {mes_proc(); goto STOP;}
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* reestimate_alloc_matvek */

/*----------------------------------------------------------------------------*/
static int reestimate_free_matvek(double **alpha, double **beta, 
				  double *scale, int T) {
# define CUR_PROC "reestimate_free_matvek"
  matrix_d_free(&alpha, T);
  matrix_d_free(&beta, T);
  m_free(scale); 
  return(0);
# undef CUR_PROC
} /* reestimate_free_matvek */   

/*----------------------------------------------------------------------------*/
static int reestimate_setlambda(local_store_t *r, model *mo) {
# define CUR_PROC "reestimate_setlambda"
  int res = -1;
  int h, i, j, m, l, j_id, positive;
  double factor, p_i;
  mes_check_0(r->pi_denom, goto STOP); 
  for (i = 0; i < mo->N; i++) {
    /* Pi */
    mo->s[i].pi =  r->pi_num[i] / r->pi_denom;
    /* A */
    /* note: denom. might be 0; never reached state? */
    p_i = 0.0;
    if (r->a_denom[i] < EPS_PREC) {     
      for (h = 0; h < mo->s[i].in_states; h++) {
	p_i += mo->s[i].in_a[h];
      }	
      if (p_i == 0) {
	char *str;
	if (mo->s[i].in_states == 0)
	  str = mprintf(NULL,0,
			" State %d can't be reached (no in_states)\n",i);
	else
	  str = mprintf(NULL,0, " State %d can't be reached (prob = 0)\n",i);
	mes_prot(str);	m_free(str);
      }
      factor = 0.0;
    }
    else 
      factor = ( 1 / r->a_denom[i] );          
    positive = 0;
    
    for (j = 0; j < mo->s[i].out_states; j++) {
      /* TEST: denom. < numerator */
      if ((r->a_denom[i] - r->a_num[i][j]) <= -EPS_PREC)
	{ mes_prot(" numerator > denom!\n"); }
      mo->s[i].out_a[j] = r->a_num[i][j] * factor;
      if (r->a_num[i][j] >= EPS_PREC) positive = 1;
      /* important: also update in_a  */
            l = 0;
      j_id = mo->s[i].out_id[j];
      while (l < mo->s[j_id].in_states)
	if  (mo->s[j_id].in_id[l] == i)
	  break;
        else 
	  l++;
      if ( l == mo->s[j_id].in_states) { 
	mes_prot("no corresponding in_a to out_a found!\n");
	goto STOP;
      }
      mo->s[j_id].in_a[l] = mo->s[i].out_a[j]; 
    }

    if (!positive) {
      char *str = 
	mprintf(NULL, 0, 
		"All numerator a[%d][j] == 0 (denom=%.4f, P(in)=%.4f)!\n", 
		i, r->a_denom[i], p_i);
      mes_prot(str);
      m_free(str);
    }

    /* if fix, continue to next state */
    if (mo->s[i].fix)
      continue;
    
    /* B */
    if (r->b_denom[i] < EPS_PREC)
      factor = 0.0;
    else 
      factor = ( 1.0 / r->b_denom[i] );      

    positive = 0;

    for (m = 0; m < mo->M; m++) {
      /* TEST: denom. < numerator */
      if ((r->b_denom[i] - r->b_num[i][m]) <= -EPS_PREC) {
	char *str = 
	  mprintf(NULL, 0, "numerator b (%.4f) > denom (%.4f)!\n", 
		  r->b_num[i][m], r->b_denom[i]); 
	mes_prot(str);
	m_free(str);
      }
      mo->s[i].b[m] = r->b_num[i][m] * factor;

      if (r->b_num[i][m] >= EPS_PREC) positive = 1;
      
    }

    if (!positive) {
      char *str = 
	mprintf(NULL, 0, "All numerator b[%d][m] == 0 (denom = %.4f)!\n", 
		i, r->b_denom[i]);
      mes_prot(str);
      m_free(str);
    }
    
  } /* for (i = 0 .. < mo->N)  */
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* reestimate_setlambda */

/*----------------------------------------------------------------------------*/
static int reestimate_one_step(model *mo, local_store_t *r, 
			       int seq_number, int *seq_length, int **O,
			       double *log_p, double *seq_w) { 
# define CUR_PROC "reestimate_one_step"
  int res = -1;  
  int k, i, j, m, t, j_id, valid;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  int T_k;
  double gamma;
  double log_p_k;

  *log_p = 0.0;
  valid = 0;
  /* loop over all sequences */
  for (k = 0; k < seq_number; k++) {
    
    T_k = seq_length[k]; /* current seq. length */
    
    /* initialization of  matrices and vector depends on T_k */
    if ( reestimate_alloc_matvek(&alpha, &beta, &scale, T_k, mo->N) == -1 ) 
      {mes_proc(); goto STOP;}
    
     if (foba_forward(mo, O[k], T_k, alpha, scale, &log_p_k) == -1) 
      {mes_proc(); goto STOP;}    
     
     if (log_p_k != +1) { /* O[k] can be generated */
       *log_p += log_p_k;
       valid = 1;
      
      if (foba_backward(mo, O[k], T_k, beta, scale) == -1) 
	{ mes_proc(); goto STOP;}
      
      /* loop over all states */
      for (i = 0; i < mo->N; i++) {
	/* Pi */
	r->pi_num[i] += seq_w[k] * alpha[0][i] * beta[0][i];
	r->pi_denom += seq_w[k] * alpha[0][i] * beta[0][i]; 
	/* A */
	for (t = 0; t < T_k - 1; t++) {
	  r->a_denom[i] += seq_w[k] * alpha[t][i] * beta[t][i]; 
	  for (j = 0; j < mo->s[i].out_states; j++) {
	    j_id = mo->s[i].out_id[j]; /* aufpassen! */
	    r->a_num[i][j] += 
	      ( seq_w[k] * alpha[t][i] 
		* mo->s[i].out_a[j] 
		* mo->s[j_id].b[ O[k][t+1] ]
		* beta[t+1][j_id] 
		* (1.0 / scale[t+1]) );  /* c[t] = 1/scale[t] */
	  }
	}
	/* ========= if state fix, continue;====================== */
	if (mo->s[i].fix)
	  continue;
	/* B */
	for (t = 0; t < T_k; t++) {
	  gamma = (seq_w[k] * alpha[t][i] * beta[t][i] );
	  r->b_denom[i] += gamma;                  
	  for (m = 0; m < mo->M; m++) {
	    if ( O[k][t] == m )
	      r->b_num[i][m] += gamma;
	  }
	}
      } /* for (i = 0 i < mo->N i++) { */

    } /* if (log_p_k != +1) */
    else {
      printf("O(%2d) can't be build from model mo!\n", k);
    } 

    reestimate_free_matvek(alpha, beta, scale, T_k);
    
  } /* for (k = 0; k < seq_number; k++) */
  
  if (valid) {
    /* new parameter lambda: set directly in model */
    if ( reestimate_setlambda(r, mo) == -1 ) { mes_proc(); goto STOP; } 

    /* only test: */
    /*    if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  }
  else { /* NO sequence can be build from model mo! */
    *log_p = +1;
  }

  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* reestimate_one_step */

/*============================================================================*/
int reestimate_baum_welch(model *mo, sequence_t *sq) {
# define CUR_PROC "reestimate_baum_welch"
  int n, k, valid;
  double log_p, log_p_old, log_p_k, diff;
  local_store_t *r = NULL;
 
  /* local store for all iterations */
  r = reestimate_alloc(mo);
  if (!r) {mes_proc(); goto STOP; };
  
  log_p_old = -DBL_MAX; 
  n = 1;

  /* main loop Baum-Welch-Alg. */
  while (n <= MAX_ITER_BW) { 
    
    if (reestimate_one_step(mo, r, sq->seq_number, sq->seq_len, sq->seq, 
			    &log_p, sq->seq_w) == -1) { 
      char *str = 
	mprintf(NULL, 0, "reestimate_one_step false (%d.step)\n",n); 
      mes_prot(str);
      m_free(str);
      goto STOP;
    }
    
    if (n == 1)
      printf("%8.5f (-log_p input model)\n", -log_p); /* */
    else
      printf("%8.5f (-log_p)\n", -log_p); /* */
      
    if (log_p == +1) {
      printf("Reestimate stoped: No sequence can be build from model mo!\n");
      break;
    }

    diff = log_p - log_p_old;
    /* error in convergence ? */
    if ( diff < -EPS_PREC) {
      char *str = 
	mprintf(NULL, 0, "No convergence: log P < log P-old! (n = %d)\n",n); 
      mes_prot(str);
      m_free(str);
      goto STOP;
    }
    else if ( log_p > EPS_PREC ) {
      char *str = 
	mprintf(NULL, 0, "No convergence: log P > 0! (n = %d)\n",n); 
      mes_prot(str);
      m_free(str);
      goto STOP;
    }
      
    /* stop iterations? */
    if ( diff < fabs((double)EPS_ITER_BW * log_p) ) {
      printf("Convergence after %d steps\n", n); 
      break;
    }
    else {
      /* for next iteration */
      log_p_old = log_p;
      reestimate_init(r, mo); /* sets all fields to zero */
      n++;
    }
  } /* while (n <= MAX_ITER) */ 
  
  /* log_p of reestimated model */
  log_p = 0.0;
  valid = 0;
  for (k = 0; k < sq->seq_number; k++) {
    if (foba_logp(mo, sq->seq[k], sq->seq_len[k], &log_p_k) == -1) 
      { mes_proc(); goto STOP; }
    if (log_p_k != +1) {
      log_p += log_p_k;
      valid = 1;
    }
  }
  if (!valid)
    log_p = +1;
  printf("%8.5f (-log_p optimized model)\n", -log_p);
  
  /* check new parameter for plausibility */
  /* if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  return(0);

STOP:
  reestimate_free(&r, mo->N);
  return(-1);
# undef CUR_PROC
} /* reestimate_baum_welch */
