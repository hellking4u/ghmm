/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/sreestimate.c
  created      : TIME: 17:11:14     DATE: Mon 15. November 1999
  $Id$

__copyright__

*******************************************************************************/

#include <math.h>
#include <float.h>
#include "mprintf.h"
#include "mes.h"
#include "sreestimate.h"
#include "smodel.h"
#include "sfoba.h"
#include "matrix.h"
#include "vector.h"
#include "randvar.h"
#include "gauss_tail.h"
#include "const.h"
#include "root_finder.h"

/* switch: turn off (0)  MESCONTR and MESINFO */
#define MCI 1
/* set control output (MES_WIN, MES_FILE, MES_FILE_WIN, ...) */
#define MESCONTR MES_FILE
/* set info output  (logP, ...) */
#define MESINFO MES_FILE

typedef struct local_store_t {
  int cos;
  double *pi_num;
  double pi_denom;
  double ***a_num;
  double **a_denom;
  double **c_num; 
  double *c_denom;
  double **mue_num;
  double **u_num;
  double **mue_u_denom; /* mue-denom. = u-denom. for sym. normal density*/
  double **sum_gt_otot; /* for truncated normal density */
  double **sum_gt_logb; /* Control according to Q-function */
} local_store_t;

/** needed for normaldensitypos (truncated normal density) */
#define ACC 1E-8
static double C_PHI;
static double CC_PHI;
/***/

static local_store_t *sreestimate_alloc(const smodel *smo);
static int sreestimate_free(local_store_t **r, int N);
static int sreestimate_init(local_store_t *r, const smodel *smo);
static int sreestimate_alloc_matvek(double ***alpha, double ***beta, 
				    double **scale, double ****b, 
				    int T, int N, int M);
static int sreestimate_precompute_b(smodel *smo, double *O, int T, double ***b);
static int sreestimate_free_matvec(double **alpha, double **beta, double *scale,
				   double ***b, int T, int N);
static int sreestimate_setlambda(local_store_t *r, smodel *smo);
static int sreestimate_one_step(smodel *smo, local_store_t *r, 
			       int seq_number, int *T, double **O,
			       double *log_p, double *seq_w);
/*----------------------------------------------------------------------------*/
/* various allocations */
static local_store_t *sreestimate_alloc(const smodel *smo) {
# define CUR_PROC "sreestimate_alloc"
  int i;
  local_store_t* r = NULL;
  if (!m_calloc(r, 1)) {mes_proc(); goto STOP;}
  r->cos=smo->cos;
  if (!m_calloc(r->pi_num, smo->N)) {mes_proc(); goto STOP;}
  if (!m_calloc(r->a_num, smo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < smo->N; i++) {
    r->a_num[i] = matrix_d_alloc(smo->cos, smo->s[i].out_states);
    if (!r->a_num[i]) {mes_proc(); goto STOP;}
  }
  r->a_denom = matrix_d_alloc(smo->N, smo->cos);
  if (!r->a_denom) {mes_proc(); goto STOP;}
  /***/
  if (!m_calloc(r->c_denom, smo->N)) {mes_proc(); goto STOP;}
  r->c_num = matrix_d_alloc(smo->N, smo->M);
  if (!(r->c_num)) {mes_proc(); goto STOP;}
  r->mue_num = matrix_d_alloc(smo->N, smo->M);
  if (!(r->mue_num)) {mes_proc(); goto STOP;}
  r->u_num = matrix_d_alloc(smo->N, smo->M);
  if (!(r->u_num)) {mes_proc(); goto STOP;}
  r->mue_u_denom = matrix_d_alloc(smo->N, smo->M);
  if (!(r->mue_u_denom)) {mes_proc(); goto STOP;}
  r->sum_gt_otot = matrix_d_alloc(smo->N, smo->M);
  if (!(r->sum_gt_otot)) {mes_proc(); goto STOP;}
  r->sum_gt_logb = matrix_d_alloc(smo->N, smo->M);
  if (!(r->sum_gt_logb)) {mes_proc(); goto STOP;}
  return(r);
STOP:
  sreestimate_free(&r, smo->N);
  return 0;
# undef CUR_PROC
} /* sreestimate_alloc */

/*----------------------------------------------------------------------------*/
static int sreestimate_free(local_store_t **r, int N) {
# define CUR_PROC "sreestimate_free"
  int i;
  mes_check_ptr(r, return(-1));
  if( !*r ) return(0);
  m_free((*r)->pi_num);  
  for (i = 0; i < N; i++)
    matrix_d_free( &((*r)->a_num[i]), (*r)->cos);
  m_free((*r)->a_num);
  matrix_d_free( &((*r)->a_denom), N);
  /***/
  m_free((*r)->c_denom);
  matrix_d_free( &((*r)->c_num), N );
  matrix_d_free( &((*r)->mue_num), N );
  matrix_d_free( &((*r)->u_num), N );
  matrix_d_free( &((*r)->mue_u_denom), N );
  matrix_d_free( &((*r)->sum_gt_otot), N );
  matrix_d_free( &((*r)->sum_gt_logb), N );
  m_free(*r);
  return(0);
# undef CUR_PROC
} /* sreestimate_free */

/*----------------------------------------------------------------------------*/
static int sreestimate_init(local_store_t *r, const smodel *smo) {
# define CUR_PROC "sreestimate_init"
  int i, j, m, osc;
  r->pi_denom = 0.0;
  for (i = 0; i < smo->N; i++) {
    r->pi_num[i] = 0.0;
    for (osc = 0; osc < smo->cos; osc++) {
      r->a_denom[i][osc] = 0.0;
      for (j = 0; j < smo->s[i].out_states; j++)
	r->a_num[i][osc][j] = 0.0;
    }
    r->c_denom[i] = 0.0;
    for (m = 0; m < smo->M; m++) {
      r->c_num[i][m] = 0.0;
      r->mue_num[i][m] = 0.0;
      r->u_num[i][m] = 0.0;
      r->mue_u_denom[i][m] = 0.0;
      r->sum_gt_otot[i][m] = 0.0;
      r->sum_gt_logb[i][m] = 0.0;
    }
  }
  return(0);
# undef CUR_PROC
} /* sreestimate_init */

/*----------------------------------------------------------------------------*/
static int sreestimate_alloc_matvek(double ***alpha, double ***beta, 
				    double **scale, double ****b, 
				    int T, int N, int M) {
# define CUR_PROC "sreestimate_alloc_matvek"
  int t, res = -1;
  *alpha = matrix_d_alloc(T, N);
  if (!(*alpha)) {mes_proc(); goto STOP;}
  *beta = matrix_d_alloc(T, N);
  if (!(*beta)) {mes_proc(); goto STOP;}
  if (!m_calloc(*scale, T)) {mes_proc(); goto STOP;}
  /* 3-dim. matrix for b[t][i][m] with m = 1..M(!): */
  if (!m_calloc(*b, T)) {mes_proc(); goto STOP;}
  for (t = 0; t < T; t++) {
    (*b)[t] = matrix_d_alloc(N, M+1);
    if (!((*b)[t])) {mes_proc(); goto STOP;}
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* sreestimate_alloc_matvek */

/*----------------------------------------------------------------------------*/
static int sreestimate_free_matvec(double **alpha,double **beta,double *scale,
				   double ***b, int T, int N) {
# define CUR_PROC "sreestimate_free_matvec"
  int t;
  matrix_d_free(&alpha, T);
  matrix_d_free(&beta, T);
  m_free(scale); 
  if (!b) return(0);
  for (t = 0; t < T; t++)
    matrix_d_free(&b[t], N);
  m_free(b);
  return(0);
# undef CUR_PROC
} /* sreestimate_free_matvec */   

/*----------------------------------------------------------------------------*/
static int sreestimate_precompute_b(smodel *smo, double *O, int T, double ***b){
# define CUR_PROC "sreestimate_precompute_b"
  int t, i, m;
  /* save sum (c_im * b_im(O_t))  in b[t][i][smo->M] */
  for (t = 0; t < T; t++)
    for (i = 0; i < smo->N; i++)
      b[t][i][smo->M] = 0.0;
  /* save c_im * b_im(O_t)  directly in  b[t][i][m] */
  for (t = 0; t < T; t++)
    for (i = 0; i < smo->N; i++)
      for (m = 0; m < smo->M; m++) {
	b[t][i][m] = smodel_calc_cmbm(smo, i, m, O[t]); 
	b[t][i][smo->M] += b[t][i][m];
      }
  return(0);
# undef CUR_PROC
} /* sreestimate_precompute_b */   


/*----------------------------------------------------------------------------*/
static int sreestimate_setlambda(local_store_t *r, smodel *smo) {
# define CUR_PROC "sreestimate_setlambda"
  int res = -1;
  int i, j, m, l, j_id, osc;
  double pi_factor, a_factor_i = 0.0, c_factor_i = 0.0, u_im, mue_im, 
    mue_left, mue_right, A, B, Atil, Btil; /* Q; */
  int a_num_pos, a_denom_pos, c_denom_pos, c_num_pos;
  char *str;

  double p_i;
  int g, h;

  if (r->pi_denom <= DBL_MIN)
    { mes_prot("pi: denominator == 0.0!\n");  goto STOP; } 
  else
    pi_factor = 1/r->pi_denom;

  for (i = 0; i < smo->N; i++) {
    /* Pi */
    smo->s[i].pi =  r->pi_num[i] * pi_factor;
    
    /* A */
    for (osc = 0; osc < smo->cos; osc++) {
      /* note: denom. might be 0; never reached state? */
      a_denom_pos = 1;
      if (r->a_denom[i][osc] <= DBL_MIN) {
	a_denom_pos = 0;
#if MCI
	if (smo->s[i].out_states > 0) { 
	  if (smo->s[i].in_states == 0)
	    mes(MESCONTR,"state %d: no in_states\n",i);
	  else {
	    /* Test: sum all ingoing == 0 ? */
	    p_i = smo->s[i].pi;
	    for (g = 0; g < smo->cos; g++)
	      for (h = 0; h < smo->s[i].in_states; h++)
		p_i += smo->s[i].in_a[g][h];
	    if (p_i == 0)
	      mes(MESCONTR,"state %d: P(in) = 0\n",i);
	    else
	      mes(MESCONTR,"state %d can be reached but a-denom. = 0.0\n",i);
	  }
	}
#endif
      }
      else
	a_factor_i = 1/r->a_denom[i][osc];
      
      a_num_pos = 0;
      for (j = 0; j < smo->s[i].out_states; j++) {
	j_id = smo->s[i].out_id[j];
	/* TEST: denom. < numerator */
	if ((r->a_denom[i][osc] - r->a_num[i][osc][j]) < -EPS_PREC) {
#if MCI
	  mes(MESCONTR,"a[%d][%d][%d]: numerator > denom.!\n", i, osc, j_id); 
#endif
	  smo->s[i].out_a[osc][j] = 1.0;
	}
	else if (a_denom_pos)
	  smo->s[i].out_a[osc][j] = r->a_num[i][osc][j] * a_factor_i;
	else 	  
	  smo->s[i].out_a[osc][j] = 0.0;
	
	if (r->a_num[i][osc][j] > 0.0) /* >= EPS_PREC ? */ 
	  a_num_pos = 1;  
      
	/* important: also update in_a  */
	l = 0;
	while (l < smo->s[j_id].in_states)
	  if  (smo->s[j_id].in_id[l] == i)
	    break;
	  else 
	    l++;
	if ( l == smo->s[j_id].in_states) { 
	  str = mprintf(NULL, 0, "no matching in_a for out_a(=a[%d][%d]) found!\n",i,j_id); 
	  mes_prot(str);
	  m_free(str);
	  goto STOP;
	}
	smo->s[j_id].in_a[osc][l] = smo->s[i].out_a[osc][j]; 
      } /* j-loop */
      
#if MCI
      if (!a_num_pos && smo->s[i].out_states > 0)
	mes(MESCONTR,"all numerators a[%d][%d][j]==0 (denom. = %.4f, P(in)=%.4f)!\n", i, osc, r->a_denom[i][osc], p_i);
#endif
    } /* osc-loop */
    
    
    /* C, Mue und U */

    /* if fix, continue to next state */
    if (smo->s[i].fix)
      continue;

    c_denom_pos = 1;
    if (r->c_denom[i] <= DBL_MIN) { /* < EPS_PREC ? */
#if MCI
      mes(MESCONTR,"c[%d][m]: denominator == 0.0!\n", i);
#endif
      c_denom_pos = 0;
    }
    else 
      c_factor_i = 1/r->c_denom[i];
    
    c_num_pos = 0;
    for (m = 0; m < smo->M; m++) {

      /* TEST: denom. < numerator */
      if ((r->c_denom[i] - r->c_num[i][m]) < 0.0) { /* < -EPS_PREC ? */
#if MCI
	mes(MESCONTR,"c[%d][%d]: numerator > denominator! (%.4f > %.4f)!\n", 
	    i, m, r->c_num[i][m], r->c_denom[i]);
#endif
	smo->s[i].c[m] = 1.0;
      }
      else if (c_denom_pos) 
	/* c_denom == 0: no change in c_im (?) */
	smo->s[i].c[m] = r->c_num[i][m] * c_factor_i;
      
      if (r->c_num[i][m] > 0.0)
	c_num_pos = 1;
      
      if ( fabs(r->mue_u_denom[i][m]) <= DBL_MIN) /* < EPS_PREC ? */
#if MCI
	mes(MESCONTR,"mue[%d][%d]: denominator == 0.0!\n", i, m);
#else
      ;
#endif
      else { 
	/* set mue_im */
	smo->s[i].mue[m] = r->mue_num[i][m] / r->mue_u_denom[i][m];		
      }
      
      /* TEST: u_denom == 0.0 ? */
      if ( fabs(r->mue_u_denom[i][m]) <= DBL_MIN) { /* < EPS_PREC ? */
#if MCI
	mes(MESCONTR,"u[%d][%d]: denominator == 0.0!\n", i, m);
#endif
	;
	/* smo->s[i].u[m]  unchanged! */
      }
      else {        	
	u_im = r->u_num[i][m] / r->mue_u_denom[i][m];		
	if (u_im <= EPS_U)
	  u_im = (double)EPS_U;
	smo->s[i].u[m] = u_im;
      }
      
      /* modification for truncated normal density:
	 1-dim optimization for mue, calculate u directly 
	 note: if denom == 0 --> mue and u not recalculated above */     
      if (smo->density == normal_pos && fabs(r->mue_u_denom[i][m]) > DBL_MIN) {
	A = smo->s[i].mue[m];
	B = r->sum_gt_otot[i][m] / r->mue_u_denom[i][m]; 

	/* A^2 ~ B -> search max at border of EPS_U */
	if (B - A*A < EPS_U) {
	  mue_left = -EPS_NDT; /* attention: only works if  EPS_NDT > EPS_U ! */
	  mue_right = A;

	  if ((pmue_umin(mue_left, A, B, EPS_NDT) > 0.0 && 
	       pmue_umin(mue_right, A, B, EPS_NDT) > 0.0) ||
	      (pmue_umin(mue_left, A, B, EPS_NDT) < 0.0 && 
	       pmue_umin(mue_right, A, B, EPS_NDT) < 0.0))
	    fprintf(stderr,"umin:fl:%.3f\tfr:%.3f\t; left %.3f\t right %3f\t A %.3f\t B %.3f\n",
		    pmue_umin(mue_left, A, B, EPS_NDT),
		    pmue_umin(mue_right, A, B, EPS_NDT),
		    mue_left, mue_right, A, B);
	  
	  mue_im = zbrent_AB(pmue_umin, mue_left, mue_right, 
			     ACC, A, B, EPS_NDT);
	  u_im = EPS_U;
	}
	else {
	  Atil = A + EPS_NDT;
	  Btil = B + EPS_NDT*A;
	  mue_left = (-C_PHI * sqrt(Btil + EPS_NDT*Atil 
				    + CC_PHI*m_sqr(Atil)/4.0) 
		      - CC_PHI*Atil/2.0 - EPS_NDT)*0.99;
	  mue_right = A;
	  if (A < Btil*randvar_normal_density_pos(-EPS_NDT,0,Btil)) 
	    mue_right = m_min(EPS_NDT,mue_right);
	  else 
	    mue_left = m_max(-EPS_NDT, mue_left);
	  mue_im = zbrent_AB(pmue_interpol, 
			     mue_left, mue_right, ACC, A, B, EPS_NDT);	 
	  u_im = Btil - mue_im*Atil;
	}
	/* set modified values of mue and u */
	smo->s[i].mue[m] = mue_im;
	if (u_im < (double)EPS_U)
	  u_im = (double)EPS_U;
	smo->s[i].u[m] = u_im;
      } /* end modifikation truncated density */
      
    } /* for (m ..) */
    
#if MCI
    if (!c_num_pos)
      mes(MESCONTR,"all numerators c[%d][m] == 0 (denominator=%.4f)!\n", 
	  i, r->c_denom[i]);
#endif
    
  } /* for (i = 0 .. < smo->N)  */
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* sreestimate_setlambda */


/*----------------------------------------------------------------------------*/
static int sreestimate_one_step(smodel *smo, local_store_t *r, int seq_number, 
				int *T, double **O, double *log_p, double *seq_w) { 
# define CUR_PROC "sreestimate_one_step"
  int res = -1;  
  int k, i, j, m, t, j_id, valid_parameter, valid_logp, osc;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  double ***b = NULL;
  int T_k = 0, T_k_max = 0, tilgphase = 0;
  double c_t, sum_alpha_a_ji, gamma, gamma_ct, f_im, quot;
  double log_p_k, osum = 0.0;
  
  *log_p = 0.0;
  valid_parameter = valid_logp = 0;
  
  /*scan for max T_k: alloc of alpha, beta, scale and b only once */
  T_k_max = T[0];
  for (k = 1; k < seq_number; k++)
    if (T[k] > T_k_max) T_k_max = T[k];
  if (sreestimate_alloc_matvek(&alpha, &beta, &scale, &b, T_k_max, smo->N, 
			       smo->M) == -1) {mes_proc(); goto STOP;}
  
  /* loop over all sequences */
  for (k = 0; k < seq_number; k++) {    
    /* Test: ignore sequences with very small weights */
    /* if  (seq_w[k] < 0.0001) 
       continue; 
    */
    /* seq. is used for calculation of log_p */
    valid_logp++;
    T_k = T[k]; 
    /* precompute output densities */    
    sreestimate_precompute_b(smo, O[k], T_k, b);
 
    if ((sfoba_forward(smo, O[k], T_k, b, alpha, scale, &log_p_k) == -1) ||
	(sfoba_backward(smo, O[k], T_k, b, beta, scale) == -1)) {
#if MCI
      mes(MESCONTR,"O(%2d) can't be build from smodel smo!\n", k);
#endif
      /* penalty costs */
      *log_p += seq_w[k] * (double)PENALTY_LOGP; 
 
      continue;
    }
    else 
      /* weighted error function */
      *log_p += log_p_k * seq_w[k];  
    /* seq. is used for parameter estimation */
    valid_parameter++;
       
    /* loop over all states*/
    for (i = 0; i < smo->N; i++) {
      
      /* Pi */
      r->pi_num[i] += seq_w[k] * alpha[0][i] * beta[0][i];
      r->pi_denom += seq_w[k] * alpha[0][i] * beta[0][i]; /* sum over all i */
      
      /* loop over t (time steps of seq.)  */
      for (t = 0; t < T_k; t++) {	
	c_t = 1/scale[t]; 
	if (t > 0) {
	  osc = sequence_d_class(O[k], t - 1, &osum); /* dummy */
	  /* A: starts at t=1 !!! */
	  r->a_denom[i][osc] += seq_w[k] * alpha[t-1][i] * beta[t-1][i];
	  for (j = 0; j < smo->s[i].out_states; j++) {
	    j_id = smo->s[i].out_id[j];
	    r->a_num[i][osc][j] += 
	      ( seq_w[k] * alpha[t-1][i]	* smo->s[i].out_a[osc][j] *
		b[t][j_id][smo->M]	* beta[t][j_id] * c_t ); /* c[t] = 1/scale[t] */
	  }
	  /* calculate sum (j=1..N){alp[t-1][j]*a_jc(t-1)i} */
	  sum_alpha_a_ji = 0.0;
	  for (j = 0; j < smo->s[i].in_states; j++) {
	    j_id = smo->s[i].in_id[j];
	    sum_alpha_a_ji += alpha[t-1][j_id] * smo->s[i].in_a[osc][j];
	  }
	}
	else {
	  /* calculate sum(j=1..N){alpha[t-1][j]*a_jci}, which is used below
	     for (t=1) = pi[i] (alpha[-1][i] not defined) !!! */
	  sum_alpha_a_ji = smo->s[i].pi; 
	}
	/*=========falls state fix, continue;======================*/	
	if (smo->s[i].fix)
	  continue;
	/* C-denominator: */
	r->c_denom[i] += seq_w[k] * alpha[t][i] * beta[t][i];  
	
	/* if sum_alpha_a_ji == 0.0, all following values are 0! */
	if ( sum_alpha_a_ji == 0.0)
	  continue; /* next t */
	
	/* loop over no of density functions for C-numer., mue and u */
	for (m = 0; m < smo->M; m++) {	  
	  /*  c_im * b_im  */
	  f_im = b[t][i][m];
	  gamma = seq_w[k] * sum_alpha_a_ji * f_im * beta[t][i];
	  gamma_ct = gamma * c_t; /* c[t] = 1/scale[t] */
	  
	  /* numerator C: */
	  r->c_num[i][m] += gamma_ct; 
	  
	  /* numerator Mue: */
	  r->mue_num[i][m] += (gamma_ct * O[k][t]);
	  /* denom. Mue/U: */
	  r->mue_u_denom[i][m] += gamma_ct;
	  /* numerator U: */
	  r->u_num[i][m] += (gamma_ct * m_sqr(O[k][t]-smo->s[i].mue[m]));
	  
	  /* sum gamma_ct * O[k][t] * O[k][t] (truncated normal density): */
	  r->sum_gt_otot[i][m] += (gamma_ct * m_sqr(O[k][t]));
	  
	  /* sum gamma_ct * log(b_im) */
	  if (gamma_ct > 0.0) {
	    quot = b[t][i][m] / smo->s[i].c[m];
	    r->sum_gt_logb[i][m] += (gamma_ct * log(quot) );
	  }
	}
	
      } /* for (t=0, t<T) */
      
    } /* for (i=0, i<smo->N) */
          
  } /* for (k = 0; k < seq_number; k++) */
  
  sreestimate_free_matvec(alpha, beta, scale, b, T_k_max, smo->N); 

  if (valid_parameter) {
    /* new parameter lambda: set directly in model */
    if ( sreestimate_setlambda(r, smo) == -1 ) { mes_proc(); return(-1); } 
    /* only test : */
    /* if (smodel_check(smo) == -1) { mes_proc(); goto STOP; } */
  }
  else { /* NO sequence can be build from smodel smo! */
    /* diskret:  *log_p = +1; */
    mes(MES_WIN, " NO sequence can be build from smodel smo!\n"); 
    return(-1);
  }
  return(valid_logp);
  /*  return(valid_parameter); */
STOP:

  sreestimate_free_matvec(alpha, beta, scale, b, T_k, smo->N);
  return(res);
# undef CUR_PROC
} /* sreestimate_one_step */


/*============================================================================*/
/* int sreestimate_baum_welch(smodel *smo, sequence_d_t *sqd) {*/
int sreestimate_baum_welch(smosqd_t *cs) {
# define CUR_PROC "sreestimate_baum_welch"
  int n, valid, valid_old, max_iter_bw;
  double log_p, log_p_old, diff, eps_iter_bw;
  local_store_t *r = NULL;
  char *str;
 
  /* truncated normal density needs static varialbles C_PHI and
     CC_PHI */
  if (cs->smo->density == normal_pos) {
    C_PHI = randvar_get_xPHIless1();
    CC_PHI = m_sqr(C_PHI);
  }

  /* local store for all iterations */
  r = sreestimate_alloc(cs->smo);
  if (!r) {mes_proc(); goto STOP; };
  sreestimate_init(r, cs->smo);
  
  log_p_old = -DBL_MAX;
  valid_old = cs->sqd->seq_number;
  n = 1;

  max_iter_bw = m_min(MAX_ITER_BW, cs->max_iter);
  eps_iter_bw = m_max(EPS_ITER_BW, cs->eps);

  while (n <= max_iter_bw) {      
    valid = sreestimate_one_step(cs->smo, r, cs->sqd->seq_number, 
				 cs->sqd->seq_len, cs->sqd->seq, &log_p,
				 cs->sqd->seq_w);
    /* to follow convergence of bw: uncomment next line */
    /*    printf("\tBW Iter %d\t log(p) %.4f\n", n, log_p); */
    if (valid == -1) { 
      str = mprintf(NULL, 0, "sreestimate_one_step false (%d.step)\n",n); 
      mes_prot(str);
      m_free(str);
      goto STOP;
    }
    
#if MCI
    if (n == 1)
      mes(MESINFO, "%8.5f (-log_p input smodel)\n", -log_p);
    else
      mes(MESINFO, "\n%8.5f (-log_p)\n", -log_p);
#endif
      
    diff = log_p - log_p_old;

    if (diff < -EPS_PREC) {
      if (valid > valid_old) {
	str = mprintf(NULL,0,"log P < log P-old (more sequences (%d) , n = %d)\n",
		      valid - valid_old, n);
	mes_prot(str); m_free(str);
      }
      /* no convergence */
      else {
	str = mprintf(NULL,0,"NO convergence: log P(%e) < log P-old(%e)! (n = %d)\n",
		      log_p, log_p_old, n);
	mes_prot(str); m_free(str);
	break;  /* goto STOP; ? */
      }
    }
    
    /* stop iteration */
    if ( diff >= 0.0 && 
	 diff < fabs(eps_iter_bw * log_p) ) {
#if MCI
      mes(MESINFO, "Convergence after %d steps\n", n); 
#endif
      break;
    }
    else {
      /* for next iteration */
      valid_old = valid;
      log_p_old = log_p;
      /* set values to zero */
      sreestimate_init(r, cs->smo);
      n++;
    }
    
  } /* while (n <= MAX_ITER_BW) */ 
  
#if MCI
  mes(MESINFO, "%8.5f (-log_p optimized smodel)\n", -log_p);
#endif
  /* log_p outside this function */
  *cs->logp = log_p;

  /* test plausibility of new parameters */    
  /*  if (smodel_check(mo) == -1) { mes_proc(); goto STOP; } */
  sreestimate_free(&r, cs->smo->N);
  return(0);
STOP:
  sreestimate_free(&r, cs->smo->N);
  return(-1);
# undef CUR_PROC
} /* sreestimate_baum_welch */

#undef ACC
#undef MCI
#undef MESCONTR
#undef MESINFO

