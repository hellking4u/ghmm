/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/ghmm/smodel.c
  created      : TIME: 21:54:32     DATE: Sun 14. November 1999
  $Id$

__copyright__

*******************************************************************************/

#include <float.h>
#include <math.h>
#include "mprintf.h"
#include "mes.h"
#include "smodel.h"
#include "sfoba.h"
#include "sgenerate.h"
#include "sreestimate.h"
#include "matrix.h"
#include "vector.h"
#include "sequence.h"
#include "const.h"
#include "rng.h"
#include "randvar.h"



/*----------------------------------------------------------------------------*/
static int smodel_state_alloc(sstate *state, int M, int in_states,
			      int out_states) {
# define CUR_PROC "smodel_state_alloc"
  int res = -1;
  if(!m_calloc(state->c, M)) {mes_proc(); goto STOP;}
  if(!m_calloc(state->mue, M)) {mes_proc(); goto STOP;}
  if(!m_calloc(state->u, M)) {mes_proc(); goto STOP;}
  if (out_states > 0) {
    if(!m_calloc(state->out_id, out_states)) {mes_proc(); goto STOP;}
    /* if COS > 1: out_a is a matrix */
    state->out_a = matrix_d_alloc(COS, out_states);
    if(!state->out_a) {mes_proc(); goto STOP;}
  }
  if (in_states > 0) {
    if(!m_calloc(state->in_id, in_states)) {mes_proc(); goto STOP;}
    /* if COS > 1: in_a is a matrix */
    state->in_a = matrix_d_alloc(COS, in_states);
    if(!state->in_a) {mes_proc(); goto STOP;}
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* smodel_state_alloc */


/*----------------------------------------------------------------------------*/
static int smodel_copy_vectors(smodel *smo, int index, double *pi, int *fix,
			       double ***a_matrix, double **c_matrix, 
			       double **mue_matrix, double **u_matrix) {
#define CUR_PROC "smodel_alloc_vectors"
  int i, c, exist, count_out = 0, count_in = 0;
  smo->s[index].pi = pi[index];
  smo->s[index].fix = fix[index];
  for (i = 0; i < smo->M; i++){
    smo->s[index].c[i] = c_matrix[index][i];
    smo->s[index].mue[i] = mue_matrix[index][i];
    smo->s[index].u[i] = u_matrix[index][i];
  }

  for (i = 0; i < smo->N; i++) {
    exist = 0;
    for (c = 0; c < COS; c++) {
      if (a_matrix[c][index][i]) { 
	exist = 1;
	break;
      }
    }
    /* transition to successor possible at least in one transition class */
    if (exist) { 
      if (count_out >= smo->s[index].out_states) {mes_proc(); return(-1);}
      smo->s[index].out_id[count_out] = i;
      for (c = 0; c < COS; c++) 
	smo->s[index].out_a[c][count_out] = a_matrix[c][index][i];
      count_out++;
    }
    exist = 0;
    for (c = 0; c < COS; c++) {
      if (a_matrix[c][i][index]) { 
	exist = 1;
	break;
      }
    }
    /* transition to predecessor possible at least in one transition class */
    if (exist) {
      if (count_in >= smo->s[index].in_states) {mes_proc(); return(-1);}
      smo->s[index].in_id[count_in] = i;
      for (c = 0; c < COS; c++)
	smo->s[index].in_a[c][count_in] = a_matrix[c][i][index];
      count_in++;
    }
  }
  return(0);
#undef CUR_PROC
} /* smodel_alloc_vectors */


/*============================================================================*/
smodel **smodel_read(char *filename, int *smo_number) {
#define CUR_PROC "smodel_read"
  int j;
  long new_models = 0;
  scanner_t *s = NULL;
  smodel **smo = NULL;
  *smo_number = 0;
  s = scanner_alloc(filename);  if(!s) {mes_proc(); goto STOP;}
  while(!s->err && !s->eof) {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP;
    if (!strcmp(s->id, "SHMM") || !strcmp(s->id, "shmm")) {
      (*smo_number)++;
      /* more mem */	 
      if (m_realloc(smo, *smo_number)) { mes_proc(); goto STOP; }
      smo[*smo_number - 1] = smodel_read_block(s, (int *) &new_models); 
      if (!smo[*smo_number - 1]) { mes_proc(); goto STOP; }
      /* copy smodel */
      if (new_models > 1) { 
	/* "-1" due to  (*smo_number)++ from above  */
	if (m_realloc(smo, *smo_number - 1 + new_models)) { 
	  mes_proc(); goto STOP; 
	}
	for (j = 1; j < new_models; j++) {
	  smo[*smo_number] = smodel_copy(smo[*smo_number - 1]);
	  if (!smo[*smo_number]) { mes_proc(); goto STOP; }
	  (*smo_number)++;
	}
      }
    }  
    else {
      scanner_error(s, "unknown identifier");
      goto STOP;
    }
    scanner_consume(s, ';'); if(s->err) goto STOP; 
  } /* while(!s->err && !s->eof) */
  if (smodel_check_compatibility(smo, *smo_number) == -1) {
    mes_proc(); goto STOP; 
  }
  return smo;
STOP:
  return NULL;
#undef CUR_PROC
} /* smodel_read */


/*============================================================================*/
smodel *smodel_read_block(scanner_t *s, int *multip){
#define CUR_PROC "smodel_read_block"
  int i,j,osc, m_read, n_read, pi_read, a_read, c_read, mue_read, 
    u_read, len, density_read, out, in, prior_read, fix_read;
  smodel *smo      = NULL;
  double *pi_vektor = NULL, **a_matrix = NULL, ***a_3dmatrix = NULL;
  double **c_matrix = NULL, **mue_matrix = NULL, **u_matrix = NULL;
  int *fix_vektor = NULL;
  
  m_read = n_read = pi_read = a_read = c_read = mue_read = u_read 
    = density_read = prior_read = fix_read = 0;
  *multip = 1; /* default */
  if (!(m_calloc(smo, 1))) { mes_proc(); goto STOP; }
  scanner_consume( s, '{' ); if(s->err) goto STOP; 
  while(!s->err && !s->eof && s->c - '}') {    
    scanner_get_name(s);
    if (strcmp(s->id, "M") && strcmp(s->id, "N") && strcmp(s->id, "Pi") &&
	strcmp(s->id, "A") && strcmp(s->id, "C") && strcmp(s->id, "Mue") &&
	strcmp(s->id, "U") && strcmp(s->id, "multip") && 
	strcmp(s->id, "prior") && strcmp(s->id, "fix_state") &&
	strcmp(s->id, "density") && strncmp(s->id, "Ak_", 3) ) {
      scanner_error(s, "unknown identifier"); goto STOP;
    }    
    scanner_consume(s, '='); if(s->err) goto STOP;
    if (!strcmp(s->id, "multip")) {
      *multip = scanner_get_int(s);
      if (*multip < 1) {/* ignore: make no sense */
	*multip = 1;
	mes_prot("Multip < 1 ignored\n");
      }
    }
    else if (!strcmp(s->id, "M")) {/* number of output components */
      if (m_read) {scanner_error(s, "identifier M twice"); goto STOP;}
      smo->M = scanner_get_int(s);
      m_read = 1;
    }
    else if (!strcmp(s->id, "N")) {/* number of states */
      if (n_read) {scanner_error(s, "identifier N twice"); goto STOP;}
      smo->N = scanner_get_int(s);
      if (!m_calloc(smo->s, smo->N)) {mes_proc(); goto STOP;}
      n_read = 1;
    }   
    else if (!strcmp(s->id, "density")) {/* which density function? */
      if (density_read) {scanner_error(s,"identifier density twice");goto STOP;}
      smo->density = (density_t)scanner_get_int(s);
      if ((int)smo->density < 0 || smo->density >= density_number)
	{ scanner_error(s, "unknown typ of density function"); goto STOP; }
      density_read = 1;
    }   
    else if (!strcmp(s->id, "prior")) {/* modelprior */
      if (prior_read) {scanner_error(s,"identifier prior twice");goto STOP;}
      smo->prior = scanner_get_double(s);
      if ((smo->prior < 0 || smo->prior > 1) && smo->prior != -1)
	{ scanner_error(s, "invalid model prior"); goto STOP; }
      prior_read = 1;
    }   
    else if (!strcmp(s->id, "Pi")) {/* Initial State Prob. */
      if (!n_read) {scanner_error(s, "need N as a range for Pi"); goto STOP;}
      if (pi_read) {scanner_error(s, "identifier Pi twice"); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "vector")) {
	scanner_consume( s, '{' ); if(s->err) goto STOP; 
	pi_vektor = scanner_get_double_array(s, &len);
	if (len != smo->N) 
	  {scanner_error(s, "wrong number of elements in PI"); goto STOP;}
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
	fix_vektor = scanner_get_int_array(s, &len);
	if (len != smo->N) 
	  {scanner_error(s, "wrong number of elements in fix_state"); goto STOP;}
	scanner_consume( s, ';' ); if(s->err) goto STOP;
	scanner_consume( s, '}' ); if(s->err) goto STOP; 
	fix_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    /* 1. Possibility: one matrix for all transition classes */
    else if (!strcmp(s->id, "A")) {
      if (!n_read) {scanner_error(s, "need N as a range for A"); goto STOP;}
      if (a_read) {scanner_error(s, "Identifier A twice"); goto STOP;}
      if (!m_calloc(a_matrix, smo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if (matrix_d_read(s, a_matrix, smo->N, smo->N)){
	  scanner_error(s, "unable to read matrix A"); goto STOP;
	}
	a_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
      /* copy transition matrix to all transition classes */
      if (!m_calloc(a_3dmatrix, COS)) {mes_proc(); goto STOP;}
      a_3dmatrix[0] = a_matrix;
      for (i = 1; i < COS; i++) {
	a_3dmatrix[i] = matrix_d_alloc_copy(smo->N, smo->N, a_matrix);
	if (!a_3dmatrix[i]) {mes_proc(); goto STOP;}
      }
    }
    /* 2. Possibility: one matrix for each transition class specified */
    else if (!strncmp(s->id, "Ak_", 3)) {
      if (!n_read) {scanner_error(s, "need N as a range for A"); goto STOP;}
      if (a_read) {scanner_error(s, "identifier A twice"); goto STOP;}
      if (!m_calloc(a_3dmatrix, COS)) {mes_proc(); goto STOP;}
      for (osc = 0; osc < COS; osc++) {
	if (!m_calloc(a_3dmatrix[osc], smo->N)) {mes_proc(); goto STOP;}
	scanner_get_name(s);
	if (!strcmp(s->id, "matrix")) {
	  if (matrix_d_read(s, a_3dmatrix[osc], smo->N, smo->N)){
	    scanner_error(s, "unable to read matrix Ak"); goto STOP;
	  }
	}
	else {
	  scanner_error(s, "unknown identifier"); goto STOP;
	}
	if (osc < COS-1) {
	  scanner_consume( s, ';' ); if(s->err) goto STOP;
	  /* read next matrix */
	  scanner_get_name(s);
	  if (strncmp(s->id, "Ak_", 3)) {
	    scanner_error(s, "to less matrices Ak"); goto STOP;
	  }
	  scanner_consume(s, '='); if(s->err) goto STOP;
	}
      }
      a_read = 1;
    }
    else if (!strcmp(s->id, "C")) {/* weight for output components */
      if ((!n_read)||(!m_read)){
	scanner_error(s, "need M and N as a range for C"); goto STOP;
      }
      if (c_read) {scanner_error(s, "identifier C twice"); goto STOP;}
      if (!m_calloc(c_matrix, smo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if(matrix_d_read(s, c_matrix, smo->N, smo->M)) {
	  scanner_error(s, "unable to read matrix C"); goto STOP;
	}
	c_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else if (!strcmp(s->id, "Mue")) {/* mean of normal density */
      if ((!n_read)||(!m_read)){
	scanner_error(s, "need M and N as a range for Mue"); goto STOP;
      }
      if (mue_read) {scanner_error(s, "identifier Mue twice"); goto STOP;}
      if (!m_calloc(mue_matrix, smo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if(matrix_d_read(s, mue_matrix, smo->N, smo->M)) {
	  scanner_error(s, "unable to read matrix Mue"); goto STOP;
	}
	mue_read = 1;
      }
      else {
	scanner_error(s, "unknown identifier"); goto STOP;
      }
    }
    else if (!strcmp(s->id, "U")) {/* variances of normal density */
      if ((!n_read)||(!m_read)){
	scanner_error(s, "need M and N as a range for U"); goto STOP;
      }
      if (u_read) {scanner_error(s, "identifier U twice"); goto STOP;}
      if (!m_calloc(u_matrix, smo->N)) {mes_proc(); goto STOP;}
      scanner_get_name(s);
      if (!strcmp(s->id, "matrix")) {
	if(matrix_d_read(s, u_matrix, smo->N, smo->M)) {
	  scanner_error(s, "unable to read matrix U"); goto STOP;
	}
	u_read = 1;
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

  if (m_read == 0 || n_read == 0 || pi_read == 0 || a_read == 0 ||
      c_read == 0 || mue_read == 0 || u_read == 0 || density_read == 0) {
    scanner_error(s, "some identifier not specified (N, M, Pi, A, C, Mue, U or density)"); 
    goto STOP;
  }
  /* set prior to -1 it none was specified */
  if (prior_read == 0)
    smo->prior = -1;
  /* default for fix is 0 */
  if (fix_read == 0) {
    if (!m_calloc(fix_vektor, smo->N)) {mes_proc(); goto STOP;}
    for (i = 0; i < smo->N; i++)
      fix_vektor[i] = 0;
  }
  /* memory alloc for all transition matrices. If a transition is possible in one
     class --> alloc memory for all classes */
  for (i = 0; i < smo->N; i++) {
    for (j = 0; j < smo->N; j++) {
      out = in = 0;
      for (osc = 0; osc < COS; osc++) {
	if (a_3dmatrix[osc][i][j] > 0.0)
	  out = 1;
	if (a_3dmatrix[osc][j][i] > 0.0)
	  in = 1;	
      }
      smo->s[i].out_states += out;
      smo->s[i].in_states += in;
    }
    if (smodel_state_alloc(smo->s + i, smo->M, smo->s[i].in_states,
			   smo->s[i].out_states)) { mes_proc(); goto STOP; }
    /* copy values read to smodel */
    if(smodel_copy_vectors(smo, i, pi_vektor, fix_vektor, a_3dmatrix, c_matrix, mue_matrix,
			   u_matrix)) {mes_proc(); goto STOP;}
  }
  if (a_3dmatrix)
    for (i = 0; i < COS; i++) matrix_d_free(&(a_3dmatrix[i]), smo->N);
  m_free(a_3dmatrix);
  matrix_d_free(&c_matrix, smo->N);
  matrix_d_free(&mue_matrix, smo->N);
  matrix_d_free(&u_matrix, smo->N);
  m_free(pi_vektor); m_free(fix_vektor);
  return(smo);
STOP:
  if (a_3dmatrix) 
    for (i = 0; i < COS; i++) matrix_d_free(&(a_3dmatrix[i]), smo->N);
  m_free(a_3dmatrix);
  matrix_d_free(&c_matrix, smo->N);
  matrix_d_free(&mue_matrix, smo->N);
  matrix_d_free(&u_matrix, smo->N);
  m_free(pi_vektor); m_free(fix_vektor);
  smodel_free(&smo);
  return NULL;
#undef CUR_PROC
} /* smodel_read_block */


/*============================================================================*/
int smodel_free(smodel **smo) {
#define CUR_PROC "smodel_free"
  int i;
  mes_check_ptr(smo, return(-1));
  if( !*smo ) return(0);
  for (i = 0; i < (*smo)->N; i++) {
    m_free((*smo)->s[i].out_id);
    m_free((*smo)->s[i].in_id);
    matrix_d_free(&((*smo)->s[i].out_a), COS);
    matrix_d_free(&((*smo)->s[i].in_a), COS);
    m_free((*smo)->s[i].c);
    m_free((*smo)->s[i].mue);
    m_free((*smo)->s[i].u);
  }
  m_free((*smo)->s);
  m_free(*smo);
  return(0);
#undef CUR_PROC
} /* smodel_free */   


/*============================================================================*/
smodel *smodel_copy(const smodel *smo) {
# define CUR_PROC "smodel_copy"
  int i, k, j, nachf, vorg, m;
  smodel *sm2 = NULL;
  if(!m_calloc(sm2, 1)) {mes_proc(); goto STOP;}
  if (!m_calloc(sm2->s, smo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < smo->N; i++) {
    nachf = smo->s[i].out_states;
    vorg = smo->s[i].in_states;
    if(!m_calloc(sm2->s[i].out_id, nachf)) {mes_proc(); goto STOP;}
    sm2->s[i].out_a = matrix_d_alloc(COS, nachf);
    if(!sm2->s[i].out_a) {mes_proc(); goto STOP;}
    if(!m_calloc(sm2->s[i].in_id, vorg)) {mes_proc(); goto STOP;}
    sm2->s[i].in_a = matrix_d_alloc(COS, vorg);
    if(!sm2->s[i].in_a) {mes_proc(); goto STOP;}
    if(!m_calloc(sm2->s[i].c, smo->M)) {mes_proc(); goto STOP;}
    if(!m_calloc(sm2->s[i].mue, smo->M)) {mes_proc(); goto STOP;}
    if(!m_calloc(sm2->s[i].u, smo->M)) {mes_proc(); goto STOP;}
    /* copy values */     
    for (j = 0; j < nachf; j++) {
      for (k = 0; k < COS; k++)
	sm2->s[i].out_a[k][j] = smo->s[i].out_a[k][j];
      sm2->s[i].out_id[j] = smo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      for (k = 0; k < COS; k++)
	sm2->s[i].in_a[k][j] = smo->s[i].in_a[k][j];
      sm2->s[i].in_id[j] = smo->s[i].in_id[j];
    }
    for (m = 0; m < smo->M; m++) {
      sm2->s[i].c[m] = smo->s[i].c[m];
      sm2->s[i].mue[m] = smo->s[i].mue[m];
      sm2->s[i].u[m] = smo->s[i].u[m];
    }
    sm2->s[i].pi = smo->s[i].pi;
    sm2->s[i].fix = smo->s[i].fix;
    sm2->s[i].out_states = nachf;
    sm2->s[i].in_states = vorg;
  }
  sm2->N = smo->N;
  sm2->M = smo->M;
  sm2->density = smo->density;
  sm2->prior = smo->prior;
  return(sm2);
STOP:
  smodel_free(&sm2);
  return(NULL);
# undef CUR_PROC
} /* smodel_copy */


/*============================================================================*/
int smodel_check(const smodel* smo) {
# define CUR_PROC "smodel_check"
  int res = -1;
  double sum;
  int i,k,j;
  /* sum  Pi[i] == 1 ? */
  sum = 0.0;
  for (i = 0; i < smo->N; i++) {
    sum += smo->s[i].pi;
  }
  if ( fabs(sum - 1.0) >= EPS_PREC )
    { mes_prot("sum Pi[i] != 1.0\n"); goto STOP; }
  /* only 0/1 in fix? */
  for (i = 0; i < smo->N; i++) {
    if (smo->s[i].fix != 0 && smo->s[i].fix != 1) {
      mes_prot("in vector fix_state only 0/1 values!\n"); goto STOP;
    }
  }
  for (i = 0; i < smo->N; i++) {
    if (smo->s[i].out_states == 0) {
      char *str = 
	mprintf(NULL,0,"out_states = 0 (state %d -> final state!)\n",i); 
      mes_prot(str);
    }
    /* sum  a[i][k][j] */
    for (k = 0; k < COS; k++) {
      sum = 0.0;
      for (j = 0; j < smo->s[i].out_states; j++) {
	sum += smo->s[i].out_a[k][j];
      }
      if ( fabs(sum - 1.0) >= EPS_PREC ) { 
	char *str = 
	  mprintf(NULL, 0, "sum out_a[j] = %.2f != 1.0 (state %d)\n", sum, i); 
	mes_prot(str);
	m_free(str);
	/* goto STOP; */
      }
    }
    /* sum c[j] */
    sum = 0.0;
    for (j = 0; j < smo->M; j++)
      sum += smo->s[i].c[j];
    if ( fabs(sum - 1.0) >= EPS_PREC ) { 
      char *str = mprintf(NULL, 0, "sum c[j] = %.2f != 1.0 (state %d)\n",sum,i);
      mes_prot(str);
      m_free(str);
      goto STOP;
    }
    /* check mue, u ? */
  }
  res = 0;
STOP:
  return(res);
# undef CUR_PROC
} /* smodel_check */


/*============================================================================*/
int smodel_check_compatibility(smodel **smo, int smodel_number) {
#define CUR_PROC "smodel_check_compatibility"
  int i, j;
  for (i = 0; i < smodel_number; i++) 
    for (j = i+1; j < smodel_number; j++) {
      if (smo[i]->N != smo[j]->N) {
	char *str = mprintf(NULL, 0, "ERROR: different number of states in smodel %d (%d) and smodel %d (%d)", i, smo[i]->N, j, smo[j]->N); 
	mes_prot(str);
	m_free(str);
	return (-1);
      }
      if (smo[i]->M != smo[j]->M) {
	char *str = mprintf(NULL, 0, "ERROR: different number of possible outputs in smodel  %d (%d) and smodel %d (%d)", i, smo[i]->M, j, smo[j]->M); 
	mes_prot(str);
	m_free(str);
	return (-1);
      }
    }
  return 0;
#undef CUR_PROC    
} /* smodel_check_compatibility */


/*============================================================================*/
double smodel_get_random_var(smodel *smo, int state, int m) {
# define CUR_PROC "smodel_get_random_var"
  switch (smo->density) {
  case normal_approx:
  case normal: 
    return(randvar_normal(smo->s[state].mue[m], smo->s[state].u[m], 0));
  case normal_pos: 
    return(randvar_normal_pos(smo->s[state].mue[m],smo->s[state].u[m],0));
  default: 
    mes(MES_WIN, "Warning: density function not specified!\n"); 
    return(-1);
  }
# undef CUR_PROC
} /* smodel_get_random_var */
  

/*============================================================================*/
int smodel_likelihood(smodel *smo, sequence_d_t *sqd, double *log_p) {
# define CUR_PROC "smodel_likelihood"
  int res = -1;
  double log_p_i;
  int matched, i;

  matched = 0;
  *log_p = 0.0;
  for (i = 0; i < sqd->seq_number; i++) {
    if (sfoba_logp(smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) { 
      *log_p += log_p_i * sqd->seq_w[i];
      matched++;
    }
    else  {
      /* Test: high costs for each unmatched Seq. */
      *log_p += PENALTY_LOGP * sqd->seq_w[i];
      matched++;
      /*      mes(MES_WIN, "sequence[%d] can't be build.\n", i); */
    }
  }
  if (!matched)
    { mes_prot("NO sequence can be build.\n"); goto STOP; }
  /* return number of "matched" sequences */
  res = matched;
STOP:
  return(res);
# undef CUR_PROC
} /* smodel_likelihood */


/*============================================================================*/
/* various print functions
/*============================================================================*/


/*============================================================================*/
void smodel_Ak_print(FILE *file, smodel *smo, int k, char *tab, 
		     char *separator, char *ending) {
  int i, j, out_state;
  for (i = 0; i < smo->N; i++) {
    out_state = 0;
    fprintf(file, "%s", tab);
    if (smo->s[i].out_states > 0 && 
	smo->s[i].out_id[out_state] == 0) {
      fprintf(file, "%.4f", smo->s[i].out_a[k][out_state]);
      out_state++;
    }
    else fprintf(file, "0.0   ");
    for (j = 1; j < smo->N; j++) {
      if (smo->s[i].out_states > out_state && 
	  smo->s[i].out_id[out_state] == j) {
	fprintf(file, "%s %.4f", separator, smo->s[i].out_a[k][out_state]);
	out_state++;
      }
      else fprintf(file, "%s 0.0   ", separator);
    }
    fprintf(file, "%s\n", ending);
  }
} /* smodel_Ak_print */


/*============================================================================*/
void smodel_C_print(FILE *file, smodel *smo, char *tab, char *separator, 
		    char *ending) {
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf(file, "%s", tab);
    fprintf(file, "%.4f", smo->s[i].c[0]);
    for (j = 1; j < smo->M; j++)
      fprintf(file, "%s %.4f", separator, smo->s[i].c[j]);
    fprintf(file, "%s\n", ending);
  }
} /* smodel_C_print */


/*============================================================================*/
void smodel_Mue_print(FILE *file, smodel *smo, char *tab, char *separator, 
		      char *ending) {
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf(file, "%s", tab);
    fprintf(file, "%.4f", smo->s[i].mue[0]);
    for (j = 1; j < smo->M; j++)
      fprintf(file, "%s %.4f", separator, smo->s[i].mue[j]);
    fprintf(file, "%s\n", ending);
  }
} /* smodel_Mue_print */


/*============================================================================*/
void smodel_U_print(FILE *file, smodel *smo, char *tab, char *separator, 
			char *ending) {
  /* attention: choose precision big enough to allow printing of  
     EPS_U in const.h */
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf(file, "%s", tab);
    fprintf(file, "%.4f", smo->s[i].u[0]);
    for (j = 1; j < smo->M; j++)
      fprintf(file, "%s %.4f", separator, smo->s[i].u[j]);
    fprintf(file, "%s\n", ending);
  }
} /* smodel_U_print */


/*============================================================================*/
void smodel_Pi_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending) {
  int i;
  fprintf(file, "%s%.4f", tab, smo->s[0].pi);
  for (i = 1; i < smo->N; i++)
    fprintf(file, "%s %.4f", separator, smo->s[i].pi);
  fprintf(file, "%s\n", ending);
} /* smodel_Pi_print */

/*============================================================================*/
void smodel_fix_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending) {
  int i;
  fprintf(file, "%s%d", tab, smo->s[0].fix);
  for (i = 1; i < smo->N; i++)
    fprintf(file, "%s %d", separator, smo->s[i].fix);
  fprintf(file, "%s\n", ending);
} 


/*============================================================================*/
void smodel_print(FILE *file, smodel *smo) {
  int k;
  fprintf(file, "SHMM = {\n\tM = %d;\n\tN = %d;\n\tdensity = %d;\n", 
	  smo->M, smo->N, (int)smo->density);
  fprintf(file, "\tprior = %.5f;\n", smo->prior);
  fprintf(file, "\tPi = vector {\n");
  smodel_Pi_print(file, smo, "\t", ",", ";");
  fprintf(file, "\t};\n");
  fprintf(file, "\tfix_state = vector {\n");
  smodel_fix_print(file, smo, "\t", ",", ";");
  fprintf(file, "\t};\n");
  for (k = 0; k < COS; k++) {
    fprintf(file, "\tAk_%d = matrix {\n", k);
    smodel_Ak_print(file, smo, k, "\t", ",", ";");
    fprintf(file, "\t};\n");
  }
  fprintf(file, "\tC = matrix {\n");
  smodel_C_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n\tMue = matrix {\n");
  smodel_Mue_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n\tU = matrix {\n");
  smodel_U_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n");
  fprintf(file, "};\n\n");
} /* smodel_print */

/*============================================================================*/
/* needed for hmm_input: only one A (=Ak_1 = Ak_2...) is written */
void smodel_print_oneA(FILE *file, smodel *smo) {
  fprintf(file, "SHMM = {\n\tM = %d;\n\tN = %d;\n\tdensity = %d;\n", 
	  smo->M, smo->N, (int)smo->density);
  fprintf(file, "\tprior = %.3f;\n", smo->prior);
  fprintf(file, "\tPi = vector {\n");
  smodel_Pi_print(file, smo, "\t", ",", ";");
  fprintf(file, "\t};\n");
  fprintf(file, "\tfix_state = vector {\n");
  smodel_fix_print(file, smo, "\t", ",", ";");
  fprintf(file, "\t};\n");
  /* Attention: assumption is: A_k are all the same */
  fprintf(file, "\tA = matrix {\n");
  smodel_Ak_print(file, smo, 0, "\t", ",", ";");
  fprintf(file, "\t};\n");
  /***/
  fprintf(file, "\tC = matrix {\n");
  smodel_C_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n\tMue = matrix {\n");
  smodel_Mue_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n\tU = matrix {\n");
  smodel_U_print(file, smo,"\t", ",", ";");
  fprintf(file, "\t};\n");
  fprintf(file, "};\n\n");
} /* smodel_print */


/*============================================================================*/
double smodel_calc_cmbm(smodel *smo, int state, int m, double omega) {
  double bm = 0.0;
  switch (smo->density) {
  case normal: bm = randvar_normal_density(omega, smo->s[state].mue[m], 
					   smo->s[state].u[m]);
    break;
  case normal_pos: bm = randvar_normal_density_pos(omega,smo->s[state].mue[m],
						   smo->s[state].u[m]);
    break;
  case normal_approx: bm = randvar_normal_density_approx(omega, 
							 smo->s[state].mue[m],
							 smo->s[state].u[m]);
    break;
  default: mes(MES_WIN, "Warning: density function not specified!\n"); 
  }
  if (bm == -1) {
    mes(MES_WIN, "Warning: density function returns -1!\n"); 
    bm = 0.0;
  }
  return(smo->s[state].c[m] * bm);  
} /* smodel_calc_cmbm */


/*============================================================================*/
/* PDF(omega) in a given state */
double smodel_calc_b(smodel *smo, int state, double omega) {
  int m;
  double b = 0.0;
  for (m = 0; m < smo->M; m++)
    b += smodel_calc_cmbm(smo, state, m, omega);
  return(b);
} /* smodel_calc_b */


/*============================================================================*/
double smodel_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, 
			    int verbose)
{
#define CUR_PROC "smodel_prob_distance"

#define STEPS 100

  double p0, p;
  double d = 0.0;
  double *d1;
  sequence_d_t *seq0 = NULL;
  sequence_d_t *tmp = NULL;
  smodel *smo1, *smo2;
  int i, t, a, k, n;
  int true_len;
  int true_number;
  int left_to_right = 0;
  long total, index;
  int step_width = 0;
  int steps = 1;

  if (verbose) {/* If we are doing it verbosely we want to have STEPS steps */ 
    step_width = maxT / STEPS;
    steps = STEPS;
  }
  else        /* else just one */ 
    step_width = maxT;

  if( !m_calloc(d1, steps) ) {mes_proc();goto STOP;} 

  smo1 = cm0;
  smo2 = cm;
 
  for (k = 0; k < 2; k++) {
    
    seq0 = sgenerate_sequences(smo1, 1998, maxT+1, 1, 0, 0);
    
    if (seq0->seq_len[0] < maxT) { /* There is an absorbing state */
            
      /* For now check that Pi puts all weight on state */
      /*
	t = 0;
	for (i = 0; i < smo1->N; i++) {
	if (smo1->s[i].pi > 0.001)
	t++;
	}    
	if (t > 1) {
	mes_prot("No proper left-to-right model. Multiple start states");
	goto STOP;
	} */
      
      left_to_right = 1;
      total = seq0->seq_len[0];
      
      while (total <= maxT) {
	
	/* create a additional sequences at once */
	a = (maxT - total) / (total / seq0->seq_number) + 1;
	/* printf("total=%d generating %d", total, a); */
	tmp = sgenerate_sequences(smo1, 0, 0, a, 0, 0);
	sequence_d_add(seq0, tmp);
	sequence_d_free(&tmp);  

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
	seq0->seq_number = index + 1;
	
	if (smodel_likelihood(smo1, seq0, &p0) == -1) {
	  /* error! */
	  mes_prot("seq0 can't be build from smo1!");
	  goto STOP;
	}
	n = smodel_likelihood(smo2, seq0, &p); // == -1: KEINE Seq. erzeugbar
	if (n < seq0->seq_number) {
	  mes(MES_WIN,
	      "problem: some seqences in seq0 can't be build from smo2\n");
	  /* what shall we do now? */
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
 
    else { /* no left to right model */
      
      true_len = seq0->seq_len[0];
      
      for (t=step_width, i=0; t <= maxT; t+= step_width, i++) {
	seq0->seq_len[0] = t;
	
	if (smodel_likelihood(smo1, seq0, &p0) == -1) {
	  /* error case */
	  mes_prot("seq0 can't be build from smo1!");
	  goto STOP;
	}
	n = smodel_likelihood(smo2, seq0, &p); // == -1: KEINE Seq. erzeugbar
	if (n < seq0->seq_number) {
	  mes(MES_WIN,
	      "problem: some sequences in seq0 can't be build from smo2\n");
	  /* what shall we do now? */
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

      }
      seq0->seq_len[0] = true_len;
    }
    
    if (symmetric) {
      sequence_d_free(&seq0);
      smo1 = cm ;
      smo2 = cm0;
    }
    else
      break;
   
  } /* k = 1,2 */ 

  sequence_d_free(&seq0);
  
  return d;

STOP:
  sequence_d_free(&seq0);
  return(-1.0);
#undef CUR_PROC
}


/*============================================================================*/
double smodel_calc_cmBm(smodel *smo, int state, int m, double omega) {
  double Bm = 0.0;
  switch (smo->density) {
  case normal_approx: 
  case normal:
    Bm = randvar_normal_cdf(omega, smo->s[state].mue[m], smo->s[state].u[m]);
    break;
  case normal_pos: 
    Bm = randvar_normal_pos_cdf(omega,smo->s[state].mue[m],smo->s[state].u[m]);
    break;
  default: mes(MES_WIN, "Warning: density function not specified!\n"); 
  }
  if (Bm == -1) {
    mes(MES_WIN, "Warning: density function returns -1!\n"); 
    Bm = 0.0;
  }
  return(smo->s[state].c[m] * Bm);  
} /* smodel_calc_Bm */


/*============================================================================*/
/* CDF(omega) in a given state */
double smodel_calc_B(smodel *smo, int state, double omega) {
  int m;
  double B = 0.0;
  for (m = 0; m < smo->M; m++)
    B += smodel_calc_cmBm(smo, state, m, omega);
  return(B);
} /* smodel_calc_B */

/*============================================================================*/
/* What is the dimension of the modell ( = dimension of the parameter vektor) ?
   count the number of free parameters in a field of models; used for calc. BIC
   Only those parameters, that can be changed during  training.
   mixture coeff from smix and priors are not counted! 
*/
int smodel_count_free_parameter(smodel **smo, int smo_number) {
  int i, k;
  int pi_counted = 0, cnt = 0;

  for (k = 0; k < smo_number; k++) {
    pi_counted = 0;
    // for states 
    for (i = 0; i < smo[k]->N; i++) {
      if (smo[k]->s[i].out_states > 1) 
	// multipl. with COS correct ???
	cnt += COS * (smo[k]->s[i].out_states - 1);
      if (smo[k]->s[i].pi != 0 && smo[k]->s[i].pi != 1) {
	pi_counted = 1;
	cnt++;
      }
      if (!smo[k]->s[i].fix) {
	if (smo[k]->M == 1) cnt += 2; // mu, sigma
	else cnt += (3 * smo[k]->M); // c, mu, sigma
      }
    } // for (i ..)
    if (pi_counted) cnt--; /* due to condition: sum(pi) = 1 */
    if (smo[k]->M > 1) cnt--; /* due to condition: sum(c) = 1 */
  }

  return cnt;
}

/*============================================================================*/
/* interval (a,b) with ~ B(a) < 0.01, B(b) > 0.99 */  
void smodel_get_interval_B(smodel *smo, int state, double *a, double *b) {
  int m;
  double mue, delta;
  switch (smo->density) {
  case normal:
  case normal_approx: 
  case normal_pos: 
    *a = DBL_MAX;
    *b = -DBL_MAX;
    for (m = 0; m < smo->M; m++) {
      mue = smo->s[state].mue[m];
      delta = 3 * sqrt(smo->s[state].u[m]);
      if (*a > mue - delta) *a = floor(mue - delta);
      if (*b < mue + delta) *b = ceil(mue + delta);
    }
    break;
  default: mes(MES_WIN, "Warning: density function not specified!\n"); 
  }
  if (smo->density == normal_pos && *a < 0.0) 
    *a = 0.0;
  return;
} /* smodel_get_interval_B */

/*============================================================================*/
double smodel_ifunc(smodel *smo, int state, double c, double x) {
  return( fabs( smodel_calc_B(smo, state, x) - c));
}

