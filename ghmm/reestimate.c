/*******************************************************************************
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/reestimate.c
  created      : TIME: 11:54:33     DATE: Mon 16. February 1998
  last-modified: TIME: 09:44:08     DATE: Wed 14. March 2001
*******************************************************************************/
/* $Id$ */

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
  double *pi_zaehler;
  double pi_nenner;
  double **a_zaehler;
  double *a_nenner;
  double **b_zaehler;
  double *b_nenner;
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
  if (!m_calloc(r->pi_zaehler, mo->N)) {mes_proc(); goto STOP;}
  /* Allozieren der a_zaehler -> individuelle Laenge */
  if (!m_calloc(r->a_zaehler, mo->N)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo->N; i++)
    if (!m_calloc(r->a_zaehler[i], mo->s[i].out_states)){mes_proc(); goto STOP;}
  if (!m_calloc(r->a_nenner, mo->N)) {mes_proc(); goto STOP;}
  r->b_zaehler = matrix_d_alloc(mo->N, mo->M);
  if (!(r->b_zaehler)) {mes_proc(); goto STOP;}
  if (!m_calloc(r->b_nenner, mo->N)) {mes_proc(); goto STOP;}
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
  m_free((*r)->pi_zaehler);  
  for (i = 0; i < N; i++)
    m_free((*r)->a_zaehler[i]);
  m_free((*r)->a_zaehler);
  m_free((*r)->a_nenner);
  matrix_d_free( &((*r)->b_zaehler), N );
  m_free((*r)->b_nenner);
  m_free(*r);
  return(0);
# undef CUR_PROC
} /* reestimate_free */

/*----------------------------------------------------------------------------*/
static int reestimate_init(local_store_t *r, const model *mo) {
# define CUR_PROC "reestimate_init"
  int i, j, m;
  r->pi_nenner = 0.0;
  for (i = 0; i < mo->N; i++) {
    r->pi_zaehler[i] = 0.0;
    r->a_nenner[i] = 0.0;
    r->b_nenner[i] = 0.0;
    for (j = 0; j < mo->s[i].out_states; j++)
      r->a_zaehler[i][j] = 0.0;
    for (m = 0; m < mo->M; m++)
      r->b_zaehler[i][m] = 0.0;
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
  int h, i, j, m, l, j_id, positiv;
  double faktor, p_i;
  mes_check_0(r->pi_nenner, goto STOP); 
  for (i = 0; i < mo->N; i++) {
    /* Pi */
    mo->s[i].pi =  r->pi_zaehler[i] / r->pi_nenner;
    /* A */
    /* Beachte: Nenner koennen 0 werden; Frage: Zustand NIE angenommen?*/
    p_i = 0.0;
    if (r->a_nenner[i] < EPS_PREC) {
      /* Test: eigentlich muesste auch Summe der in_a == 0 sein: */
      for (h = 0; h < mo->s[i].in_states; h++) {
	p_i += mo->s[i].in_a[h];
      }	
      if (p_i == 0) {
	char *str;
	if (mo->s[i].in_states == 0)
	  str = mprintf(NULL,0,
	    "Zustand %d kann nicht angenommen werden (keine in_states)\n",i);
	else
	  str = 
	    mprintf(NULL,0,"Zustand %d kann nicht angenommen werden (P=0)\n",i);
	mes_prot(str);
	m_free(str);
      }
      /* erstmal lassen (beser: aij gar nicht neu bestimmen?): */
      faktor = 0.0;
    }
    else 
      faktor = ( 1 / r->a_nenner[i] );      
    
    positiv = 0;
    
    for (j = 0; j < mo->s[i].out_states; j++) {
      /* TEST: Check, ob Zeahler < Nenner */
      if ((r->a_nenner[i] - r->a_zaehler[i][j]) <= -EPS_PREC)
	{ mes_prot("Zahler a groesser als Nenner!\n"); }
      mo->s[i].out_a[j] = r->a_zaehler[i][j] * faktor;

      if (r->a_zaehler[i][j] >= EPS_PREC) positiv = 1;

      /* Wichtig: in_a des entspr. Zustands auch aendern */
      l = 0;
      j_id = mo->s[i].out_id[j];
      while (l < mo->s[j_id].in_states)
	if  (mo->s[j_id].in_id[l] == i)
	  break;
        else 
	  l++;
      if ( l == mo->s[j_id].in_states) { 
	mes_prot("passendes in_a zu out_a nicht gefunden!\n"); 
	goto STOP;
      }
      mo->s[j_id].in_a[l] = mo->s[i].out_a[j]; 
    }

    if (!positiv) {
      char *str = 
	mprintf(NULL, 0, 
		"Alle Zaehler a[%d][j] == 0 (Nenner=%.4f, P(in)=%.4f)!\n", 
		i, r->a_nenner[i], p_i);
      mes_prot(str);
      m_free(str);
    }
    
    /* B */
    if (r->b_nenner[i] < EPS_PREC)
      faktor = 0.0;
    else 
      faktor = ( 1.0 / r->b_nenner[i] );      

    positiv = 0;

    for (m = 0; m < mo->M; m++) {
      /* TEST: Check, ob Zaehler < Nenner */
      if ((r->b_nenner[i] - r->b_zaehler[i][m]) <= -EPS_PREC) {
	char *str = 
	  mprintf(NULL, 0, "Zaehler b (%.4f) groesser als Nenner (%.4f)!\n", 
		  r->b_zaehler[i][m], r->b_nenner[i]); 
	mes_prot(str);
	m_free(str);
      }
      mo->s[i].b[m] = r->b_zaehler[i][m] * faktor;

      if (r->b_zaehler[i][m] >= EPS_PREC) positiv = 1;
      
    }

    if (!positiv) {
      char *str = 
	mprintf(NULL, 0, "Alle Zaehler b[%d][m] == 0 (Nenner=%.4f)!\n", 
		i, r->b_nenner[i]);
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
  int k, i, j, m, t, j_id, gueltig;
  double **alpha = NULL;
  double **beta = NULL;
  double *scale = NULL;
  int T_k;
  double gamma;
  double log_p_k;

  *log_p = 0.0;
  gueltig = 0;
  /* Schleife ueber alle Ausgabestrings */
  for (k = 0; k < seq_number; k++) {
    
    T_k = seq_length[k]; /* aktuelle Seq.laenge -> max. Laenge EINMAL alloz.? */
    
    /* Initialisieren der Matritzen und des Vektors abhaengig von T_k */
    if ( reestimate_alloc_matvek(&alpha, &beta, &scale, T_k, mo->N) == -1 ) 
      {mes_proc(); goto STOP;}
    
    /* Forward fuer aktuelles O(k) */
    if (foba_forward(mo, O[k], T_k, alpha, scale, &log_p_k) == -1) 
      {mes_proc(); goto STOP;}    

    if (log_p_k != +1) { /* d.h. O(k) kann von mo erzeugt werden */ 

      *log_p += log_p_k;
      /* printf("log P(O|lambda) = %8.5f  (vor reestimate)\n", log_p_k); */
      gueltig = 1;
      
      /* Backward fuer aktuelles O(k) */
      if (foba_backward(mo, O[k], T_k, beta, scale) == -1) 
	{ mes_proc(); goto STOP;}
      
      /* Schleife ueber alle Zustaende */
      for (i = 0; i < mo->N; i++) {
	/* Pi */
	r->pi_zaehler[i] += seq_w[k] * alpha[0][i] * beta[0][i];
	r->pi_nenner += seq_w[k] * alpha[0][i] * beta[0][i]; /* Summe ueber alle i */
	/* A */
	for (t = 0; t < T_k - 1; t++) {
	  r->a_nenner[i] += seq_w[k] * alpha[t][i] * beta[t][i]; /* unabh. von j ! */
	  for (j = 0; j < mo->s[i].out_states; j++) {
	    j_id = mo->s[i].out_id[j]; /* aufpassen! */
	    r->a_zaehler[i][j] += 
	      ( seq_w[k] * alpha[t][i] 
		* mo->s[i].out_a[j] 
		* mo->s[j_id].b[ O[k][t+1] ]
		* beta[t+1][j_id] 
		* (1.0 / scale[t+1]) );  /* c[t] = 1/scale[t] */
	  }
	}
	/* B */
	for (t = 0; t < T_k; t++) {
	  gamma = (seq_w[k] * alpha[t][i] * beta[t][i] );
	  r->b_nenner[i] += gamma;                   /* unabh. von m ! */
	  for (m = 0; m < mo->M; m++) {
	    if ( O[k][t] == m )
	      r->b_zaehler[i][m] += gamma;
	  }
	}
      } /* for (i = 0 i < mo->N i++) { */

    } /* if (log_p_k != +1) */
    else {
      printf("O(%2d) can't be build from model mo!\n", k);
    } 

    reestimate_free_matvek(alpha, beta, scale, T_k);
    
  } /* for (k = 0; k < seq_number; k++) */
  
  if (gueltig) {
    /* neues Lambda setzen: Modell mo DIREKT aendern !!! */
    if ( reestimate_setlambda(r, mo) == -1 ) { mes_proc(); goto STOP; } 

    /* TESTPHASE: Plausibilitaets-Check des neuen Modells
    if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
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
  int n, k, gueltig;
  double log_p, log_p_old, log_p_k, diff;
  local_store_t *r = NULL;
 
  /* einmal local_store_t *r bereitstellen fuer alle Iterationen */
  r = reestimate_alloc(mo);
  if (!r) {mes_proc(); goto STOP; };
  
  log_p_old = -DBL_MAX; 
  n = 1;
  while (n <= MAX_ITER_BW) {  /* maximale Iterationsanzahl */
    
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

    /* Konvergenzfehler */
    diff = log_p - log_p_old;
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
      
    /* Abbruchbedingung */
    if ( diff < fabs((double)EPS_ITER_BW * log_p) ) {
      printf("Convergence after %d steps\n", n); 
      break;
    }
    else {
      log_p_old = log_p;
      /* alle Felder in r wieder auf 0.0 setzen fuer naechste Iteration */
      reestimate_init(r, mo);
      n++;
    }

  } /* while (n <= MAX_ITER) */ 
  
  /* log_p des endgueltigen Modells berechnen */
  log_p = 0.0;
  gueltig = 0;
  for (k = 0; k < sq->seq_number; k++) {
    if (foba_logp(mo, sq->seq[k], sq->seq_len[k], &log_p_k) == -1) 
      { mes_proc(); goto STOP; }
    if (log_p_k != +1) {
      log_p += log_p_k;
      gueltig = 1;
    }
  }
  if (!gueltig)
    log_p = +1;
  /* printf("log P(O|lambda) = %8.5f  (n = %d)\n", log_p, n); */
  printf("%8.5f (-log_p optimized model)\n", -log_p);
  
  /* Plausibilitaetscheck der neuen Werte
  if (model_check(mo) == -1) { mes_proc(); goto STOP; } */
  return(0);

STOP:
  reestimate_free(&r, mo->N);
  return(-1);
# undef CUR_PROC
} /* reestimate_baum_welch */
