/*******************************************************************************
  author       : Bernd Wichern
  filename     : /homes/hmm/wichern/hmm/src/smap_classify.c
  created      : TIME: 15:57:29     DATE: Mon 03. January 2000
  last-modified: TIME: 16:40:45     DATE: Wed 03. January 2001
*******************************************************************************/

#include <math.h>
#include <float.h>
#include "mes.h"
#include "sfoba.h"
#include "const.h"
#include "matrix.h"
#include "smap_classify.h"

typedef struct local_store_t {
  double ***alpha;
  double **scale;
  double *prob;
  double **p;
  double **sum_alpha;
  double *alpha_1; 
  double *prior;
  int    *error;
} local_store_t;

/*----------------------------------------------------------------------------*/
static local_store_t *smap_classify_alloc(int mo_number, int states, int T) {
#define CUR_PROC "smap_classify_alloc"
  int i;
  local_store_t *map = NULL;
  if (!m_calloc(map, 1)) {mes_proc(); goto STOP;}
  if (!m_calloc(map->alpha, mo_number)) {mes_proc(); goto STOP;}
  for (i = 0; i < mo_number; i++) {
    map->alpha[i] = matrix_d_alloc(T, states);
    if (!map->alpha[i]) {mes_proc(); goto STOP;}
  }
  map->scale = matrix_d_alloc(mo_number, T);
  if (!map->scale) {mes_proc(); goto STOP;}
  map->p = matrix_d_alloc(mo_number, T + 1);
  if (!map->p) {mes_proc(); goto STOP;}
  map->sum_alpha = matrix_d_alloc(mo_number, T);
  if (!map->sum_alpha) {mes_proc(); goto STOP;}
  if (!m_calloc(map->prob, mo_number)) {mes_proc(); goto STOP;}
  if (!m_calloc(map->alpha_1, mo_number)) {mes_proc(); goto STOP;}
  if (!m_calloc(map->error, mo_number)) {mes_proc(); goto STOP;}
  if (!m_calloc(map->prior, mo_number)) {mes_proc(); goto STOP;}
  return map;
 STOP:
  return NULL;
#undef CUR_PROC
}
/*----------------------------------------------------------------------------*/
static int smap_classify_free(local_store_t **map, int mo_number, int T) {
# define CUR_PROC "smap_classify_free"
  int i;
  mes_check_ptr(map, return(-1));
  if( !*map ) return(0);
  m_free((*map)->prob);
  m_free((*map)->alpha_1);
  m_free((*map)->error);
  m_free((*map)->prior);
  matrix_d_free(&((*map)->scale), mo_number);
  matrix_d_free(&((*map)->p), mo_number);
  matrix_d_free(&((*map)->sum_alpha), mo_number);
  for (i = 0; i < mo_number; i++)
    matrix_d_free(&((*map)->alpha[i]), T);
  m_free(*map);

  return 0;
#undef CUR_PROC
}

/*============================================================================*/

/* Maximum A Posteriori Classification Algorithm (MAPCA): 
   given a field of models and one sequence and suppose the sequence
   has been produced by one of these models. This algorithm calculates
   for each model the probability, that the seq. comes from the model.
   This bayesian approach uses a prior for the models. If none is specified
   equal prob. is assumed.
   The maps are copied into "result", which has to be of dimension "smo_number"
   Ref.: A. Kehagias: Bayesian Classification of HMM, Math. Comp. Modelling
   (1995)
*/
  
int smap_classify(smodel **smo, double *result, int smo_number,
		      double *O, int T) {
#define CUR_PROC "smap_classify"
  int t, n, m, max_model = -1, build = 0; /* time, states, models */
  double log_p, denom_26, sum = 0.0, max_result = 0;
  local_store_t *map = NULL;


  if(smo == NULL || smo_number <= 0 || O == NULL || T < 0)
    {mes_proc(); goto STOP;}

  map = smap_classify_alloc(smo_number, smo[0]->N, T);

  /* no prior defined --> same prior for all models */
  for (m = 0; m < smo_number; m++)
    if (smo[m]->prior == -1)
      map->prior[m] = 1 / (double) smo_number;
    else
      map->prior[m] = smo[m]->prior;
  
  for (m = 0; m < smo_number; m++)
    sum += map->prior[m];
  if (fabs(1 - sum) > 0.0001) {
    mes_prot("Sum of model priors != 1 or mixing of priors and non-priors \n");
    goto STOP;
  }

  
  /* for all models */
  for (m = 0; m < smo_number; m++) {
    map->p[m][0] = map->prior[m];                    /* (22) */
    map->alpha_1[m] = 1 / (double) smo_number;       /* (23) */
    /* forward alg.                                     (24) */
    if (sfoba_forward(smo[m], O, T, NULL, map->alpha[m], map->scale[m],
		      &log_p) == -1) {
      /* Sequenz kann nicht von dem Modell nicht generiert werden */
      map->error[m] = 1;
      continue;
    }
    build = 1;
    for (t = 0; t < T; t++) {
      for (n = 0; n < smo[0]->N; n++)
	map->sum_alpha[m][t] += map->alpha[m][t][n]; /* Sums (25) */
      if (map->sum_alpha[m][t] == 0)
	{mes_prot("sum_alpha = 0\n"); goto STOP;} /* andere Loesung ??? */
    }
  }
  if (!build) {
    mes_prot("Prob. = 0 for all models\n"); goto STOP;
  }
  
  for (t = 0; t < T; t++) {
    denom_26 = 0.0;
    for (m = 0; m < smo_number; m++) {
      if (!map->error[m]) {
	if (t != 0)
	  map->prob[m] = map->sum_alpha[m][t] * map->scale[m][t] / /* (25) */
	    map->sum_alpha[m][t - 1];
	else
	  map->prob[m] = map->sum_alpha[m][t] * map->scale[m][t];  /* (25) */
	denom_26 += map->prob[m] * map->p[m][t];
      }
      else
	map->prob[m] = 0.0; /* Seq. cannot be generated by this model */
    }
    for (m = 0; m < smo_number; m++) {
      if (denom_26 != 0)
	map->p[m][t + 1] = map->prob[m] * map->p[m][t] / denom_26; /* (26) */
      else 
	{mes_prot("denom_26 == 0!\n"); goto STOP;}
    }
  }
  /* copy into result and find best model */
  max_result = 0.0;
  for (m = 0; m < smo_number; m++) {
    result[m] = map->p[m][T];
    if (result[m] > max_result) {
      max_result = result[m];
      max_model = m;
    }
  }
			     
 STOP:
  smap_classify_free(&map, smo_number, T);
  return max_model;

#undef CUR_PROC
}

/* Alternative zu MAPCA (smap_classify) direkte Berechnung von p[m] ueber 
   Bayes Formel  anstatt Rekursion ueber t. 
   p(m | O) = p(O | m) * p(m) / (sum_i p(O | i) * p(i))
   Gleiches Resultat?
   
*/

int smap_bayes(smodel **smo, double *result, int smo_number, double *O, int T) {
#define CUR_PROC "smap_bayes"
  double *prior, *log_p;
  double sum = 0.0, p_von_O = 0.0, max_result = 0.0;
  int *error;
  int found_model = 0, err = 0;
  int m, max_model = -1;

  /* nothing to do here */
  if (smo_number == 1) {
    result[0] = 1.0;
    return 0;
  }

  for (m = 0; m < smo_number; m++)
    result[m] = 0;

  if (!m_calloc(prior, smo_number)) {mes_proc(); goto STOP;}
  if (!m_calloc(log_p, smo_number)) {mes_proc(); goto STOP;}
  if (!m_calloc(error, smo_number)) {mes_proc(); goto STOP;}

 if(smo == NULL || smo_number <= 0 || O == NULL || T < 0)
    {mes_proc(); goto STOP;}

 /* no prior defined --> same prior for all models */
 for (m = 0; m < smo_number; m++)
   if (smo[m]->prior == -1)
     prior[m] = 1 / (double) smo_number;
   else
     prior[m] = smo[m]->prior;
 
 for (m = 0; m < smo_number; m++)
   sum += prior[m];
 if (fabs(1 - sum) > 0.0001) {
   mes_prot("Sum of model priors != 1 or mixing of priors and non-priors \n");
   for (m = 0; m < smo_number; m++)
     printf("%.6f  ", prior[m]);
   printf("\n");
   goto STOP;
 }
 
 /* log_p berechnen, jedes Modell; verknuepft mit prior ergibt
  die Likelihood von O, gegeben alle Modelle */
 for (m = 0; m < smo_number; m++)
   if (sfoba_logp(smo[m], O, T, &log_p[m]) == -1) {
     error[m] = 1;
   }
   else {
     error[m] = 0;
     p_von_O += exp(log_p[m]) * prior[m];
     found_model = 1;
   }
 
 if (p_von_O == 0) {
   mes_prot("P(O) = 0!\n"); 
   err = 1;
 }
 if (!found_model) {
   mes_prot("-1 from sfoba_logp for all models\n"); 
   err = 1;
 }
 if (err)  goto STOP;
 
 /* Likelihood der Modelle, gegeben O */
 for (m = 0; m < smo_number; m++) {
   if (!error[m]) {
     result[m] = exp(log_p[m]) * prior[m] / p_von_O;
     if (result[m] > max_result) {
       max_model = m;
       max_result = result[m];
     }
   }   
 }
 /* TEST */
 if (max_model == -1) {
   printf("smap_bayes: max_model == -1\n");
   for (m = 0; m < smo_number; m++) 
     if (!error[m])
       printf("%f %.4f %.4f;  ", log_p[m], prior[m], p_von_O);
   printf("\n");
 }
  STOP:

  m_free(prior);
  m_free(log_p);
  m_free(error);
  return max_model;
#undef CUR_PROC
}
