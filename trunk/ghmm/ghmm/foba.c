/*******************************************************************************
  author       : Bernd Wichern
  filename     : /homes/hmm/bknab/src/foba.c
  created      : TIME: 00:00:00     DATE: Tue 00. xxx 0000
  last-modified: TIME: 10:47:02     DATE: Tue 16. November 1999
*******************************************************************************/
#include "foba.h"
#include "math.h"
#include "const.h"
#include "matrix.h"
#include "mes.h"

static int foba_initforward(model *mo, double *alpha_1, int symb, 
			    double *scale) {
# define CUR_PROC "foba_initforward"
  int i;
  double c_0;
  scale[0] = 0.0;
  for (i = 0; i < mo->N; i++) {
    alpha_1[i] = mo->s[i].pi * mo->s[i].b[symb];
    scale[0] += alpha_1[i];
  }
  if (scale[0] >= EPS_PREC) {
    c_0 = 1/scale[0];
    for (i = 0; i < mo->N; i++) 
      alpha_1[i] *= c_0;
  }
  return(0); /* BEACHTE: scale[0] kann jetzt 0 sein! */
# undef CUR_PROC
} /* foba_initforward */

/*----------------------------------------------------------------------------*/

static double foba_stepforward(state *s, double *alpha_t, double b_symb) {
  int i, id;
  double value = 0.0;

  for (i = 0; i < s->in_states; i++) {
    id = s->in_id[i];
    value += s->in_a[i] * alpha_t[id];
  }
  value *= b_symb; /* b_symb vor die Summe ziehen */
  return(value);

} /* foba_stepforward */

/*============================================================================*/

int foba_forward(model *mo, const int *O, int len, double **alpha, 
		 double *scale, double *log_p) {
# define CUR_PROC "foba_forward"
  int i, t;
  double c_t;
  foba_initforward(mo, alpha[0], O[0], scale);
  if (scale[0] < EPS_PREC) {
    /* d.h. 1.Zeichen kann vom HMM nicht erzeugt werden */
    *log_p = +1;
  }
  else {
    *log_p = - log(1/scale[0]);
    for (t = 1; t < len; t++) {
      scale[t] = 0.0;
      for (i = 0; i < mo->N; i++) {
	alpha[t][i] = foba_stepforward(&mo->s[i], alpha[t-1], mo->s[i].b[O[t]]);
	scale[t] +=  alpha[t][i];
      }
      if (scale[t] < EPS_PREC) {
	/* d.h. O-string kann vom HMM nicht erzeugt werden */
	*log_p = +1.0;
	break;
      }
      c_t = 1/scale[t];
      for (i = 0; i < mo->N; i++) 
	alpha[t][i] *= c_t;
      /* aufsummieren der log(c[t]) zur Berechnung log( P(O|lambda) ) */
      *log_p -= log(c_t);
    }
  }
  return 0;
# undef CUR_PROC
} /* foba_forward */

/*============================================================================*/

int foba_backward(model *mo, const int *O, int len, double **beta, const double *scale) {
# define CUR_PROC "foba_backward"
  double *beta_tmp, sum;
  int i, j, j_id, t;
  int res = -1;
  if (!m_calloc(beta_tmp, mo->N)) {mes_proc(); goto STOP;}
  for (t = 0; t < len; t++)
    mes_check_0(scale[t], goto STOP);
  /* Initialisierung */
  for (i = 0; i < mo->N; i++) {
    beta[len-1][i] = 1;
    beta_tmp[i] = 1/scale[len-1];
  }

  /* Backward Step for t = T-2, ..., 0 */
  /* beta_tmp: Vektor zum Zwischenspeichern der skalierten beta in einem 
     Zeitschritt */
  for (t = len-2; t >= 0; t--) {
    for (i = 0; i < mo->N; i++) {
      sum = 0.0;
      for (j = 0; j < mo->s[i].out_states; j++) {
	j_id = mo->s[i].out_id[j];
	sum += mo->s[i].out_a[j] * mo->s[j_id].b[O[t+1]] * beta_tmp[j_id];
      }
      beta[t][i] = sum;
    }
    for (i = 0; i < mo->N; i++) 
      beta_tmp[i] = beta[t][i]/scale[t];
  }
  res = 0;
STOP:
  m_free(beta_tmp);
  return(res);
# undef CUR_PROC
} /* foba_backward */


/*============================================================================*/
int foba_logp(model *mo, const int *O, int len, double *log_p) {
# define CUR_PROC "foba_logp"
  int res = -1;
  double **alpha, *scale = NULL;
  alpha = matrix_d_alloc(len, mo->N);
  if (!alpha) {mes_proc(); goto STOP;}
  if (!m_calloc(scale, len)) {mes_proc(); goto STOP;}
  /* forward durchlaufen lassen */
  if (foba_forward(mo, O, len, alpha, scale, log_p) == -1 )
    { mes_proc(); goto STOP; }
  res = 0;
STOP:
  matrix_d_free(&alpha, len);
  m_free(scale);
  return(res);
# undef CUR_PROC
} /* foba_logp */
