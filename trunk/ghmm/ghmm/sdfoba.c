/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sdfoba.c
*       Authors:  Utz J. Pape
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


#include <math.h>

#include "ghmm.h"
#include "sdmodel.h"
#include "matrix.h"
#include "mes.h"
#include "ghmm_internals.h"


static int sdfoba_initforward (sdmodel * mo, double *alpha_1, int symb,
                               double *scale)
{
# define CUR_PROC "sdfoba_initforward"
  int i, j, id, in_id;
  int class = 0;
  double c_0;
  scale[0] = 0.0;
  /*iterate over non-silent states*/
  for (i = 0; i < mo->N; i++) {
    if (!(mo->silent[i])) {
      if (symb != mo->M) {
        alpha_1[i] = mo->s[i].pi * mo->s[i].b[symb];
        /*      printf("\nalpha1[%i]=%f\n",i,alpha_1[i]);*/
      }

      else {
        alpha_1[i] = mo->s[i].pi;
      }
      scale[0] += alpha_1[i];

    }
  }
  /*iterate over silent states*/
  for (i = 0; i < mo->topo_order_length; i++) {
    id = mo->topo_order[i];
    alpha_1[id] = mo->s[id].pi;
    /*      printf("\nsilent_start alpha1[%i]=%f\n",id,alpha_1[id]);*/
    for (j = 0; j < mo->s[id].in_states; j++) {
      in_id = mo->s[id].in_id[j];
      alpha_1[id] += mo->s[id].in_a[class][j] * alpha_1[in_id];
    }
    /*      printf("\n\tsilent_run alpha1[%i]=%f\n",id,alpha_1[id]);*/
    scale[0] += alpha_1[id];
  }
  /*  printf("\n%f\n",scale[0]);*/
  if (scale[0] >= EPS_PREC) {
    c_0 = 1 / scale[0];
    for (i = 0; i < mo->N; i++)
      alpha_1[i] *= c_0;
  }

  return (0);                   /* attention scale[0] might be 0 */
# undef CUR_PROC
}                               /* sdfoba_initforward */

/*----------------------------------------------------------------------------*/

static double sdfoba_stepforward (sdstate * s, double *alpha_t,
                                  const double b_symb, int class)
{
  int i, id;
  double value = 0.0;

  for (i = 0; i < s->in_states; i++) {
    id = s->in_id[i];
    /*    printf("state:\t%s, instate:\t%d, prob:\t%f\n", s->label, id, s->in_a[class][i]);*/
    value += s->in_a[class][i] * alpha_t[id];
  }
  value *= b_symb;
  return (value);

}                               /* sdfoba_stepforward */

/*============================================================================*/

int sdfoba_forward (sdmodel * mo, const int *O, int len, double **alpha,
                    double *scale, double *log_p)
{
# define CUR_PROC "sdfoba_forward"
  int i, t, id;
  double c_t, dblems;
  int class = 0;

  /*if (mo->model_type & kSilentStates)
     sdmodel_topo_ordering(mo);
   */
  sdfoba_initforward (mo, alpha[0], O[0], scale);
  if (scale[0] < EPS_PREC) {
    /* means: first symbol can't be generated by hmm */
    printf ("\nnach init gestoppt\n");
    *log_p = +1;
  }
  else {
    *log_p = -log (1 / scale[0]);
    for (t = 1; t < len; t++) {
      scale[t] = 0.0;
      /*      printf("\nStep t=%i mit len=%i, O[i]=%i\n",t,len,O[t]);*/
      if (mo->cos > 1)
        class = mo->get_class (&(mo->N), t-1);
      /*iterate over non-silent states*/
      /*printf("\nnach Class\n");*/
      for (i = 0; i < mo->N; i++) {
        if (!(mo->model_type & kSilentStates) || !(mo->silent[i])) {
          if (O[t] != mo->M) {
            dblems = mo->s[i].b[O[t]];
          }
          else {
            dblems = 1.0;
          }
          alpha[t][i] =
            sdfoba_stepforward (&mo->s[i], alpha[t - 1], dblems, class);

          /*           printf("alpha[%i][%i] = %f\n", t, i, alpha[t][i]);*/
          scale[t] += alpha[t][i];
          /*printf("\nalpha[%d][%d] = %e, scale[%d] = %e", t,i, alpha[t][i], t, scale[t]);*/

          /*printf("scale[%i] = %f\n", t, scale[t]);*/
        }
      }
      /*printf("\nvor silent states\n");*/
      /*iterate over silent state*/
      if (mo->model_type & kSilentStates) {
        for (i = 0; i < mo->topo_order_length; i++) {
          /*printf("\nget id\n");*/
          id = mo->topo_order[i];
          /*printf("\nin stepforward\n");*/
          alpha[t][id] = sdfoba_stepforward (&mo->s[id], alpha[t], 1, class);
          /*   printf("alpha[%i][%i] = %f\n", t, id, alpha[t][id]);*/
          /*printf("\nnach stepforward\n");*/
          scale[t] += alpha[t][id];
          /*printf("\nalpha[%d][%d] = %e, scale[%d] = %e", t,id, alpha[t][id], t, scale[t]);*/
          /*printf("silent state: %d\n", id);*/
        }
      }
      /*printf("\nnach silent states\n");*/
      if (scale[t] < EPS_PREC) {
        /*char *str = */
        printf ("numerically questionable: ");
        /*mes_prot(str);*/
        /*printf("\neps = %e\n", EPS_PREC);*/
        /*printf("\nscale kleiner als eps HUHU, Zeit t = %d, scale = %e\n", t, scale[t]);*/
        /*printf("\nlabel: %s, char: %d, ems: %f\n", mo->s[92].label,O[t], mo->s[4].b[O[t]]);*/
        /*vector_d_print(stdout, alpha[t],mo->N, "\t", " ", "\n");*/
        /* O-string  can't be generated by hmm */
        /*      *log_p = +1.0;*/
        /*break;*/
      }
      c_t = 1 / scale[t];
      for (i = 0; i < mo->N; i++)
        alpha[t][i] *= c_t;
      /* sum log(c[t]) to get  log( P(O|lambda) ) */
      *log_p -= log (c_t);
    }
  }

  return 0;
# undef CUR_PROC
}                               /* sdfoba_forward */

/*============================================================================*/
static int sdfobau_initforward (sdmodel * mo, double *alpha_1, int symb,
                                double *scale)
{
# define CUR_PROC "sdfoba_initforward"
  int i, j, id, in_id;
  int class = 0;
  double c_0;
  scale[0] = 0.0;
  /*iterate over non-silent states*/
  for (i = 0; i < mo->N; i++) {
    if (!(mo->silent[i])) {
      alpha_1[i] = mo->s[i].pi * mo->s[i].b[symb];
      /*printf("\nalpha1[%i]=%f\n",i,alpha_1[i]);*/
      scale[0] += alpha_1[i];
    }
  }
  /*iterate over silent states*/
  for (i = 0; i < mo->topo_order_length; i++) {
    id = mo->topo_order[i];
    alpha_1[id] = mo->s[id].pi;
    /*printf("\nsilent_start alpha1[%i]=%f\n",id,alpha_1[id]);*/
    for (j = 0; j < mo->s[id].in_states; j++) {
      in_id = mo->s[id].in_id[j];
      alpha_1[id] += mo->s[id].in_a[class][j] * alpha_1[in_id];
      /*printf("\n\tsilent_run alpha1[%i]=%f\n",id,alpha_1[id]);*/
    }
    /*scale[0] += alpha_1[id];*/
  }
  /*printf("\n%f\n",scale[0]);*/
  if (scale[0] >= EPS_PREC) {
    c_0 = 1 / scale[0];
    for (i = 0; i < mo->N; i++)
      alpha_1[i] *= c_0;
  }
  return (0);                   /* attention scale[0] might be 0 */
# undef CUR_PROC
}                               /* sdfoba_initforward */

/*----------------------------------------------------------------------------*/

int sdfobau_forward (sdmodel * mo, const int *O, int len, double **alpha,
                     double *scale, double *log_p)
{
# define CUR_PROC "sdfoba_forward"
  int i, t, id;
  double c_t;
  int class = 0;

  if (mo->model_type & kSilentStates)
    sdmodel_topo_ordering (mo);

  sdfobau_initforward (mo, alpha[0], O[0], scale);
  if (scale[0] < EPS_PREC) {
    /* means: first symbol can't be generated by hmm */
    *log_p = +1;
  }
  else {
    *log_p = -log (1 / scale[0]);
    for (t = 1; t < len; t++) {
      scale[t] = 0.0;
      if (mo->cos > 1)
        class = mo->get_class (&(mo->N), t - 1);
      /*iterate over non-silent states*/
      for (i = 0; i < mo->N; i++) {
        if (!(mo->model_type & kSilentStates) || !(mo->silent[i])) {
          alpha[t][i] =
            sdfoba_stepforward (&mo->s[i], alpha[t - 1], mo->s[i].b[O[t]],
                                class);
          scale[t] += alpha[t][i];
        }
      }
      /*iterate over silent state       */
      if (mo->model_type & kSilentStates) {
        for (i = 0; i < mo->topo_order_length; i++) {
          id = mo->topo_order[i];
          alpha[t][id] = sdfoba_stepforward (&mo->s[id], alpha[t], 1, class);
          /*scale[t] += alpha[t][id];*/
        }
      }
      if (scale[t] < EPS_PREC) {
        /* O-string  can't be generated by hmm */
        *log_p = +1.0;
        break;
      }
      c_t = 1 / scale[t];
      for (i = 0; i < mo->N; i++)
        alpha[t][i] *= c_t;
      /* sum log(c[t]) to get  log( P(O|lambda) ) */
      *log_p -= log (c_t);
    }
  }

  return 0;
# undef CUR_PROC
}                               /* sdfoba_forward */

/*============================================================================*/

int sdfoba_descale (double **alpha, double *scale, int t, int n,
                    double **newalpha)
{
# define CUR_PROC "sdfoba_descale"
  int i, j, k;
  /*printf("\nAngekommen, t=%i, n=%i\n",t,n);*/
  for (i = 0; i < t; i++) {
    /*printf("i=%i\n",i);*/
    for (j = 0; j < n; j++) {
      /*printf("\tj=%i\n",j);*/
      newalpha[i][j] = alpha[i][j];
      /*newalpha[i][j] *= scale[j];*/
      /*printf(".");*/
      for (k = 0; k <= i; k++) {
        /*printf(",");*/
        newalpha[i][j] *= scale[k];
      }
    }
  }
  /*printf("\ndescale geschafft\n");*/
  return 0;
# undef CUR_PROC
}                               /* sdfoba_descale */

/*============================================================================*/

int sdfoba_backward (sdmodel * mo, const int *O, int len, double **beta,
                     const double *scale)
{
# define CUR_PROC "sdfoba_backward"
  double *beta_tmp, sum;
  int i, j, j_id, t;
  int res = -1;
  ARRAY_CALLOC (beta_tmp, mo->N);
  for (t = 0; t < len; t++)
    mes_check_0 (scale[t], goto STOP);
  /* initialize */
  for (i = 0; i < mo->N; i++) {
    beta[len - 1][i] = 1;
    beta_tmp[i] = 1 / scale[len - 1];
  }

  /* Backward Step for t = T-2, ..., 0 */
  /* beta_tmp: Vector for storage of scaled beta in one time step */
  for (t = len - 2; t >= 0; t--) {
    for (i = 0; i < mo->N; i++) {
      sum = 0.0;
      for (j = 0; j < mo->s[i].out_states; j++) {
        j_id = mo->s[i].out_id[j];
        /*sum += mo->s[i].out_a[j] * mo->s[j_id].b[O[t+1]] * beta_tmp[j_id];*/
      }
      beta[t][i] = sum;
    }
    for (i = 0; i < mo->N; i++)
      beta_tmp[i] = beta[t][i] / scale[t];
  }
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  m_free (beta_tmp);
  return (res);
# undef CUR_PROC
}                               /* sdfoba_backward */


/*============================================================================*/
int sdfoba_logp (sdmodel * mo, const int *O, int len, double *log_p)
{
#define CUR_PROC "sdfoba_logp"
  int res = -1;
  double **alpha, *scale = NULL;
  alpha = matrix_d_alloc (len, mo->N);
  if (!alpha) {
    mes_proc ();
    goto STOP;
  }
  ARRAY_CALLOC (scale, len);
  /* run foba_forward */
  if (sdfoba_forward (mo, O, len, alpha, scale, log_p) == -1) {
    mes_proc ();
    goto STOP;
  }
  res = 0;
  matrix_d_free (&alpha, len);
  m_free (scale);
  return (res);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  matrix_d_free (&alpha, len);
  m_free (scale);
  return (res);
# undef CUR_PROC
}                               /* sdfoba_logp */
