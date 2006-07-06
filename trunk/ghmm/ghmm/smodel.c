/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/smodel.c
*       Authors:  Bernhard Knab, Benjamin Georgi
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

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "../config.h"
#endif

#include <float.h>
#include <math.h>
#include <assert.h>

#include "ghmm.h"
#include "mprintf.h"
#include "mes.h"
#include "smodel.h"
#include "sfoba.h"
#include "sreestimate.h"
#include "matrix.h"
#include "vector.h"
#include "sequence.h"
#include "rng.h"
#include "randvar.h"
#include "string.h"
#include "ghmm_internals.h"
#include "xmlreader.h"
#include "xmlwriter.h"

#include "obsolete.h"

/*----------------------------------------------------------------------------*/
int ghmm_c_state_alloc (ghmm_cstate * s,
                        int M, int in_states, int out_states, int cos)
{
# define CUR_PROC "ghmm_c_state_alloc"
  int res = -1;
  int i;
  ARRAY_CALLOC (s->c, M);
  ARRAY_CALLOC (s->mue, M);
  ARRAY_CALLOC (s->u, M);
  ARRAY_CALLOC (s->a,M);

  ARRAY_CALLOC (s->mixture_fix, M);
  ARRAY_CALLOC (s->density,M);
  s->xPosition = 0.0;
  s->yPosition = 0.0;
   

  /* mixture component fixing deactivated by default */
  for (i = 0; i < M; i++) {
    s->mixture_fix[i] = 0;
  }

  if (out_states > 0) {
    ARRAY_CALLOC (s->out_id, out_states);
    s->out_a = ighmm_cmatrix_alloc (cos, out_states);
    if (!s->out_a) {
      mes_proc ();
      goto STOP;
    }
  }
  if (in_states > 0) {
    ARRAY_CALLOC (s->in_id, in_states);
    s->in_a = ighmm_cmatrix_alloc (cos, in_states);
    if (!s->in_a) {
      mes_proc ();
      goto STOP;
    }
  }
  res = 0;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* ghmm_c_state_alloc */

/*----------------------------------------------------------------------------*/
ghmm_cmodel * ghmm_cmodel_calloc(int N, int modeltype) {
#define CUR_PROC "ghmm_cmodel_calloc"
  ghmm_cmodel * mo;

  assert(modeltype & GHMM_kContinuousHMM);
  ARRAY_CALLOC(mo, 1);
  mo->N = N;
  mo->M = 0;
  ARRAY_CALLOC(mo->s, N);
  return mo;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  ghmm_c_free(&mo);
  return NULL;
#undef CUR_PROC
}

int ghmm_c_class_change_alloc (ghmm_cmodel * smo)
{
#define CUR_PROC "ghmm_c_class_change_alloc"
  ghmm_c_class_change_context *c = NULL;
  ARRAY_CALLOC (c, 1);

  c->k = -1;
  c->get_class = NULL;

  smo->class_change = c;
  return (0);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (-1);
# undef CUR_PROC
}

#ifdef GHMM_OBSOLETE
/*----------------------------------------------------------------------------*/
static int smodel_copy_vectors (ghmm_cmodel * smo, int index, double *pi, int *fix,
                                double ***a_matrix, double **c_matrix,
                                double **mue_matrix, double **u_matrix)
{
#define CUR_PROC "smodel_alloc_vectors"
/* XXX - this function is only used in the smo file reading functions. It should be removed 
since it does not take care of het .mixtures densities  */
  int i, c, exist, count_out = 0, count_in = 0;
  smo->s[index].pi = pi[index];
  smo->s[index].fix = fix[index];
  for (i = 0; i < smo->s[index].M; i++) {
    smo->s[index].c[i] = c_matrix[index][i];
    smo->s[index].mue[i] = mue_matrix[index][i];
    smo->s[index].u[i] = u_matrix[index][i];
  }

  for (i = 0; i < smo->N; i++) {
    exist = 0;
    for (c = 0; c < smo->cos; c++) {
      if (a_matrix[c][index][i]) {
        exist = 1;
        break;
      }
    }
    /* transition to successor possible at least in one transition class */
    if (exist) {
      if (count_out >= smo->s[index].out_states) {
        mes_proc ();
        return (-1);
      }
      smo->s[index].out_id[count_out] = i;
      for (c = 0; c < smo->cos; c++)
        smo->s[index].out_a[c][count_out] = a_matrix[c][index][i];
      count_out++;
    }
    exist = 0;
    for (c = 0; c < smo->cos; c++) {
      if (a_matrix[c][i][index]) {
        exist = 1;
        break;
      }
    }
    /* transition to predecessor possible at least in one transition class */
    if (exist) {
      if (count_in >= smo->s[index].in_states) {
        mes_proc ();
        return (-1);
      }
      smo->s[index].in_id[count_in] = i;
      for (c = 0; c < smo->cos; c++)
        smo->s[index].in_a[c][count_in] = a_matrix[c][i][index];
      count_in++;
    }
  }
  return (0);
#undef CUR_PROC
}                               /* smodel_alloc_vectors */


/*============================================================================*/
ghmm_cmodel **ghmm_c_read (const char *filename, int *smo_number)
{
#define CUR_PROC "ghmm_c_read"
/* XXX - old function. No support to heterogeneous densities  */
  int j;
  long new_models = 0;
  scanner_t *s = NULL;
  ghmm_cmodel **smo = NULL;
  *smo_number = 0;
  s = ighmm_scanner_alloc (filename);
  if (!s) {
    mes_proc ();
    goto STOP;
  }
  while (!s->err && !s->eof) {
    ighmm_scanner_get_name (s);
    ighmm_scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    if (!strcmp (s->id, "SHMM") || !strcmp (s->id, "shmm")) {
      (*smo_number)++;
      /* more mem */
      ARRAY_REALLOC (smo, *smo_number);
      smo[*smo_number - 1] = ghmm_c_read_block (s, (int *) &new_models);
      if (!smo[*smo_number - 1]) {
        mes_proc ();
        goto STOP;
      }
      /* copy ghmm_cmodel */
      if (new_models > 1) {
        /* "-1" due to  (*smo_number)++ from above  */
        ARRAY_REALLOC (smo, *smo_number - 1 + new_models);
        for (j = 1; j < new_models; j++) {
          smo[*smo_number] = ghmm_c_copy (smo[*smo_number - 1]);
          if (!smo[*smo_number]) {
            mes_proc ();
            goto STOP;
          }
          (*smo_number)++;
        }
      }
    }
    else {
      ighmm_scanner_error (s, "unknown identifier");
      goto STOP;
    }
    ighmm_scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(!s->err && !s->eof) */
  /*if (ghmm_c_check_compatibility(smo, *smo_number) == -1) {
     mes_proc(); goto STOP;
     } */

  ighmm_scanner_free (&s);             
  return smo;
  
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
    ighmm_scanner_free (&s);
    return NULL;
#undef CUR_PROC
}    /* ghmm_c_read */


/*============================================================================*/
ghmm_cmodel *ghmm_c_read_block (scanner_t * s, int *multip)
{
#define CUR_PROC "ghmm_c_read_block"
/* XXX - old function. No support to heterogeneous densities  */
  int i, j, osc, m_read, n_read, pi_read, a_read, c_read, mue_read, cos_read,
    u_read, len, density_read, out, in, prior_read, fix_read, M;
  ghmm_cmodel *smo = NULL;
  double *pi_vektor = NULL, **a_matrix = NULL, ***a_3dmatrix = NULL;
  double **c_matrix = NULL, **mue_matrix = NULL, **u_matrix = NULL;
  int *fix_vektor = NULL;
  ghmm_density_t  density = 0;

  m_read = n_read = pi_read = a_read = c_read = mue_read = u_read
    = density_read = prior_read = fix_read = cos_read = 0;
  *multip = 1;                  /* default */
  ARRAY_CALLOC (smo, 1);

  ighmm_scanner_consume (s, '{');
  if (s->err)
    goto STOP;
  while (!s->err && !s->eof && s->c - '}') {
    ighmm_scanner_get_name (s);
    if (strcmp (s->id, "M") && strcmp (s->id, "N") && strcmp (s->id, "Pi") &&
        strcmp (s->id, "A") && strcmp (s->id, "C") && strcmp (s->id, "Mue") &&
        strcmp (s->id, "U") && strcmp (s->id, "multip")
        && strcmp (s->id, "cos") && strcmp (s->id, "prior")
        && strcmp (s->id, "fix_state") && strcmp (s->id, "density")
        && strncmp (s->id, "Ak_", 3)) {
      ighmm_scanner_error (s, "unknown identifier");
      goto STOP;
    }
    ighmm_scanner_consume (s, '=');
    if (s->err)
      goto STOP;
    if (!strcmp (s->id, "multip")) {
      *multip = ighmm_scanner_get_int (s);
      if (*multip < 1) {        /* ignore: makes no sense */
        *multip = 1;
        mes_prot ("Multip < 1 ignored\n");
      }
    }
    else if (!strcmp (s->id, "M")) {    /* number of output components */
      if (m_read) {
        ighmm_scanner_error (s, "identifier M twice");
        goto STOP;
      }
      M = ighmm_scanner_get_int (s);
      m_read = 1;
    }
    else if (!strcmp (s->id, "N")) {    /* number of states */
      if (n_read) {
        ighmm_scanner_error (s, "identifier N twice");
        goto STOP;
      }
      smo->N = ighmm_scanner_get_int (s);
      ARRAY_CALLOC (smo->s, smo->N);
      /*ARRAY_CALLOC (smo->density, smo->N);*/
      n_read = 1;
    }
    else if (!strcmp (s->id, "density")) {      /* which density function? */
      if (density_read) {
        ighmm_scanner_error (s, "identifier density twice");
        goto STOP;
      }
      /*current smo format only suport models with same densities*/      
      density = (ghmm_density_t) ighmm_scanner_get_int (s);
      if ((int) density < 0 || density >= density_number) {
        ighmm_scanner_error (s, "unknown typ of density function");
        goto STOP;
      }
      density_read = 1;
    }
    else if (!strcmp (s->id, "prior")) {        /* modelprior */
      if (prior_read) {
        ighmm_scanner_error (s, "identifier prior twice");
        goto STOP;
      }
      smo->prior = ighmm_scanner_get_edouble (s);
      if ((smo->prior < 0 || smo->prior > 1) && smo->prior != -1) {
        ighmm_scanner_error (s, "invalid model prior");
        goto STOP;
      }
      prior_read = 1;
    }
    else if (!strcmp (s->id, "cos")) {  /* number of transition classes */
      if (cos_read) {
        ighmm_scanner_error (s, "identifier cos twice");
        goto STOP;
      }
      smo->cos = ighmm_scanner_get_int (s);
      if (smo->cos <= 0) {
        ighmm_scanner_error (s, "invalid model cos");
        goto STOP;
      }
      cos_read = 1;
    }
    else if (!strcmp (s->id, "Pi")) {   /* Initial State Prob. */
      if (!n_read) {
        ighmm_scanner_error (s, "need N as a range for Pi");
        goto STOP;
      }
      if (pi_read) {
        ighmm_scanner_error (s, "identifier Pi twice");
        goto STOP;
      }
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "vector")) {
        ighmm_scanner_consume (s, '{');
        if (s->err)
          goto STOP;
        pi_vektor = scanner_get_double_earray (s, &len);
        if (len != smo->N) {
          ighmm_scanner_error (s, "wrong number of elements in PI");
          goto STOP;
        }
        ighmm_scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        ighmm_scanner_consume (s, '}');
        if (s->err)
          goto STOP;
        pi_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "fix_state")) {
      if (!n_read) {
        ighmm_scanner_error (s, "need N as a range for fix_state");
        goto STOP;
      }
      if (fix_read) {
        ighmm_scanner_error (s, "identifier fix_state twice");
        goto STOP;
      }
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "vector")) {
        ighmm_scanner_consume (s, '{');
        if (s->err)
          goto STOP;
        fix_vektor = scanner_get_int_array (s, &len);
        if (len != smo->N) {
          ighmm_scanner_error (s, "wrong number of elements in fix_state");
          goto STOP;
        }
        ighmm_scanner_consume (s, ';');
        if (s->err)
          goto STOP;
        ighmm_scanner_consume (s, '}');
        if (s->err)
          goto STOP;
        fix_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    /* 1. Possibility: one matrix for all transition classes */
    else if (!strcmp (s->id, "A")) {
      if (!cos_read) {
        ighmm_scanner_error (s, "cos unknown (needed for dim(A))");
        goto STOP;
      }
      if (!n_read) {
        ighmm_scanner_error (s, "need N as a range for A");
        goto STOP;
      }
      if (a_read) {
        ighmm_scanner_error (s, "Identifier A twice");
        goto STOP;
      }
      ARRAY_CALLOC (a_matrix, smo->N);
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (ighmm_cmatrix_read (s, a_matrix, smo->N, smo->N)) {
          ighmm_scanner_error (s, "unable to read matrix A");
          goto STOP;
        }
        a_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
      /* copy transition matrix to all transition classes */
      ARRAY_CALLOC (a_3dmatrix, smo->cos);
      a_3dmatrix[0] = a_matrix;
      for (i = 1; i < smo->cos; i++) {
        a_3dmatrix[i] = ighmm_cmatrix_alloc_copy (smo->N, smo->N, a_matrix);
        if (!a_3dmatrix[i]) {
          mes_proc ();
          goto STOP;
        }
      }
    }
    /* 2. Possibility: one matrix for each transition class specified */
    else if (!strncmp (s->id, "Ak_", 3)) {
      if (!cos_read) {
        ighmm_scanner_error (s, "cos unknown (needed for dim(A))");
        goto STOP;
      }
      if (!n_read) {
        ighmm_scanner_error (s, "need N as a range for A");
        goto STOP;
      }
      if (a_read) {
        ighmm_scanner_error (s, "identifier A twice");
        goto STOP;
      }
      ARRAY_CALLOC (a_3dmatrix, smo->cos);
      for (osc = 0; osc < smo->cos; osc++) {
        ARRAY_CALLOC (a_3dmatrix[osc], smo->N);
        ighmm_scanner_get_name (s);
        if (!strcmp (s->id, "matrix")) {
          if (ighmm_cmatrix_read (s, a_3dmatrix[osc], smo->N, smo->N)) {
            ighmm_scanner_error (s, "unable to read matrix Ak");
            goto STOP;
          }
        }
        else {
          ighmm_scanner_error (s, "unknown identifier");
          goto STOP;
        }
        if (osc < smo->cos - 1) {
          ighmm_scanner_consume (s, ';');
          if (s->err)
            goto STOP;
          /* read next matrix */
          ighmm_scanner_get_name (s);
          if (strncmp (s->id, "Ak_", 3)) {
            ighmm_scanner_error (s, "not enough matrices Ak");
            goto STOP;
          }
          ighmm_scanner_consume (s, '=');
          if (s->err)
            goto STOP;
        }
      }
      a_read = 1;
    }
    else if (!strcmp (s->id, "C")) {    /* weight for output components */
      if ((!n_read) || (!m_read)) {
        ighmm_scanner_error (s, "need M and N as a range for C");
        goto STOP;
      }
      if (c_read) {
        ighmm_scanner_error (s, "identifier C twice");
        goto STOP;
      }
      ARRAY_CALLOC (c_matrix, smo->N);
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (ighmm_cmatrix_read (s, c_matrix, smo->N, M)) {
          ighmm_scanner_error (s, "unable to read matrix C");
          goto STOP;
        }
        c_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "Mue")) {  /* mean of normal density */
      if ((!n_read) || (!m_read)) {
        ighmm_scanner_error (s, "need M and N as a range for Mue");
        goto STOP;
      }
      if (mue_read) {
        ighmm_scanner_error (s, "identifier Mue twice");
        goto STOP;
      }
      ARRAY_CALLOC (mue_matrix, smo->N);
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (ighmm_cmatrix_read (s, mue_matrix, smo->N, M)) {
          ighmm_scanner_error (s, "unable to read matrix Mue");
          goto STOP;
        }
        mue_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else if (!strcmp (s->id, "U")) {    /* variances of normal density */
      if ((!n_read) || (!m_read)) {
        ighmm_scanner_error (s, "need M and N as a range for U");
        goto STOP;
      }
      if (u_read) {
        ighmm_scanner_error (s, "identifier U twice");
        goto STOP;
      }
      ARRAY_CALLOC (u_matrix, smo->N);
      ighmm_scanner_get_name (s);
      if (!strcmp (s->id, "matrix")) {
        if (ighmm_cmatrix_read (s, u_matrix, smo->N, M)) {
          ighmm_scanner_error (s, "unable to read matrix U");
          goto STOP;
        }
        u_read = 1;
      }
      else {
        ighmm_scanner_error (s, "unknown identifier");
        goto STOP;
      }
    }
    else {
      ighmm_scanner_error (s, "unknown identifier");
      goto STOP;
    }
    ighmm_scanner_consume (s, ';');
    if (s->err)
      goto STOP;
  }                             /* while(..s->c-'}') */
  ighmm_scanner_consume (s, '}');
  if (s->err)
    goto STOP;

  if (m_read == 0 || n_read == 0 || pi_read == 0 || a_read == 0 ||
      c_read == 0 || mue_read == 0 || u_read == 0 || density_read == 0) {
    ighmm_scanner_error (s,
                   "some identifier not specified (N, M, Pi, A, C, Mue, U or density)");
    goto STOP;
  }
  /* set prior to -1 it none was specified */
  if (prior_read == 0)
    smo->prior = -1;
  /* default for fix is 0 */
  if (fix_read == 0) {
    ARRAY_CALLOC (fix_vektor, smo->N);
    for (i = 0; i < smo->N; i++)
      fix_vektor[i] = 0;
  }
  /* memory alloc for all transition matrices. If a transition is possible in one
     class --> alloc memory for all classes */
  for (i = 0; i < smo->N; i++) {
    /*smo->density[i] = density;*/
    for (j = 0; j < smo->N; j++) {
      out = in = 0;
      for (osc = 0; osc < smo->cos; osc++) {
        if (a_3dmatrix[osc][i][j] > 0.0)
          out = 1;
        if (a_3dmatrix[osc][j][i] > 0.0)
          in = 1;
      }
      smo->s[i].out_states += out;
      smo->s[i].in_states += in;
    }
    if (ghmm_c_state_alloc (smo->s + i, M, smo->s[i].in_states,
                            smo->s[i].out_states, smo->cos)) {
      mes_proc ();
      goto STOP;
    }
    /* given no support of the smo file, make all densities of same type and
       the cutoff value for the tail gaussian equal to the previously 
       used constant */
    smo->M = M;
    for (j=0; j < smo->M; j++){
      smo->s[i].density[j] = (ghmm_density_t)density;
      smo->s[i].a[j] = -GHMM_EPS_NDT;
      smo->s[i].M = M;
    }
    /* copy values read to smodel */
    if (smodel_copy_vectors
        (smo, i, pi_vektor, fix_vektor, a_3dmatrix, c_matrix, mue_matrix,
         u_matrix)) {
      mes_proc ();
      goto STOP;
    }

  }
  if (a_3dmatrix)
    for (i = 0; i < smo->cos; i++)
      ighmm_cmatrix_free (&(a_3dmatrix[i]), smo->N);
  m_free (a_3dmatrix);
  ighmm_cmatrix_free (&c_matrix, smo->N);
  ighmm_cmatrix_free (&mue_matrix, smo->N);
  ighmm_cmatrix_free (&u_matrix, smo->N);
  m_free (pi_vektor);
  m_free (fix_vektor);
  return (smo);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  if (a_3dmatrix)
    for (i = 0; i < smo->cos; i++)
      ighmm_cmatrix_free (&(a_3dmatrix[i]), smo->N);
  m_free (a_3dmatrix);
  ighmm_cmatrix_free (&c_matrix, smo->N);
  ighmm_cmatrix_free (&mue_matrix, smo->N);
  ighmm_cmatrix_free (&u_matrix, smo->N);
  m_free (pi_vektor);
  m_free (fix_vektor);
  ghmm_c_free (&smo);
  return NULL;
#undef CUR_PROC
}                               /* ghmm_c_read_block */

#endif /* GHMM_OBSOLETE */


/*============================================================================*/
/*============================================================================*/
ghmm_cmodel **ghmm_c_xml_read (const char *filename, int *smo_number)
{
#define CUR_PROC "ghmm_c_read"
  ghmm_fileData_s *f;
  ghmm_cmodel * * smo;

  f = ghmm_parseHMMDocument(filename);
  if (f) {
    assert(f->modelType &  GHMM_kContinuousHMM);
    *smo_number = (f->noModels);
    smo = f->model.c;
    
    free(f); /* XXX - by now, we free f */  
    return smo;
  } else
    return NULL;
#undef CUR_PROC
}                               /* ghmm_c_read_block */


/*============================================================================*/
/*============================================================================*/
int ghmm_c_xml_write (const char *file, ghmm_cmodel ** smo, int smo_number)
{
#define CUR_PROC "ghmm_c_read"
  ghmm_fileData_s *f;

  ARRAY_MALLOC(f,1);
  f->noModels = smo_number;
  f->modelType = GHMM_kContinuousHMM;
  f->model.c = smo;
  ighmm_writeHMMDocument(f, file);
  free(f);  
  return 0;
STOP:  
  return -1;
#undef CUR_PROC
}                               /* ghmm_c_read_block */

/*============================================================================*/
int ghmm_c_free (ghmm_cmodel ** smo)
{
#define CUR_PROC "ghmm_c_free"
  int i;
  mes_check_ptr (smo, return (-1));
  for (i = 0; i < (*smo)->N; i++) {
    
    /* if there are no out_states field was never allocated */ 
    if ((*smo)->s[i].out_states > 0){
      m_free ((*smo)->s[i].out_id);
    }  
    /* if there are no in_states field was never allocated */ 
    if ((*smo)->s[i].in_states > 0){
      m_free ((*smo)->s[i].in_id);
    }  
    ighmm_cmatrix_free (&((*smo)->s[i].out_a), (*smo)->cos);
    ighmm_cmatrix_free (&((*smo)->s[i].in_a), (*smo)->cos);
    m_free ((*smo)->s[i].c);
    m_free ((*smo)->s[i].mue);
    m_free ((*smo)->s[i].u);
    m_free ((*smo)->s[i].a);
    m_free ((*smo)->s[i].mixture_fix);
    m_free ((*smo)->s[i].density);
    
  }
  m_free ((*smo)->s);

  if ((*smo)->class_change) {
    if ((*smo)->class_change->user_data) {
      m_free ((*smo)->class_change->user_data);
    }
    m_free ((*smo)->class_change);
  }
  /*m_free ((*smo)->density);*/
  m_free (*smo);
  return (0);
#undef CUR_PROC
}                               /* ghmm_c_free */


/*============================================================================*/
ghmm_cmodel *ghmm_c_copy (const ghmm_cmodel * smo)
{
# define CUR_PROC "ghmm_c_copy"
  int i, k, j, nachf, vorg, m;
  ghmm_cmodel *sm2 = NULL;
  ARRAY_CALLOC (sm2, 1);
  ARRAY_CALLOC (sm2->s, smo->N);
  for (i = 0; i < smo->N; i++) {
    nachf = smo->s[i].out_states;
    vorg = smo->s[i].in_states;
    ARRAY_CALLOC (sm2->s[i].out_id, nachf);
    sm2->s[i].out_a = ighmm_cmatrix_alloc (smo->cos, nachf);
    if (!sm2->s[i].out_a) {
      mes_proc ();
      goto STOP;
    }
    ARRAY_CALLOC (sm2->s[i].in_id, vorg);
    sm2->s[i].in_a = ighmm_cmatrix_alloc (smo->cos, vorg);
    if (!sm2->s[i].in_a) {
      mes_proc ();
      goto STOP;
    }
    ARRAY_CALLOC (sm2->s[i].c, smo->s[i].M);
    ARRAY_CALLOC (sm2->s[i].mue, smo->s[i].M);
    ARRAY_CALLOC (sm2->s[i].u, smo->s[i].M);
    ARRAY_CALLOC (sm2->s[i].mixture_fix, smo->s[i].M);
    ARRAY_CALLOC (sm2->s[i].a, smo->s[i].M);
    ARRAY_CALLOC (sm2->s[i].density, smo->s[i].M);
    
    /* copy values */
    for (j = 0; j < nachf; j++) {
      for (k = 0; k < smo->cos; k++)
        sm2->s[i].out_a[k][j] = smo->s[i].out_a[k][j];
      sm2->s[i].out_id[j] = smo->s[i].out_id[j];
    }
    for (j = 0; j < vorg; j++) {
      for (k = 0; k < smo->cos; k++)
        sm2->s[i].in_a[k][j] = smo->s[i].in_a[k][j];
      sm2->s[i].in_id[j] = smo->s[i].in_id[j];
    }
    for (m = 0; m < smo->s[i].M; m++) {
      sm2->s[i].c[m] = smo->s[i].c[m];
      sm2->s[i].mue[m] = smo->s[i].mue[m];
      sm2->s[i].u[m] = smo->s[i].u[m];
      sm2->s[i].mixture_fix[m] = smo->s[i].mixture_fix[m];
      sm2->s[i].a[m] = smo->s[i].a[m];
      sm2->s[i].density[m] = smo->s[i].density[m];
    }
    sm2->s[i].M = smo->s[i].M;
    sm2->s[i].pi = smo->s[i].pi;
    sm2->s[i].fix = smo->s[i].fix;
    sm2->s[i].out_states = nachf;
    sm2->s[i].in_states = vorg;
  }
  sm2->cos = smo->cos;
  sm2->N = smo->N;
  sm2->M = smo->M;
  sm2->prior = smo->prior;
  return (sm2);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  ghmm_c_free (&sm2);
  return (NULL);
# undef CUR_PROC
}                               /* ghmm_c_copy */


/*============================================================================*/
int ghmm_c_check (const ghmm_cmodel * smo)
{
# define CUR_PROC "ghmm_c_check"
  int valid = 0;
  double sum;
  int i, k, j;
  char * str;
  /* sum  Pi[i] == 1 ? */
  sum = 0.0;
  
  for (i = 0; i < smo->N; i++) {
    sum += smo->s[i].pi;
  }
  if (fabs (sum - 1.0) >= GHMM_EPS_PREC) {
    mes_prot ("sum Pi[i] != 1.0\n");
    valid = -1;
    /*goto STOP; */
  }
  /* only 0/1 in fix? */
  for (i = 0; i < smo->N; i++) {
    if (smo->s[i].fix != 0 && smo->s[i].fix != 1) {
      mes_prot ("in vector fix_state only 0/1 values!\n");
      valid = -1;
      /*goto STOP;*/
    }
  }
  for (i = 0; i < smo->N; i++) {
    if (smo->s[i].out_states == 0) {
      str = ighmm_mprintf (NULL, 0, "out_states = 0 (state %d -> final state!)\n", i);
      mes_prot (str);
      m_free (str);
    }
    /* sum  a[i][k][j] */
    for (k = 0; k < smo->cos; k++) {
      sum = 0.0;
      for (j = 0; j < smo->s[i].out_states; j++) {
        sum += smo->s[i].out_a[k][j];
      }
      if (fabs (sum - 1.0) >= GHMM_EPS_PREC) {
        str = ighmm_mprintf (NULL, 0, "sum out_a[j] = %.4f != 1.0 (state %d, class %d)\n", sum, i,k);
        mes_prot (str);
        m_free (str);
        valid = -1;
        /*goto STOP; */
      }
    }
    /* sum c[j] */
    sum = 0.0;
    for (j = 0; j < smo->s[i].M; j++)
      sum += smo->s[i].c[j];
    if (fabs (sum - 1.0) >= GHMM_EPS_PREC) {
      str = ighmm_mprintf (NULL, 0, "sum c[j] = %.2f != 1.0 (state %d)\n", sum, i);
      mes_prot (str);
      m_free (str);
      valid = -1;            
      /* goto STOP; */
    }
    /* check mue, u ? */
  }

/* for lazy-evaluation-like model checking
   uncomment 'goto STOP' statements and 'STOP:' line  */  
/* STOP: */
  return (valid);
# undef CUR_PROC
}                               /* ghmm_c_check */


/*============================================================================*/
int ghmm_c_check_compatibility (ghmm_cmodel ** smo, int smodel_number)
{
#define CUR_PROC "ghmm_c_check_compatibility"
/* XXX - old function not used any more !!! */
  int i, j;
  char * str;

  for (i = 0; i < smodel_number; i++)
    for (j = i + 1; j < smodel_number; j++) {
      if (smo[i]->N != smo[j]->N) {
        str = ighmm_mprintf (NULL, 0,
                "ERROR: different number of states in smodel %d (%d) and smodel %d (%d)",
                i, smo[i]->N, j, smo[j]->N);
        mes_prot (str);
        m_free (str);
        return (-1);
      }
      if (smo[i]->s[0].M != smo[j]->s[0].M) {
        str = ighmm_mprintf (NULL, 0,
                "ERROR: different number of possible outputs in smodel  %d (%d) and smodel %d (%d)",
                i, smo[i]->s[0].M, j, smo[j]->s[0].M);
        mes_prot (str);
        m_free (str);
        return (-1);
      }
    }
  return 0;
#undef CUR_PROC
}                               /* ghmm_c_check_compatibility */


/*============================================================================*/
double ghmm_c_get_random_var (ghmm_cmodel * smo, int state, int m)
{
# define CUR_PROC "ghmm_c_get_random_var"
  switch (smo->s[state].density[m]) {
  case normal_approx:
  case normal:
    return (ighmm_rand_normal (smo->s[state].mue[m], smo->s[state].u[m], 0));
  case normal_right:
    return (ighmm_rand_normal_right (smo->s[state].a[m],smo->s[state].mue[m], smo->s[state].u[m], 0));
  case normal_left:
    return -(ighmm_rand_normal_right (-smo->s[state].a[m],-smo->s[state].mue[m], smo->s[state].u[m], 0));
  case uniform:
	return (ighmm_rand_uniform_cont (0, smo->s[state].mue[m] ,smo->s[state].u[m]));
  default:
    ighmm_mes (MES_WIN, "Warning: density function not specified!\n");
    return (-1);
  }
# undef CUR_PROC
}                               /* ghmm_c_get_random_var */



/*============================================================================*/

ghmm_cseq *ghmm_c_generate_sequences (ghmm_cmodel * smo, int seed,
                                         int global_len, long seq_number,
                                         long label, int Tmax)
{
# define CUR_PROC "ghmm_c_generate_sequences"

  /* An end state is characterized by not having an output probabiliy. */

  ghmm_cseq *sq = NULL;
  int pos, n, i, j, m, reject_os, reject_tmax, badseq, class;
  double p, sum;
  int len = global_len, up = 0, stillbadseq = 0, reject_os_tmp = 0;

  sq = ghmm_cseq_calloc (seq_number);
  if (!sq) {
    mes_proc ();
    goto STOP;
  }

  /* A specific length of the sequences isn't given. As a model should have
     an end state, the konstant MAX_SEQ_LEN is used. */
  if (len <= 0)
    len = (int) GHMM_MAX_SEQ_LEN;

  /* Maximum length of a sequence not given */
  if (Tmax <= 0)
    Tmax = (int) GHMM_MAX_SEQ_LEN;


  /* rng is also used by ighmm_rand_std_normal 
     seed == -1: Initialization, has already been done outside the function */
  if (seed >= 0) {
    ghmm_rng_init ();
    if (seed > 0)
      GHMM_RNG_SET (RNG, seed);
    else                        /* Random initialization! */
      ghmm_rng_timeseed (RNG);
  }

  n = 0;
  reject_os = reject_tmax = 0;

  while (n < seq_number) {
    /* Test: A new seed for each sequence */
    /*   ghmm_rng_timeseed(RNG); */
    stillbadseq = badseq = 0;
    ARRAY_CALLOC (sq->seq[n], len);

    /* Get a random initial state i */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (i = 0; i < smo->N; i++) {
      sum += smo->s[i].pi;
      if (sum >= p)
        break;
    }
    if (i == smo->N) {          /* Can happen by a rounding error in the input */
      i--;
      while (i > 0 && smo->s[i].pi == 0.0)
        i--;
    }

    /* Get a random initial output
       -> get a random m and then respectively a pdf omega. */
    p = GHMM_RNG_UNIFORM (RNG);
    sum = 0.0;
    for (m = 0; m < smo->s[i].M; m++) {
      sum += smo->s[i].c[m];
      if (sum >= p)
        break;
    }
    if (m == smo->s[i].M)
      m--;
    /* Get random numbers according to the density function */
    sq->seq[n][0] = ghmm_c_get_random_var (smo, i, m);
    pos = 1;

    /* The first symbol chooses the start class */
    if (smo->cos == 1) {
      class = 0;
    }
    else {
      /*printf("1: cos = %d, k = %d, t = %d\n",smo->cos,n,state);*/

      if (!smo->class_change->get_class) {
        printf ("ERROR: get_class not initialized\n");
        return (NULL);
      }
      class = smo->class_change->get_class (smo, sq->seq[n], n, 0);
      if (class >= smo->cos){
        printf("ERROR: get_class returned index %d but model has only %d classes !\n",class,smo->cos);
        goto STOP;
      }
    }
    while (pos < len) {
      /* Get a new state */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (j = 0; j < smo->s[i].out_states; j++) {
        sum += smo->s[i].out_a[class][j];
        if (sum >= p)
          break;
      }

      if (j == smo->s[i].out_states) {  /* Can happen by a rounding error */
        j--;
        while (j > 0 && smo->s[i].out_a[class][j] == 0.0)
          j--;
      }
      if (sum == 0.0) {
        if (smo->s[i].out_states > 0) {
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
          else if (class < smo->cos - 1) {
            class++;
            up = 1;
            continue;
          }
          else {
            stillbadseq = 1;
            break;
          }
        }
        else {
          /* Final state reached, out of while-loop */
          break;
        }
      }

      i = smo->s[i].out_id[j];

      /* fprintf(stderr, "%d\n", i); */
      /*      fprintf(stderr, "%d\n", i); */

      /* Get output from state i */
      p = GHMM_RNG_UNIFORM (RNG);
      sum = 0.0;
      for (m = 0; m < smo->s[i].M; m++) {
        sum += smo->s[i].c[m];
        if (sum >= p)
          break;
      }

      if (m == smo->s[i].M) {
        m--;
        while (m > 0 && smo->s[i].c[m] == 0.0)
          m--;
      }
      /* Get a random number from the corresponding density function */
      sq->seq[n][pos] = ghmm_c_get_random_var (smo, i, m);

      /* Decide the class for the next step */
      if (smo->cos == 1) {
        class = 0;
      }
      else {
        /*printf("2: cos = %d, k = %d, t = %d\n",smo->cos,n,state);*/
        if (!smo->class_change->get_class) {
          printf ("ERROR: get_class not initialized\n");
          return (NULL);
        }
        class = smo->class_change->get_class (smo, sq->seq[n], n, pos);
        printf ("class = %d\n", class);
        if (class >= smo->cos){
          printf("ERROR: get_class returned index %d but model has only %d classes !\n",class,smo->cos);
          goto STOP;
        }
      }
      up = 0;
      pos++;

    }                           /* while (state < len) */
    if (badseq) {
      reject_os_tmp++;
    }

    if (stillbadseq) {
      reject_os++;
      m_free (sq->seq[n]);
      /*      printf("cl %d, s %d, %d\n", class, i, n); */
    }
    else if (pos > Tmax) {
      reject_tmax++;
      m_free (sq->seq[n]);
    }
    else {
      if (pos < len)

        ARRAY_REALLOC (sq->seq[n], pos);
      sq->seq_len[n] = pos;
#ifdef GHMM_OBSOLETE
      sq->seq_label[n] = label;
#endif /* GHMM_OBSOLETE */
      /* ighmm_cvector_print(stdout, sq->seq[n], sq->seq_len[n]," "," ",""); */
      n++;
    }
    /*    printf("reject_os %d, reject_tmax %d\n", reject_os, reject_tmax); */

    if (reject_os > 10000) {
      mes_prot ("Reached max. no. of rejections\n");
      break;
    }
    if (!(n % 1000))
      printf ("%d Seqs. generated\n", n);
  }                             /* n-loop */

  if (reject_os > 0)
    printf ("%d sequences rejected (os)!\n", reject_os);
  if (reject_os_tmp > 0)
    printf ("%d sequences changed class\n", reject_os_tmp - reject_os);
  if (reject_tmax > 0)
    printf ("%d sequences rejected (Tmax, %d)!\n", reject_tmax, Tmax);

  /* printf ("End of smodel_generate_sequences.\n"); */

  return (sq);
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  ghmm_cseq_free (&sq);
  return (NULL);
# undef CUR_PROC
}                               /* ghmm_c_generate_sequences */


/*============================================================================*/
int ghmm_c_likelihood (ghmm_cmodel * smo, ghmm_cseq * sqd, double *log_p)
{
# define CUR_PROC "ghmm_c_likelihood"
  int res = -1;
  double log_p_i;
  int matched, i;

  matched = 0;
  *log_p = 0.0;
  for (i = 0; i < sqd->seq_number; i++) {

    /* setting sequence number in class_change struct if necessary */
    if (smo->cos > 1) {
      if (!smo->class_change) {
        printf ("cos = %d but class_change not initialized !\n", smo->cos);
        goto STOP;
      }
      smo->class_change->k = i;
    }

    if (ghmm_c_logp (smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) {
      *log_p += log_p_i * sqd->seq_w[i];
      matched++;
    }
    else {
      /* Test: high costs for each unmatched Seq. */
      *log_p += GHMM_PENALTY_LOGP * sqd->seq_w[i];
      matched++;
      ighmm_mes (MES_WIN, "sequence[%d] can't be build.\n", i);
    }
  }
  if (!matched) {
    mes_prot ("NO sequence can be build.\n");
    goto STOP;
  }
  /* return number of "matched" sequences */
  res = matched;

  /* resetting sequence number to default value */
  if (smo->cos > 1) {
    smo->class_change->k = -1;
  }
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return (res);
# undef CUR_PROC
}                               /* ghmm_c_likelihood */

int ghmm_c_individual_likelihoods (ghmm_cmodel * smo, ghmm_cseq * sqd,
                                   double *log_ps)
{
  int matched, res;
  double log_p_i;
  int i;
  res = -1;
  matched = 0;

  for (i = 0; i < sqd->seq_number; i++) {

    /* setting sequence number in class_change struct if necessary */
    if (smo->cos > 1) {
      if (!smo->class_change) {
        printf ("cos = %d but class_change not initialized !\n", smo->cos);
        goto STOP;
      }
      smo->class_change->k = i;
    }
    if (ghmm_c_logp (smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) {
      log_ps[i] = log_p_i;
      matched++;
    }
    else {
      /* Test: very small log score for sequence cannot be produced. */
      log_ps[i] = -DBL_MAX;
      /* fprintf(stderr,"sequence[%d] cannot be build.\n", i); */
    }
  }

  res = matched;
  if (matched == 0) {
    fprintf (stderr, "smodel_likelihood: NO sequence can be build.\n");
  }

  /* resetting sequence number to default value */
  if (smo->cos > 1) {
    smo->class_change->k = -1;
  }

  /* return number of "matched" sequences */
  return res;
STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  return -1;
}

/*============================================================================*/
/* various print functions                                                    */
/*============================================================================*/


/*============================================================================*/
void ghmm_c_Ak_print (FILE * file, ghmm_cmodel * smo, int k, char *tab,
                      char *separator, char *ending)
{
  int i, j, out_state;
  for (i = 0; i < smo->N; i++) {
    out_state = 0;
    fprintf (file, "%s", tab);
    if (smo->s[i].out_states > 0 && smo->s[i].out_id[out_state] == 0) {
      fprintf (file, "%.4f", smo->s[i].out_a[k][out_state]);
      out_state++;
    }
    else
      fprintf (file, "0.0   ");
    for (j = 1; j < smo->N; j++) {
      if (smo->s[i].out_states > out_state &&
          smo->s[i].out_id[out_state] == j) {
        fprintf (file, "%s %.4f", separator, smo->s[i].out_a[k][out_state]);
        out_state++;
      }
      else
        fprintf (file, "%s 0.0   ", separator);
    }
    fprintf (file, "%s\n", ending);
  }
}                               /* ghmm_c_Ak_print */


/*============================================================================*/
void ghmm_c_C_print (FILE * file, ghmm_cmodel * smo, char *tab, char *separator,
                     char *ending)
{
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.4f", smo->s[i].c[0]);
    for (j = 1; j < smo->s[i].M; j++)
      fprintf (file, "%s %.4f", separator, smo->s[i].c[j]);
    fprintf (file, "%s\n", ending);
  }
}                               /* ghmm_c_C_print */


/*============================================================================*/
void ghmm_c_Mue_print (FILE * file, ghmm_cmodel * smo, char *tab, char *separator,
                       char *ending)
{
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.4f", smo->s[i].mue[0]);
    for (j = 1; j < smo->s[i].M; j++)
      fprintf (file, "%s %.4f", separator, smo->s[i].mue[j]);
    fprintf (file, "%s\n", ending);
  }
}                               /* ghmm_c_Mue_print */


/*============================================================================*/
void ghmm_c_U_print (FILE * file, ghmm_cmodel * smo, char *tab, char *separator,
                     char *ending)
{
  /* attention: choose precision big enough to allow printing of  
     EPS_U in const.h */
  int i, j;
  for (i = 0; i < smo->N; i++) {
    fprintf (file, "%s", tab);
    fprintf (file, "%.4f", smo->s[i].u[0]);
    for (j = 1; j < smo->s[i].M; j++)
      fprintf (file, "%s %.4f", separator, smo->s[i].u[j]);
    fprintf (file, "%s\n", ending);
  }
}                               /* ghmm_c_U_print */


/*============================================================================*/
void ghmm_c_Pi_print (FILE * file, ghmm_cmodel * smo, char *tab, char *separator,
                      char *ending)
{
  int i;
  fprintf (file, "%s%.4f", tab, smo->s[0].pi);
  for (i = 1; i < smo->N; i++)
    fprintf (file, "%s %.4f", separator, smo->s[i].pi);
  fprintf (file, "%s\n", ending);
}                               /* ghmm_c_Pi_print */

/*============================================================================*/
void ghmm_c_fix_print (FILE * file, ghmm_cmodel * smo, char *tab, char *separator,
                       char *ending)
{
  int i;
  fprintf (file, "%s%d", tab, smo->s[0].fix);
  for (i = 1; i < smo->N; i++)
    fprintf (file, "%s %d", separator, smo->s[i].fix);
  fprintf (file, "%s\n", ending);
}


/*============================================================================*/
void ghmm_c_print (FILE * file, ghmm_cmodel * smo)
{
/* old function - No support to heterogeneous densities */
  int k;
  fprintf (file,
           "SHMM = {\n\tM = %d;\n\tN = %d;\n\tdensity = %d;\n\tcos = %d;\n",
           smo->M, smo->N, (int) smo->s[0].density[0], smo->cos);
  /* smo files support only models with a single density */
  fprintf (file, "\tprior = %.5f;\n", smo->prior);
  fprintf (file, "\tPi = vector {\n");
  ghmm_c_Pi_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  fprintf (file, "\tfix_state = vector {\n");
  ghmm_c_fix_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  for (k = 0; k < smo->cos; k++) {
    fprintf (file, "\tAk_%d = matrix {\n", k);
    ghmm_c_Ak_print (file, smo, k, "\t", ",", ";");
    fprintf (file, "\t};\n");
  }
  fprintf (file, "\tC = matrix {\n");
  ghmm_c_C_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n\tMue = matrix {\n");
  ghmm_c_Mue_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n\tU = matrix {\n");
  ghmm_c_U_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  fprintf (file, "};\n\n");
}                               /* ghmm_c_print */

/*============================================================================*/
/* needed for hmm_input: only one A (=Ak_1 = Ak_2...) is written */
void ghmm_c_print_oneA (FILE * file, ghmm_cmodel * smo)
{
/* old function - No support to heterogeneous densities */
  fprintf (file,
           "SHMM = {\n\tM = %d;\n\tN = %d;\n\tdensity = %d;\n\tcos = %d;\n",
           smo->M, smo->N, (int) smo->s[0].density[0], smo->cos);
  fprintf (file, "\tprior = %.3f;\n", smo->prior);
  fprintf (file, "\tPi = vector {\n");
  ghmm_c_Pi_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  fprintf (file, "\tfix_state = vector {\n");
  ghmm_c_fix_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  /* Attention: assumption is: A_k are all the same */
  fprintf (file, "\tA = matrix {\n");
  ghmm_c_Ak_print (file, smo, 0, "\t", ",", ";");
  fprintf (file, "\t};\n");
  /***/
  fprintf (file, "\tC = matrix {\n");
  ghmm_c_C_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n\tMue = matrix {\n");
  ghmm_c_Mue_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n\tU = matrix {\n");
  ghmm_c_U_print (file, smo, "\t", ",", ";");
  fprintf (file, "\t};\n");
  fprintf (file, "};\n\n");
}                               /* ghmm_c_print */


/*============================================================================*/
double ghmm_c_calc_cmbm (ghmm_cmodel * smo, int state, int m, double omega)
{
  double bm = 0.0;
  switch (smo->s[state].density[m]) {
  case normal:
    bm = ighmm_rand_normal_density (omega, smo->s[state].mue[m],
                                 smo->s[state].u[m]);
    break;
  case normal_right:
    bm = ighmm_rand_normal_density_trunc (omega, smo->s[state].mue[m],
                                     smo->s[state].u[m],smo->s[state].a[m]);
    break;
  case normal_left:
    bm = ighmm_rand_normal_density_trunc (-omega, -smo->s[state].mue[m],
                                     smo->s[state].u[m],-smo->s[state].a[m]);
    break;
  case normal_approx:
    bm = ighmm_rand_normal_density_approx (omega,
                                        smo->s[state].mue[m],
                                        smo->s[state].u[m]);
    break;
  case uniform:
    bm = ighmm_rand_uniform_density (omega, 
                                  smo->s[state].mue[m],
                                  smo->s[state].u[m]);
    break;
  default:
    ighmm_mes (MES_WIN, "Warning: density function not specified!\n");
  }
  if (bm == -1) {
    ighmm_mes (MES_WIN, "Warning: density function returns -1!\n");
    bm = 0.0;
  }
  return (smo->s[state].c[m] * bm);
}                               /* ghmm_c_calc_cmbm */


/*============================================================================*/
/* PDF(omega) in a given state */
double ghmm_c_calc_b (ghmm_cmodel * smo, int state, double omega)
{
  int m;
  double b = 0.0;
  for (m = 0; m < smo->s[state].M; m++)
    b += ghmm_c_calc_cmbm (smo, state, m, omega);
  return (b);
}                               /* ghmm_c_calc_b */


/*============================================================================*/
double ghmm_c_prob_distance (ghmm_cmodel * cm0, ghmm_cmodel * cm, int maxT,
                             int symmetric, int verbose)
{
#define CUR_PROC "ghmm_c_prob_distance"

#define STEPS 100

  double p0, p;
  double d = 0.0;
  double *d1;
  ghmm_cseq *seq0 = NULL;
  ghmm_cseq *tmp = NULL;
  ghmm_cmodel *smo1, *smo2;
  int i, t, a, k, n;
  int true_len;
  int true_number;
  int left_to_right = 0;
  long total, index;
  int step_width = 0;
  int steps = 1;

  if (verbose) {                /* If we are doing it verbosely we want to have STEPS steps */
    step_width = maxT / STEPS;
    steps = STEPS;
  }
  else                          /* else just one */
    step_width = maxT;

  ARRAY_CALLOC (d1, steps);

  smo1 = cm0;
  smo2 = cm;

  for (k = 0; k < 2; k++) {

    seq0 = ghmm_c_generate_sequences (smo1, 0, maxT + 1, 1, 0, maxT + 1);

    /*ghmm_cseq_print(stdout,seq0,0);*/

    if (seq0->seq_len[0] < maxT) {      /* There is an absorbing state */

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
        tmp = ghmm_c_generate_sequences (smo1, 0, 0, a, 0, maxT + 1);
        ghmm_cseq_add (seq0, tmp);
        ghmm_cseq_free (&tmp);

        total = 0;
        for (i = 0; i < seq0->seq_number; i++)
          total += seq0->seq_len[i];
      }
    }

    if (left_to_right) {

      for (t = step_width, i = 0; t <= maxT; t += step_width, i++) {

        index = 0;
        total = seq0->seq_len[0];

        /* Determine how many sequences we need to get a total of t
           and adjust length of last sequence to obtain total of 
           exactly t */

        while (total < t) {
          index++;
          total += seq0->seq_len[index];
        }

        true_len = seq0->seq_len[index];
        true_number = seq0->seq_number;

        if ((total - t) > 0)
          seq0->seq_len[index] = total - t;
        seq0->seq_number = index + 1;

        if (ghmm_c_likelihood (smo1, seq0, &p0) == -1) {
          /* error! */
          mes_prot ("seq0 can't be build from smo1!");
          goto STOP;


        }
        n = ghmm_c_likelihood (smo2, seq0, &p); /* ==-1: KEINE Seq. erzeugbar */
        if (n < seq0->seq_number) {
          ighmm_mes (MES_WIN,
               "problem: some seqences in seq0 can't be build from smo2\n");
          /* what shall we do now? */
          goto STOP;
        }

        /*printf("1/%d *(%f - %f)\n",t,p0,p);*/

        d = 1.0 / t * (p0 - p);

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
          printf ("%d\t%f\t%f\t%f\n", t, p0, p, d);

        seq0->seq_len[index] = true_len;
        seq0->seq_number = true_number;
      }
    }

    else {                      /* no left to right model */

      true_len = seq0->seq_len[0];

      for (t = step_width, i = 0; t <= maxT; t += step_width, i++) {
        seq0->seq_len[0] = t;

        if (ghmm_c_likelihood (smo1, seq0, &p0) == -1) {
          /* error case */
          mes_prot ("seq0 can't be build from smo1!");
          goto STOP;
        }
        n = ghmm_c_likelihood (smo2, seq0, &p);/*== -1: KEINE Seq. erzeugbar*/
        if (n < seq0->seq_number) {
          ighmm_mes (MES_WIN,
               "problem: some sequences in seq0 can't be build from smo2\n");
          /* what shall we do now? */
          goto STOP;
        }

        /*printf("1/%d *(%f - %f)\n",t,p0,p);*/

        d = 1.0 / t * (p0 - p);

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
          printf ("%d\t%f\t%f\t%f\n", t, p0, p, d);

      }
      seq0->seq_len[0] = true_len;
    }

    if (symmetric) {
      ghmm_cseq_free (&seq0);
      smo1 = cm;
      smo2 = cm0;
    }
    else
      break;

  }                             /* k = 1,2 */

  ghmm_cseq_free (&seq0);

  return d;

STOP:     /* Label STOP from ARRAY_[CM]ALLOC */
  ghmm_cseq_free (&seq0);
  return (-1.0);
#undef CUR_PROC
}


/*============================================================================*/
double ghmm_c_calc_cmBm (ghmm_cmodel * smo, int state, int m, double omega)
{
  double Bm = 0.0;
  switch (smo->s[state].density[m]) {
  case normal_approx:
  case normal:
    Bm = ighmm_rand_normal_cdf (omega, smo->s[state].mue[m], smo->s[state].u[m]);
    break;
  case normal_right:
    Bm = ighmm_rand_normal_right_cdf (omega, smo->s[state].mue[m],
                              smo->s[state].u[m],smo->s[state].a[m]);
    break;
  case normal_left:
    Bm = (ighmm_rand_normal_right_cdf (omega, -smo->s[state].mue[m], smo->s[state].u[m],-smo->s[state].a[m]));
    break;
  case uniform:
    Bm = (ighmm_rand_uniform_cdf (omega, smo->s[state].mue[m], smo->s[state].u[m]));
    break;
  default:
    ighmm_mes (MES_WIN, "Warning: density function not specified!\n");
  }
  if (Bm == -1) {
    ighmm_mes (MES_WIN, "Warning: density function returns -1!\n");
    Bm = 0.0;
  }
  return (smo->s[state].c[m] * Bm);
}                               /* smodel_calc_Bm */


/*============================================================================*/
/* CDF(omega) in a given state */
double ghmm_c_calc_B (ghmm_cmodel * smo, int state, double omega)
{
  int m;
  double B = 0.0;
  for (m = 0; m < smo->s[state].M; m++)
    B += ghmm_c_calc_cmBm (smo, state, m, omega);
  return (B);
}                               /* ghmm_c_calc_B */

/*============================================================================*/
/* What is the dimension of the modell ( = dimension of the parameter vektor) ?
   count the number of free parameters in a field of models; used for calc. BIC
   Only those parameters, that can be changed during  training.
   mixture coeff from smix and priors are not counted! 
*/
int ghmm_c_count_free_parameter (ghmm_cmodel ** smo, int smo_number)
{
  int i, k;
  int pi_counted = 0, cnt = 0;

  for (k = 0; k < smo_number; k++) {
    pi_counted = 0;
    /* for states */
    for (i = 0; i < smo[k]->N; i++) {
      if (smo[k]->s[i].out_states > 1)
        /* multipl. with COS correct ??? */
        cnt += smo[k]->cos * (smo[k]->s[i].out_states - 1);
      if (smo[k]->s[i].pi != 0 && smo[k]->s[i].pi != 1) {
        pi_counted = 1;
        cnt++;
      }
      if (!smo[k]->s[i].fix) {
        if (smo[k]->s[i].M == 1)
          cnt += 2;             /* mu, sigma */
        else
          cnt += (3 * smo[k]->s[i].M);       /* c, mu, sigma */
      }
    }                           /* for (i ..) */
    if (pi_counted)
      cnt--;                    /* due to condition: sum(pi) = 1 */
    if (smo[k]->s[0].M > 1)
      cnt--;                    /* due to condition: sum(c) = 1 */
  }

  return cnt;
}

/*============================================================================*/
/* interval (a,b) with ~ B(a) < 0.01, B(b) > 0.99 */
void ghmm_c_get_interval_B (ghmm_cmodel * smo, int state, double *a, double *b)
{
  int m;
  double mue, delta, max, min;
  *a = DBL_MAX;
  *b = -DBL_MAX;
  for (m = 0; m < smo->s[state].M; m++) {
    switch (smo->s[state].density[m]) {
    case normal:
    case normal_approx:
      mue = smo->s[state].mue[m];
      delta = 3 * sqrt (smo->s[state].u[m]);
      if (*a > mue - delta)
        *a = floor (mue - delta);
      if (*b < mue + delta)
        *b = ceil (mue + delta);
    case normal_right:
    case normal_left:
      if (smo->s[state].density[state] == normal_right && *a < smo->s[state].a[m])
        *a = smo->s[state].a[m];
      if (smo->s[state].density[state] == normal_left && *b > smo->s[state].a[m])
        *b = smo->s[state].a[m] ;
    break;
    case uniform:
      max = smo->s[state].mue[m];
      min = smo->s[state].u[m];
      *a = floor ((0.01 * (max - min))+min);
      *b = ceil ((0.99 * (max - min))+min);
    break;      
    default:
      ighmm_mes (MES_WIN, "Warning: density function not specified!\n");
    }
  }

  return;
}                               /* ghmm_c_get_interval_B */

/*============================================================================*/
double ghmm_c_ifunc (ghmm_cmodel * smo, int state, double c, double x)
{
  return (fabs (ghmm_c_calc_B (smo, state, x) - c));
}

/*============================ TEST =========================================*/
#ifdef XXX
int smodel_test_callback(int pos){
   char* ModuleName = "class_change";
   char* FunctionName = "getClass";
   int class;
   PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *pValue;
    
   /*Py_Initialize();*/      /* Init Python Interpreter*/
   
   /*PyRun_SimpleString("import sys\n");*/
   /*PyRun_SimpleString("sys.stdout.write('Hello from an embedded Python Script\\n')\n"); */

   
   printf("C: Importing Python module ... ");
   pName = PyString_FromString(ModuleName);
   pModule = PyImport_Import(pName);       /* Import module*/
   pDict = PyModule_GetDict(pModule);
   printf("done.\n");    
    
   printf("C: Calling Python with value %d\n",pos);
   pFunc = PyDict_GetItemString(pDict, FunctionName);

   pArgs = PyTuple_New(1);
   pValue = PyInt_FromLong(pos);

   PyTuple_SetItem(pArgs, 0, pValue); 
   pValue = PyObject_CallObject(pFunc, pArgs); /* Calling Python */

   /* parsing the result from Python to C data type
   class = PyInt_AsLong(pValue);
   printf("C: The returned class is %d\n",class); */
     
   /*Py_Finalize();         */
   
   return class; 
 
}
#endif
