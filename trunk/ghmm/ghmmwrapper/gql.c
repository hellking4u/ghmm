/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: gql.c
*       Authors:   Wasinee Rungsarityotin, Benjamin Georgi
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

/* XXX FIXME: breaks out of tree build of ghmmwrapper */
#include "../config.h"

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/mes.h>
#include <ghmm/matrix.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/smodel.h>
#include <ghmm/sfoba.h>
#include <ghmm/obsolete.h>
/*#include <ghmm++/GHMM_convertXMLtoC.h>*/
#include <float.h>
#include <assert.h>


static double *log_p_pt;
static int seq_rank(const void *a1, const void *a2);
  
int seq_rank(const void *a1, const void *a2) {
  
  int i, j;

  i = *((int*)a1);
  j = *((int*)a2);

  if ( log_p_pt[i] == log_p_pt[j] ) return 0;
  else
    if ( log_p_pt[i] < log_p_pt[j] ) return -1;
    else
      return 1;
}

static int ghmm_c_state_alloc(sstate *state,
			      int M,
			      int in_states,
			      int out_states,
			      int cos) {
#define CUR_PROC "ghmm_c_state_alloc"
  int res = -1;
  if (!((state->c) = ighmm_calloc(sizeof(*(state->c)) * M)))
    {mes_proc(); goto STOP;}
  if (!((state->mue) = ighmm_calloc(sizeof(*(state->mue)) * M)))
    {mes_proc(); goto STOP;}
  if (!((state->u) = ighmm_calloc(sizeof(*(state->u)) * M)))
    {mes_proc(); goto STOP;}
  if (!((state->a) = ighmm_calloc(sizeof(*(state->a)) * M)))
    {mes_proc(); goto STOP;}
  if (!((state->density) = ighmm_calloc(sizeof(*(state->density)) * M)))
    {mes_proc(); goto STOP;}
    
  if (out_states > 0) {
    if (!((state->out_id) = ighmm_calloc(sizeof(*(state->out_id)) * out_states)))
      {mes_proc(); goto STOP;}
    state->out_a = ighmm_cmatrix_alloc(cos, out_states);
    if(!state->out_a) {mes_proc(); goto STOP;}
  }
  if (in_states > 0) {
    if (!((state->in_id) = ighmm_calloc(sizeof(*(state->in_id)) * (in_states)))) 
      {mes_proc(); goto STOP;}
    state->in_a = ighmm_cmatrix_alloc(cos, in_states);
    if(!state->in_a) {mes_proc(); goto STOP;}
  }
  res = 0;
STOP:
  return(res);
#undef CUR_PROC
}

smodel *smodel_alloc_fill(int N, int M, int cos, double prior, int density) {
#define CUR_PROC "smodel_alloc_fill"
  int i;
  smodel *smo=NULL;
  if (!(smo = ighmm_malloc (sizeof (smodel)))) {mes_proc(); goto STOP;}  
  smo->M   = M;
  smo->N   = N;
  smo->cos = cos;
  smo->prior = prior;
  if (!(smo->s = ighmm_calloc(sizeof(*(smo->s)) * (smo->N))))
    {mes_proc(); goto STOP;}

  for(i=0; i < smo->N; i++) {
    ghmm_c_state_alloc(&smo->s[i], smo->M, smo->N, smo->N, cos);
  }
  return smo;
STOP:
  fprintf (stderr, "smodel_alloc_fill(1): out of memory\n");
  return NULL;
#undef CUR_PROC
}

void smodel_set_pivector(smodel *smo, int i, double piv) {
  if (smo->s != NULL) {
    smo->s[i].pi = piv;
  }
}

void smodel_set_fixvector(smodel *smo, int i, double fixv) {
  if (smo->s != NULL) {
    smo->s[i].fix = fixv;
  }
}

void smodel_set_transition(smodel *smo, int i, int j, int cos, double prob) {
  int in, out;
  if (cos >= smo->cos) {
    fprintf(stderr, "smodel_set_transition(cos): cos > state->cos\n");
    exit(-1);	
  }
  if (smo->s && smo->s[i].out_a && smo->s[j].in_a) {
    for(out=0; out < smo->s[i].out_states; out++) {
      if ( smo->s[i].out_id[out] == j ) {
	smo->s[i].out_a[cos][out] = prob;
	fprintf(stderr, "smodel_set_transition(0):State %d, %d, = %f\n", i, j, prob);
	break;
      }
    }

    for(in=0; in < smo->s[j].in_states; in++) {
      if ( smo->s[j].in_id[in] == i ) {
	smo->s[j].in_a[cos][in] = prob;
	break;
      }
    }
  }
}


double smodel_get_transition(smodel *smo, int i, int j, int cos) {
  int out;
  if (cos >= smo->cos) {
    fprintf(stderr, "smodel_get_transition(0): cos > state->cos\n");
    exit(-1);	
  }
  if (smo->s && smo->s[i].out_a && smo->s[j].in_a) {
    for(out=0; out < smo->s[i].out_states; out++) {
      if ( smo->s[i].out_id[out] == j ) {
	return 	smo->s[i].out_a[cos][out];
      }
    }
  }
  fprintf(stderr, "smodel_get_transition(1): data structure not initialized\n");
  return -1.0;
}

void smodel_set_mean(smodel *smo, int i, double *mu) {
  int m;
  if (smo->s != NULL) {
    for(m = 0; m < smo->M; m++)
      smo->s[i].mue[m] = mu[m];
  }
}

void smodel_set_variance(smodel *smo, int i, double *variance) {
  int m;
  if (smo->s != NULL) {
    for(m = 0; m < smo->M; m++) {
      smo->s[i].u[m] = variance[m];
      assert( smo->s[i].u[m] > 0.0 );
    }
  }
}

/*
void call_smodel_print(char *filename, smodel *smo) {
  FILE *fp=fopen(filename, "w");
  if (fp == NULL) {
    fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);    
  } else {
    ghmm_c_print(fp, smo);
    fclose(fp);
  }
} */

int smodel_sorted_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps, int *seq_rank) {
  int matched, res;
  double log_p_i;
  int i;
  res=-1;
  matched=0;
  for ( i = 0; i < sqd->seq_number; i++) {
    seq_rank[i] = i;
    if (ghmm_c_logp(smo, sqd->seq[i], sqd->seq_len[i], &log_p_i) != -1) { 
      log_ps[i] = log_p_i;
      matched++;
    }
    else  {
      /* Test: very small log score for sequence cannot be produced */
      log_ps[i] = -DBL_MAX;
      /*      ighmm_mes(MES_WIN, "sequence[%d] can't be build.\n", i); */
    }
  }

  res=matched;
  if (matched == 0) { 
    fprintf(stderr, "smodel_likelihood: NO sequence can be build.\n"); 
  } else  qsort( seq_rank, sqd->seq_number, sizeof(int), (int (*)(const void *, const void *))seq_rank);

  /* return number of "matched" sequences */
  return res;
}
