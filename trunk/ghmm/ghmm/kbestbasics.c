/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/kbestbasics.c
*       Authors:  Alexander Riemer, Janne Grunau
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


#include <stdlib.h>
#include <math.h>
#include "mes.h"
#include "model.h"
#include "kbestbasics.h"

/* inserts new hypothesis into list at position indicated by pointer plist */
inline void hlist_insertElem(hypoList** plist, int newhyp, hypoList* parlist) {
#define CUR_PROC "hlist_insertElem"
    hypoList* newlist;

    if (!m_calloc(newlist, 1)) {mes_proc(); goto STOP;}
    newlist->hyp_c = newhyp;
    if (parlist)
      parlist->refcount += 1;
    newlist->parent = parlist;
    newlist->next = *plist;

    *plist = newlist;
    return;
STOP:
    mes_prot("hlist_insertElem failed\n");
    exit(1);
#undef CUR_PROC
}

/* removes hypothesis at position indicated by pointer plist from the list
   removes recursively parent hypothesis with refcount==0 */
inline void hlist_removeElem(hypoList** plist) {
#define CUR_PROC "hlist_removeElem"
    hypoList* tempPtr = (*plist)->next;

    free((*plist)->gamma_a);
    free((*plist)->gamma_id);
    if ((*plist)->parent) {
      (*plist)->parent->refcount -= 1;
      if (0==(*plist)->parent->refcount)
	hlist_removeElem(&((*plist)->parent));
    }
    free(*plist);

    *plist = tempPtr;
#undef CUR_PROC
}

/**
   Propagates list of hypotheses forward by extending each old hypothesis to
   the possible new hypotheses depending on the states in which the old
   hypothesis could end and the reachable labels
   @return number of old hypotheses
   @param mo:         pointer to the model
   @param h:          pointer to list of hypotheses
   @param hplus:      address of a pointer to store the propagated hypotheses
   @param labels:     number of labels
   @param nr_s:       number states which have assigned a label 
 */
int hlist_propFwd(model* mo, hypoList* h, hypoList** hplus, int labels, int* nr_s) {
#define CUR_PROC "hlist_propFwd"
  int i, j, c, k;
  int i_id, j_id, g_nr;
  int no_oldHyps=0;
  hypoList* hP=h;
  hypoList** created;

  if (!m_malloc(created, labels)) {mes_proc(); goto STOP;}

  /* extend the all hypotheses with the labels of out_states of all states in the hypotesis */
  while (hP != NULL) {

    /* lookup table for labels, created[i]!=0 iff the current hypotheses was propagated forward with label i */
    for (c=0; c<labels; c++) created[c] = NULL;

    /* extend the current hypothesis and add all states which may have probability greater null */
    for (i=0; i<hP->gamma_states; i++) {
      i_id = hP->gamma_id[i];
      for (j=0; j<mo->s[i_id].out_states; j++) {
	j_id = mo->s[i_id].out_id[j];
	c = mo->s[j_id].label;

	/* create a new hypothesis with label c */
	if (!created[c]) {
	  hlist_insertElem(hplus, c, hP);
	  created[c] = *hplus;
	  /* initiallize gamma-array with safe size (number of states */
	  if (!m_malloc((*hplus)->gamma_a,  nr_s[c])) {mes_proc(); goto STOP;}
	  if (!m_malloc((*hplus)->gamma_id, nr_s[c])) {mes_proc(); goto STOP;}
	  (*hplus)->gamma_id[0]  = j_id;
	  (*hplus)->gamma_a[0]   = 1.0;
	  (*hplus)->gamma_states = 1; 
	}
	/* add a new gamma state to the existing hypothesis with c */
	else {
	  g_nr = created[c]->gamma_states;
	  /* search for state j_id in the gamma list */
	  for (k=0; k<g_nr; k++) 
	    if (j_id == created[c]->gamma_id[k])
	      break;
	  /* add the state to the gamma list */
	  if (k==g_nr) {
	    created[c]->gamma_id[g_nr] = j_id;
	    created[c]->gamma_a[g_nr]  = 1.0;
	    created[c]->gamma_states   = g_nr+1;
	  }
	}
      }
    }
    hP = hP->next;
    no_oldHyps++;
  }

  /* reallocating gamma-array to the correct size */
  for (c=0; c<labels; c++) {
    if (created[c]) {
      if (m_realloc(created[c]->gamma_a,  created[c]->gamma_states)) {mes_proc(); goto STOP;}
      if (m_realloc(created[c]->gamma_id, created[c]->gamma_states)) {mes_proc(); goto STOP;}
    }
  }

  free(created);
  return (no_oldHyps);
STOP:
  mes_prot("hlist_propFwd failed\n");
  exit(1);
#undef CUR_PROC
}


/**
   Calculates the logarithm of sum(exp(log_a[j,a_pos])+exp(log_gamma[j,g_pos]))
   which corresponds to the logarithm of the sum of a[j,a_pos]*gamma[j,g_pos]
   @return logSum for products of a row from gamma and a row from matrix A
   @param log_a:      row of the transition matrix with logarithmic values (1.0 for log(0))
   @param s:          state whose gamma-value is calculated
   @param parent:     a pointer to the parent hypothesis
*/
inline double logGammaSum(double* log_a, state* s, hypoList* parent) {
#define CUR_PROC "logGammaSum"
  double result;
  int j, j_id, k;
  double max=1.0;
  int argmax=0;
  double* logP;

  /* shortcut for the trivial case */
  if (parent->gamma_states==1)
    for (j=0; j<s->in_states; j++)
      if (parent->gamma_id[0]==s->in_id[j])
	return parent->gamma_a[0] + log_a[j];

  if (!m_malloc(logP, s->in_states)) {mes_proc(); goto STOP;}

  /* calculate logs of a[k,l]*gamma[k,hi] as sums of logs and find maximum: */
  for (j=0; j<s->in_states; j++) {
    j_id = s->in_id[j];
    /* search for state j_id in the gamma list */
    for (k=0; k<parent->gamma_states; k++)
      if (parent->gamma_id[k] == j_id)
	break;
    if (k==parent->gamma_states)
      logP[j] = 1.0;
    else {
      logP[j] = log_a[j] + parent->gamma_a[k];
      if (max==1.0 || logP[j]>max) {
	max = logP[j];
	argmax=j;
      }
    }
  }

  /* calculate max+log(1+sum[j!=argmax; exp(logP[j]-max)])  */
  result = 1.0;
  for (j=0; j<s->in_states; j++)
    if (j!=argmax && logP[j]!=1.0) result+= exp(logP[j]-max);
  
  result  = log(result);
  result += max;

  free(logP);
  return result;
STOP:
  mes_prot("logGammaSum failed\n");
  exit(1);
#undef CUR_PROC
}


/**
   Calculates the logarithm of the sum of the probabilities whose logarithms are
   stored in the given array
   @return log of sum of exp(a[i])
   @param a:          array of logarithms of probabilities (a[i] < 0 for all i)
   @param N:          length of a
*/
inline double logSum(double* a, int N) {
#define CUR_PROC "logSum"
  int i;
  double max=1.0;
  int argmax=0;
  double result;

  /* find maximum value in a: */
  for (i=0; i<N; i++)
    if (a[i]<KBEST_EPS && (max==1.0 || a[i]>max)) {
      max=a[i];
      argmax=i;
    }

  /* calculate max+log(1+sum[i!=argmax; exp(a[i]-max)])  */
  result = 1.0;
  for (i=0; i<N; i++) {
    if (a[i]<KBEST_EPS && i!= argmax) result+= exp(a[i]-max);
  }
  result=log(result);
  result+=max;
  return result;
#undef CUR_PROC
}
