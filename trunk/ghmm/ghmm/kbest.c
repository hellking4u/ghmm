/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/kbest.c
*       Authors:  Anyess von Bock, Alexander Riemer, Janne Grunau
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
#include "model.h"
#include "kbest.h"
#include "kbestbasics.h"
#include <math.h>


/**
  Builds logarithmic transition matrix from the states' in_a values
  @return transition matrix with logarithmic values, 1.0 if a[i,j] = 0
  @param s:           array of all states of the model
  @param N:           number of states in the model
 */
inline double* buildLogMatrix(state* s, int N) {
  int i,j;
  double* log_a;      /* log(a(i,j)) => log_a[i*N+j] */

  /* create & initialize matrix: */
  log_a = (double*)malloc(sizeof(double)*N*N);
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      log_a[i*N+j]=1.0; /* positive values cannot occur as logarithms */

  for (i=0; i<N; i++)
    for (j=0; j<s[i].in_states; j++) {
      if (s[i].in_a[j] > KBEST_EPS)  /* only accept values greater than EPS */
	log_a[s[i].in_id[j]*N+i] = log(s[i].in_a[j]);
    }

  return log_a;
}


/**
   Calculates the most probable labeling for the given sequence in the given
   model using k-best decoding.
   Labels must be from interval [0:max_label] without gaps!!! (not checked)
   Model must not have silent states. (checked in Python wrapper)
   @return array of labels (internal representation)
   @param mo:         pointer to a model
   @param o_seq:      output sequence (array of internal representation chars)
   @param seq_len:    length of output sequence
   @param k:          number of hypotheses to keep for each state
   @param log_p:      variable reference to store the log prob. of the labeling
 */
int* kbest(model* mo, int* o_seq, int seq_len, int k, double* log_p) {
  int i, t, c, l, m;		/* counters */
  int no_oldHyps;		/* number of hypotheses until position t-1 */
  int b_index;			/* index for addressing states' b arrays */
  int no_labels=0;

  /* logarithmized transition matrix A, log(a(i,j)) => log_a[i*N+j],
       1.0 for zero probability */
  double* log_a;

  /* matrix of hypotheses, holds for every position in the sequence a list
     of hypotheses */
  hypoList **h;
  hypoList *hP;
  
  /* vectors for rows in the matrices */
  int *hypothesis;
  double* gammaRow;
  
  /* pointer & prob. of the k most probable hypotheses for each state
       - matrices of dimensions #states x k:  argm(i,l) => argmaxs[i*k+l] */
  double* maxima;
  hypoList** argmaxs;
  
  /* pointer to & probability of most probable hypothesis in a certain state */
  hypoList* argmax;
  double sum;

  /* break if sequence empty or k<1: */
  if (seq_len <= 0 || k<=0) return NULL;

  h = (hypoList**)calloc(seq_len, sizeof(hypoList*));

  /** 1. Initialization (extend empty hypothesis to #labels hypotheses of
         length 1): */

  /* get number of labels (= maximum label + 1): */
  for (i=0; i < mo->N; i++) {
    if (mo->s[i].label > no_labels)
      no_labels = mo->s[i].label;
  }
  no_labels++;

  /* initialize h: */
  h[0] = NULL;
  for (c=no_labels-1; c>=0; c--) {
    /* create #labels hypotheses of max. length seq_len: */
    hlist_insertElem(&(h[0]), c, NULL);
    /* create #hypotheses rows in gamma of length #states: */
    gammaRow = (double*)malloc(sizeof(double)*mo->N);
    for(i=0; i < mo->N; i++) {
      /* if hypothesis c ends with label of state i:
	   gamma(i,c):= log(pi[i]*b[i](o_seq[0]))
	 else: gamma(i,c):= -INF (represented by 1.0)
      */
      if (c == mo->s[i].label) {
      	b_index = get_emission_index(mo,i,o_seq[0],0);
      	if (b_index < 0 || mo->s[i].pi < KBEST_EPS || mo->s[i].b[b_index] < KBEST_EPS)
	  gammaRow[i] = 1.0;
        else
	  gammaRow[i] = log(mo->s[i].pi) + log(mo->s[i].b[b_index]);
      }
      else gammaRow[i] = 1.0;
    }
    /* save the gamma values in the hypotheses and
       chose all initial hypotheses */
    h[0]->gamma = gammaRow;
    h[0]->chosen = 1;
  }

  /* calculate transition matrix with logarithmic values: */
  log_a = buildLogMatrix(mo->s, mo->N);

  /* initialize temporary arrays: */
  maxima=(double*)malloc(sizeof(double)*mo->N*k);     /* for each state save k */
  argmaxs=(hypoList**)malloc(mo->N*k*sizeof(hypoList*));

  /*------ Main loop: Cycle through the sequence: ------*/
  for (t=1; t<seq_len; t++) {
    
    /* put o_seq[t-1] in emission history: */
    update_emission_history(mo,o_seq[t-1]);
    
    /** 2. Propagate hypotheses forward and update gamma: */
    no_oldHyps = propFwd(h[t-1], &(h[t]), no_labels, seq_len);
    
    /*-- calculate new gamma: --*/
    hP = h[t];
    /* cycle through list of hypotheses */
    while (hP != NULL) {
      hP->gamma = (double*)malloc(mo->N*sizeof(double));
      for (i=0; i < mo->N; i++) {
      	if (mo->s[i].label == hP->hyp_c) {
	  /* if hypothesis hP ends with label of state i:
	       gamma(i,c):= log(sum(exp(a(j,i)*exp(oldgamma(j,old_c)))))
	                      + log(b[i](o_seq[t]))
	     else: gamma(i,c):= -INF (represented by 1.0)
	  */
      	  hP->gamma[i] = logGammaSum(log_a, i, mo->N, hP->parent->gamma);
	  b_index=get_emission_index (mo, i, o_seq[t], t);
      	  if (b_index<0 || hP->gamma[i]==1.0 || mo->s[i].b[b_index]<KBEST_EPS)
	    hP->gamma[i] = 1.0;
	  else 
	    hP->gamma[i] += log(mo->s[i].b[b_index]);
	}
	else hP->gamma[i] = 1.0;
      }
      hP = hP->next;
    }
    
    /** 3. Choose the k most probable hypotheses for each state and discard all
	   hypotheses that were not chosen: */
    
    /* initialize temporary arrays: */
    for (i=0; i< mo->N*k; i++) {
      maxima[i]  = 1.0;
      argmaxs[i] = NULL;
    }
    
    /* cycle through hypotheses & calculate the k most probable hypotheses for
       each state: */
    hP=h[t];
    while (hP != NULL) {
      for (i=0; i < mo->N; i++) {
	if (hP->gamma[i]>KBEST_EPS) continue;
	/* find first best hypothesis that is worse than current hypothesis: */
	for (l=0; l<k && maxima[i*k+l] < KBEST_EPS && maxima[i*k+l] > hP->gamma[i]; l++);
	if (l<k) {
	  /* for each m>l: m'th best hypothesis becomes (m+1)'th best */
	  for (m=k-1; m>l; m--) {
	    argmaxs[i*k+m] = argmaxs[i*k+m-1];
	    maxima[i*k+m] = maxima[i*k+m-1];
	  }
	  /* save new l'th best hypothesis: */
	  maxima[i*k+l]  = hP->gamma[i];
	  argmaxs[i*k+l] = hP;
	}
      }
      hP=hP->next;
    }
    
    /* set 'chosen' for all hypotheses from argmaxs array: */
    for (i=0; i < mo->N*k; i++)
      /* only choose hypotheses whose prob. is at least threshold*max_prob */
      if (maxima[i]!=1.0 && maxima[i] >= KBEST_THRESHOLD+maxima[(i%mo->N)*k])
	argmaxs[i]->chosen=1;

    /* remove hypotheses that were not chosen from the lists: */
    /* remove all hypotheses till the first chosen one */
    while (h[t] != NULL && 0==h[t]->chosen)
      hlist_removeElem(&(h[t]));
    /* remove all other not chosen hypotheses */
    hP=h[t];
    while (hP->next != NULL) {
      if (1==hP->next->chosen)
      	hP = hP->next;
      else
	hlist_removeElem(&(hP->next));
    }
  }
  /* dispose of temporary arrays: */
  free(argmaxs);
  free(maxima);
  free(log_a);    /* transition matrix is no longer needed from here on */
  
  /** 4. Save the hypothesis with the highest probability over all states: */
  hP=h[seq_len-1];
  argmax=NULL;
  *log_p = 1.0;   /* log_p will store log of maximum summed probability */
  while (hP != NULL) {
    /* sum probabilities for each hypothesis over all states: */
    sum=logSum(hP->gamma,mo->N);
    /* and select maximum sum */
    if (sum < KBEST_EPS && (*log_p==1.0 || sum>*log_p)) {
      *log_p=sum;
      argmax=hP;
    }
    hP = hP->next;
  }

  /* found a valid path? */
  if (*log_p < KBEST_EPS) {
    /* yes: extract chosen hypothesis: */
    hypothesis=(int*)malloc(sizeof(int)*seq_len);
    for (i=seq_len-1; i>=0; i--) {
      hypothesis[i] = argmax->hyp_c;
      argmax = argmax->parent;
    }
  }
  else
    /* no: return 1.0 representing -INF and an empty hypothesis */
    hypothesis=NULL;

  /* dispose of calculation matrices: */
  hP=h[seq_len-1];
  while (hP!=NULL)
    hlist_removeElem(&hP);
  free(h);
  return hypothesis;
}


