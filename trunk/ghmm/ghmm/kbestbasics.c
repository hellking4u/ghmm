/*******************************************************************************
  author       : Alexander Riemer
  filename     : ghmm/ghmm/kbestbasics.c
  created      : TIME: 17:56:38     DATE: Tue 22. June 2004
  $Id$

__copyright__

*******************************************************************************/


#include "kbestbasics.h"
#include <math.h>
#include <stdlib.h>

/* inserts new hypothesis into list at position indicated by pointer plist */
inline void hlist_insertElem(hypoList** plist, int* newhyp) {
    hypoList* newlist;
    newlist = (hypoList*)malloc(sizeof(hypoList));
    newlist->hyp = newhyp;
    newlist->next = *plist;
    *plist = newlist;
}

/* removes hypothesis at position indicated by pointer plist from the list */
inline void hlist_removeElem(hypoList** plist) {
    hypoList* tempPtr=(*plist)->next;
    free((*plist)->hyp);
    free(*plist);
    *plist=tempPtr;
}

/* deletes the whole list of hypotheses */
inline void hlist_delete(hypoList** list){
    hypoList* temp_ptr;
    while (*list != NULL){
        free((*list)->hyp);
        temp_ptr = (*list)->next;
        free(*list);
        *list = temp_ptr;
    }
}

/* inserts new row into gamma list at position indicated by pointer plist */
inline void glist_insertElem(gammaList** plist, double* newG) {
    gammaList* newlist;
    newlist = (gammaList*)malloc(sizeof(gammaList));
    newlist->g = newG;
    newlist->next = *plist;
    *plist = newlist;
}

/* removes row at position indicated by pointer plist from gamma list */
inline void glist_removeElem(gammaList** plist) {
    gammaList* tempPtr=(*plist)->next;
    free((*plist)->g);
    free(*plist);
    *plist=tempPtr;
}

/* deletes the whole gamma list */
inline void glist_delete(gammaList** list){
    gammaList* temp_ptr;
    while (*list != NULL){
    	free((*list)->g);
        temp_ptr = (*list)->next;
        free(*list);
        *list = temp_ptr;
    }
}


/**
   Propagates list of hypotheses forward by extending each old hypothesis to
   #labels new hypotheses
   @return number of old hypotheses
   @param h:          pointer to list of hypotheses
   @param labels:     number of labels
   @param t:          last position in current sub-sequence
   @param seq_len:    total sequence length
 */
inline int propFwd(hypoList* h, int labels, int t, int seq_len) {
  int c,d,i;
  int no_oldHyps=0;
  int* hypothesis;
  hypoList* list = h;
  hypoList* end;
  
  while (list != NULL) {
    /* extend original hypotheses by first label (=0) */
    list->hyp[t]=0;
    no_oldHyps++;
    end = list;
    list = list->next;
  }
  for (c=1; c<labels; c++) {
    /* start at beginning of original list */
    list = h;
    for (d=0; d<no_oldHyps; d++) {
      /* create new hypothesis by extending an old one: */
      hypothesis = (int*)malloc(sizeof(int)*seq_len);
      for (i=0; i<t; i++) {
        hypothesis[i] = list->hyp[i];
      }
      hypothesis[t] = c;
      hlist_insertElem(&(end->next),hypothesis);
      end = end->next;
      list = list->next;
    }
  }
  return no_oldHyps;
}


/**
   Calculates the logarithm of sum(exp(log_a[j,a_pos])+exp(log_gamma[j,g_pos]))
   which corresponds to the logarithm of the sum of a[j,a_pos]*gamma[j,g_pos]
   @return logSum for products of a row from gamma and a row from matrix A
   @param log_a:      transition matrix with logarithmic values (1.0 for log(0))
   @param a_pos:      number of row from log_a
   @param N:          width of matrix log_a
   @param gamma:      matrix gamma with logarithmic values (1.0 for log(0))
   @param g_pos:      number of row in gamma
   @param no_oldHyps: width of matrix gamma
*/
inline double logGammaSum(double* log_a, int a_pos, int N, double* gamma,
			  int g_pos, int no_oldHyps) {
  double result;
  int j;
  double max=1.0;
  int argmax=0;

  double* logP = (double*)malloc(sizeof(double)*N);

  /* calculate logs of a[k,l]*gamma[k,hi] as sums of logs and find maximum: */
  for (j=0; j<N; j++) {
    if (log_a[j*N+a_pos]<KBEST_EPS && gamma[j*no_oldHyps+g_pos]<KBEST_EPS) {
      logP[j] = log_a[j*N+a_pos] + gamma[j*no_oldHyps+g_pos];
      if (logP[j]<KBEST_EPS && (max==1.0 || logP[j]>max)) {
	max = logP[j];
	argmax=j;
      }
    }
    else logP[j] = 1.0;
  }

  /* calculate max+log(1+sum[j!=argmax; exp(logP[j]-max)])  */
  result = 1.0;
  for (j=0; j<N; j++) {
    if (logP[j]<KBEST_EPS && j!=argmax) result+= exp(logP[j]-max);
  }
  result=log(result);
  result+=max;

  free(logP);
  return result;
}


/**
   Calculates the logarithm of the sum of the probabilities whose logarithms are
   stored in the given array
   @return log of sum of exp(a[i])
   @param a:          array of logarithms of probabilities (a[i] < 0 for all i)
   @param N:          length of a
*/
inline double logSum(double* a, int N) {
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
}
