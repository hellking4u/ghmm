/*******************************************************************************
  author       : Alexander Riemer, Janne Grunau
  filename     : ghmm/ghmm/kbestbasics.c
  created      : TIME: 17:56:38     DATE: Tue 22. June 2004
  $Id$

__copyright__

*******************************************************************************/


#include "kbestbasics.h"
#include <math.h>
#include <stdlib.h>

/* inserts new hypothesis into list at position indicated by pointer plist */
inline void hlist_insertElem(hypoList** plist, int newhyp, hypoList* parlist) {

    hypoList* newlist;

    newlist = (hypoList*)calloc(1, sizeof(hypoList));
    newlist->hyp_c = newhyp;
    if (parlist)
      parlist->refcount += 1;
    newlist->parent = parlist;
    newlist->next = *plist;

    *plist = newlist;
}

/* removes hypothesis at position indicated by pointer plist from the list
   removes recursively parent hypothesis with refcount==0 */
inline void hlist_removeElem(hypoList** plist) {

    hypoList* tempPtr = (*plist)->next;

    if ((*plist)->gamma) free((*plist)->gamma);
    if ((*plist)->parent) {
      (*plist)->parent->refcount -= 1;
      if (0==(*plist)->parent->refcount)
	hlist_removeElem(&((*plist)->parent));
    }
    free(*plist);

    *plist = tempPtr;
}

/**
   Propagates list of hypotheses forward by extending each old hypothesis to
   #labels new hypotheses
   @return number of old hypotheses
   @param h:          pointer to list of hypotheses
   @param labels:     number of labels
   @param seq_len:    total sequence length
 */
int propFwd(hypoList* h, hypoList** hplus, int labels, int seq_len) {
  int c;
  int no_oldHyps = 1;
  hypoList* list = h;
  hypoList* end;

  hlist_insertElem(hplus, 0, list);
  end = *hplus;
  list = list->next;

  /* extend original hypotheses */
  for (c=0; c<labels; c++) {
    while (list != NULL) {
      no_oldHyps++;
      hlist_insertElem(&(end->next), c, list);
      list = list->next;
      end = end->next;
    }
    list = h;
  }
  
  return (no_oldHyps/labels);
}


/**
   Calculates the logarithm of sum(exp(log_a[j,a_pos])+exp(log_gamma[j,g_pos]))
   which corresponds to the logarithm of the sum of a[j,a_pos]*gamma[j,g_pos]
   @return logSum for products of a row from gamma and a row from matrix A
   @param log_a:      transition matrix with logarithmic values (1.0 for log(0))
   @param a_pos:      number of row from log_a
   @param N:          width of matrix log_a
   @param gamma:      a row of matrix gamma with logarithmic values (1.0 for log(0))
*/
inline double logGammaSum(double* log_a, int a_pos, int N, double* gamma) {
  double result;
  int j;
  double max=1.0;
  int argmax=0;

  double* logP = (double*)malloc(sizeof(double)*N);

  /* calculate logs of a[k,l]*gamma[k,hi] as sums of logs and find maximum: */
  for (j=0; j<N; j++) {
    if (log_a[j*N+a_pos]<KBEST_EPS && gamma[j]<KBEST_EPS) {
      logP[j] = log_a[j*N+a_pos] + gamma[j];
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
