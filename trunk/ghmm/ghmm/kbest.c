/*******************************************************************************
  author       : Anyess von Bock, Alexander Riemer
  filename     : ghmm/ghmm/kbest.c
  created      : TIME: 16:41:02     DATE: Mon 26. April 2004
  $Id$

__copyright__

*******************************************************************************/


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
  int i, j, t, c, l, m;		/* counters */
  int no_oldHyps;		/* number of hypotheses until position t-1 */
  int b_index;			/* index for addressing states' b arrays */
  int no_labels;

  /* logarithmized transition matrix A, log(a(i,j)) => log_a[i*N+j],
       1.0 for zero probability */
  double* log_a;

  /* list of hypotheses */
  hypoList *h, *hP;
  
  /* matrix of dimensions (#states x #hypotheses) that stores probabilities
      of hypotheses for the labeling of subsequences of o_seq,
      gamma(i,c) ==> c'th element of gammaList -> g[i] */
  gammaList *gamma, *gammaP;
  
  /* matrix gamma, optimized with regard to calculation speed,
      gamma(i,c) = oldgamma[i * #(old hypotheses) + c] */
  double* oldgamma;
  
  /* vectors for rows in the matrices */
  int *hypothesis;
  double* gammaRow;
  
  /* indices & prob. of the k most probable hypotheses for each state
       - matrices of dimensions #states x k:  argm(i,l) => argmaxs[i*k+l] */
  int* argmaxs;
  double* maxima;
  
  /* array of booleans specifying which hypotheses were chosen */
  int* chosen;
  
  /* pointer to & probability of most probable hypothesis in a certain state */
  hypoList* argmax;
  double sum;
  
  /* break if sequence empty or k<1: */
  if (seq_len <= 0 || k<=0) return NULL;
  
  /** 1. Initialization (extend empty hypothesis to #labels hypotheses of
         length 1): */

  /* get number of labels (= maximum label + 1): */
  no_labels = -1;  
  for (i=0; i < mo->N; i++) {
    if (mo->s[i].label > no_labels)
      no_labels = mo->s[i].label;
  }
  no_labels++;
  
  /* initialize h & gamma: */
  h = NULL;
  gamma = NULL;
  for (c=no_labels-1; c>=0; c--) {
    /* create #labels hypotheses of max. length seq_len: */
    hypothesis = (int*)malloc(sizeof(int)*seq_len);
    hypothesis[0] = c;    /* initial hypotheses consist of a single label */
    hlist_insertElem(&h,hypothesis);
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
	  gammaRow[i]=1.0;
        else gammaRow[i] = log(mo->s[i].pi)+log(mo->s[i].b[b_index]);
      }
      else gammaRow[i]=1.0;
    }
    glist_insertElem(&gamma,gammaRow);
  }

  /* calculate transition matrix with logarithmic values: */
  log_a = buildLogMatrix(mo->s, mo->N);

  /* initialize temporary arrays: */
  maxima=(double*)malloc(sizeof(double)*mo->N*k); /* for each state save k  */
  argmaxs=(int*)malloc(sizeof(int)*mo->N*k);      /* most prob. hypotheses  */

  /*------ Main loop: Cycle through the sequence: ------*/
  for (t=1; t<seq_len; t++) {
    /* put o_seq[t-1] in emission history: */
    update_emission_history(mo,o_seq[t-1]);
  	
    /** 2. Propagate hypotheses forward and update gamma: */
    no_oldHyps = propFwd(h, no_labels, t, seq_len);
    
    /* build gamma matrix for calculations from the gamma list: */
    oldgamma = (double*)malloc(sizeof(double)*no_oldHyps*mo->N);
    gammaP = gamma;
    for (c=0; c<no_oldHyps; c++) {
      for (i=0; i < mo->N; i++) {
      	oldgamma[i*no_oldHyps+c] = gammaP->g[i];
      }
      gammaP = gammaP->next;
    }

    /*--    update gamma:    --*/
    c=0;
    hP = h;
    gammaP = gamma;
    /* cycle through list of hypotheses and gamma list simultaneously */
    while (hP != NULL) {
      for (i=0; i < mo->N; i++) {
      	if (mo->s[i].label == hP->hyp[t]) {
	  /* if hypothesis c ends with label of state i:
	       gamma(i,c):= log(sum(exp(a(j,i)*exp(oldgamma(j,old_c)))))
	                      + log(b[i](o_seq[t]))
	     else: gamma(i,c):= -INF (represented by 1.0)
	   */
	  gammaP->g[i] = logGammaSum(log_a, i, mo->N,
				     oldgamma, c%no_oldHyps, no_oldHyps);
      	  b_index=get_emission_index (mo,i,o_seq[t],t);
      	  if (b_index<0 || mo->s[i].b[b_index]<KBEST_EPS) gammaP->g[i]=1.0;
	  if (gammaP->g[i]!=1.0) {
	    gammaP->g[i]+= log(mo->s[i].b[b_index]);
	  }
      	}
      	else gammaP->g[i] = 1.0;
      }
      c++;    /* count hypotheses while cycling through the list */
      hP = hP->next;
      /* extend gamma list if necessary: */
      if (gammaP->next == NULL && hP != NULL) {
      	  gammaP->next = (gammaList*)malloc(sizeof(gammaList));
      	  gammaP->next->g = (double*)malloc(sizeof(double)*mo->N);
      	  gammaP->next->next = NULL;
      }
      gammaP = gammaP->next;
    }
    /* dispose of gamma calculation matrix: */
    free(oldgamma);

    /** 3. Choose the k most probable hypotheses for each state and discard all
           hypotheses that were not chosen: */
    
    /* initialize temporary arrays: */
    chosen=(int*)malloc(sizeof(int)*c);    /* #hypotheses is stored in c now */
    for (i=0; i<c; i++) chosen[i]=0;

    for (i=0; i< mo->N*k; i++){
      maxima[i]=1.0;
      argmaxs[i]=0;
    }

    /* cycle through gamma list & calculate the k most probable hypotheses for
       each state: */
    gammaP=gamma;
    c=0;
    while (gammaP != NULL) {
      for (i=0; i < mo->N; i++) {
	if (gammaP->g[i]>KBEST_EPS) continue;
	/* find first best hypothesis that is worse than current hypothesis: */
	for (l=0; l<k && maxima[i*k+l] < KBEST_EPS
	              && maxima[i*k+l] > gammaP->g[i]; l++);
	if (l<k) {
	  /* for each m>l: m'th best hypothesis becomes (m+1)'th best */
	  for (m=k-1; m>l; m--) {
	    argmaxs[i*k+m] = argmaxs[i*k+m-1];
	    maxima[i*k+m] = maxima[i*k+m-1];
	  }
	  /* save new l'th best hypothesis: */
	  maxima[i*k+l]=gammaP->g[i];
	  argmaxs[i*k+l]=c;
	}
      }
      c++;
      gammaP=gammaP->next;
    }

    /* get 'chosen' hypotheses from argmaxs array: */
    for (i=0; i < mo->N*k; i++) {
      /* only choose hypotheses whose prob. is at least threshold*max_prob */
      if (maxima[i]!=1.0 && maxima[i] >= KBEST_THRESHOLD+maxima[(i%mo->N)*k])
	chosen[argmaxs[i]]=1;
    }
    /* remove hypotheses that were not chosen from the lists: */
    for (c=0; !chosen[c]; c++) {
      /* remove hypotheses until the first hypothesis is a chosen one: */
      hlist_removeElem(&h);
      glist_removeElem(&gamma);
    }
    hP=h;
    gammaP=gamma;
    /* remove all other hypotheses that were not chosen: */
    while (hP->next !=NULL) {
      if (!chosen[++c]) {
        hlist_removeElem(&(hP->next));
        glist_removeElem(&(gammaP->next));
      }
      else {
      	hP=hP->next;
      	gammaP=gammaP->next;
      }
    }
    free(chosen);   /* dispose of temporary vector */
  }    
  /* dispose of temporary arrays: */
  free(argmaxs);
  free(maxima);
  free(log_a);    /* transition matrix is no longer needed from here on */
  
  /** 4. Save the hypothesis with the highest probability over all states: */
  gammaP=gamma;
  hP=h;
  argmax=NULL;
  *log_p = 1.0;   /* log_p will store log of maximum summed probability */
  while (gammaP != NULL) {
    /* sum probabilities for each hypothesis over all states: */
    sum=logSum(gammaP->g,mo->N);
    /* and select maximum sum */
    if (sum < KBEST_EPS && (*log_p==1.0 || sum>*log_p)) {
      *log_p=sum;
      argmax=hP;
    }
    gammaP = gammaP->next;
    hP = hP->next;
  }

  /* found a valid path? */
  if (*log_p < KBEST_EPS) {
    /* yes: extract chosen hypothesis: */
    hypothesis=(int*)malloc(sizeof(int)*seq_len);
    for (i=0; i < seq_len; i++)
      hypothesis[i] = argmax->hyp[i];
  }
  else
    /* no: return 1.0 representing -INF and an empty hypothesis */
    hypothesis=NULL;

  /* dispose of calculation matrices: */
  hlist_delete(&h);
  glist_delete(&gamma);

  return hypothesis;
}
