/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/kbestbasics.h
*       Authors:  Alexander Riemer, Janne Grunau
*
*       Copyright (C) 1998-2004 Alexander Schliep 
*       Copyright (C) 1998-2001 ZAIK/ZPR, Universitaet zu Koeln
*	Copyright (C) 2002-2004 Max-Planck-Institut fuer Molekulare Genetik, 
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
#ifndef KBESTBASICS_H
#define KBESTBASICS_H


#ifdef __cplusplus
extern "C" {
#endif


/* threshold probability (logarithmized) */
#define KBEST_THRESHOLD -3.50655789732
           /* log(0.03) => threshold: 3% of most probable partial hypothesis */
#define KBEST_EPS 1E-15

/** Data type for linked list of hypotheses
    Stores the actual label, a link to the parent hypothesis, a counter of the
    links to this hypothesis, a the gamm values */
typedef struct hypo_List {
  int hyp_c;
  int refcount;
  int chosen;
  double* gamma;
  struct hypo_List* next;
  struct hypo_List* parent;
} hypoList;

/* inserts new hypothesis into list at position indicated by pointer plist */
inline void hlist_insertElem(hypoList** plist, int newhyp, hypoList* parlist);

/* removes hypothesis at position indicated by pointer plist from the list */
inline void hlist_removeElem(hypoList** plist);

/**
   Propagates list of hypotheses forward by extending each old hypothesis to
   #labels new hypotheses
   @return number of old hypotheses
   @param h:          pointer to list of hypotheses
   @param labels:     number of labels
   @param seq_len:    total sequence length
 */
int propFwd(hypoList* h, hypoList** hplus, int labels, int seq_len);


/**
   Calculates the logarithm of sum(exp(log_a[j,a_pos])+exp(log_gamma[j]))
   which corresponds to the logarithm of the sum of a[j,a_pos]*gamma[j]
   @return logSum for products of a row from gamma and a row from matrix A
   @param log_a:      transition matrix with logarithmic values (1.0 for log(0))
   @param a_pos:      number of row from log_a
   @param N:          width of matrix log_a
   @param gamma:      a row of matrix gamma with logarithmic values (1.0 for log(0))
*/
inline double logGammaSum(double* log_a, int a_pos, int N, double* gamma);


/**
   Calculates the logarithm of the sum of the probabilities whose logarithms are
   stored in the given array
   @return log of sum of exp(a[i])
   @param a:          array of logarithms of probabilities (a[i] < 0 for all i)
   @param N:          length of a
*/
inline double logSum(double* a, int N);

#ifdef __cplusplus
}
#endif

#endif
