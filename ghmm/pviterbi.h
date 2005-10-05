/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/viterbi.h
*       Authors:  Bernhard Knab, Benjamin Georgi
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

#ifndef PVITERBI_H
#define PVITERBI_H

#ifdef __cplusplus
extern "C" {
#endif

#include "pmodel.h"
#include "psequence.h"
#include "linkedlist.h"


typedef struct plocal_store_t {
  /** precomputed log probabilities for transitions into the states 
      for each transition class **/
  double *** log_in_a;
  /** precomputed log probabilities for each state for the emissions  **/
  double ** log_b;
  /** lookback matrix for the last offset steps **/
  double *** phi;
  /** log probabilities for the current u, v and every state **/
  double *phi_new;
  /** traceback matrix **/
  int ***psi;
  /** for convinience store a pointer to the model **/
  pmodel * mo;
  /** for the debug mode store information of matrix sizes **/
  /** length of sequence X determines size of psi **/
  int len_x;
  /** length of sequence Y determines size of phi and psi **/
  int len_y;
  /** non functional stuff **/
  int    *topo_order;
  int    topo_order_length;
} plocal_store_t;

void ghmm_dp_print_viterbi_store(plocal_store_t * pv);

/**@name Viterbi-Algorithmus */
/*@{ (Doc++-Group: viterbi) */

/**
  Viterbi algorithm. Calculates the Viterbi path (the optimal path trough
  the model) and the Viterbi probability to a given model and a given 
  sequence. The matrices in the local_store struct are allocated using
  stat_matrix_d_alloc.
  @return Viterbi path
  @param mo:    model
  @param o:     sequence
  @param len:   length of the sequence
  @param log_p: probability of the sequence in the Viterbi path
  */
int *ghmm_dp_viterbi(pmodel *mo, psequence * X, psequence * Y, double *log_p, int *path_length);

int *ghmm_dp_viterbi_variable_tb(pmodel *mo, psequence * X, psequence * Y, double *log_p, int *path_length, int start_traceback_with);

int *ghmm_dp_viterbi_test(pmodel *mo, psequence * X, psequence * Y, double *log_p, int *path_length);

/**
  Calculates the logarithmic probability to a given path through the 
  states (does not have to be the Viterbi path), given sequence and
  a model.
  @param mo:        model
  @param o:         sequence
  @param len:       length of the sequence
  @param state_seq: path through the states
  @return log P
  */

double ghmm_dp_viterbi_logp(pmodel *mo, psequence * X, psequence * Y, int *state_seq, int state_seq_len);


#ifdef __cplusplus
}
#endif


#endif

/*@} (Doc++-Group: viterbi) */
