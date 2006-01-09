/*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmm/ghmm/sdmodel.h
*       Authors:  Wasinee Rungsarityotin, Benjamin Georgi
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
#ifndef GHMM_SDMODEL_H
#define GHMM_SDMODEL_H

#include "ghmm.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@name HMM-Modell */
/*@{ (Doc++-Group: model) */

/** @name ghmm_dsstate
    The basic structure, keeps all parameters that belong to a state. 
*/
  struct ghmm_dsstate {
  /** Initial probability */
    double pi;
  /** Output probability */
    double *b;
  /** ID of the following state */
    int *out_id;
  /** ID of the previous state */
    int *in_id;

  /** transition probs to successor states. It is a
   matrix in case of mult. transition matrices (COS > 1)*/
    double **out_a;
  /** transition probs from predecessor states. It is a
   matrix in case of mult. transition matrices (COS > 1) */
    double **in_a;

  /** Number of successor states */
    int out_states;
  /** Number of precursor states */
    int in_states;
  /** if fix == 1 --> b stays fix during the training */
    int fix;
    /* XXX Specific variable for ProfileHMM to count the number of
       match states. Not elegant solution.
       WS: if 1 then counts me, 0 don't count me */
    int countme;

  /** Position for graphical editing */
    int xPosition;
    int yPosition;
  };
  typedef struct ghmm_dsstate ghmm_dsstate;

/** @name ghmm_dmodel
    The complete HMM. Contains all parameters, that define a HMM.
*/
  struct ghmm_dsmodel {
  /** Number of states */
    int N;
  /** Number of outputs */
    int M;
  /** ghmm_dsmodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
    int cos;
  /** Vector of the states */
    ghmm_dsstate *s;
  /** Prior for the a priori probability for the model.
      A value of -1 indicates that no prior is defined. */
    double prior;

  /** Contains bit flags for various model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */

  /** pointer to class function   */
    int (*get_class) (int *, int);

  /** Contains bit flags for various model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
    int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
    int *silent;

    int topo_order_length; /*AS*/
    int *topo_order; /*WR*/

  /** Store for each state a class label. Limits the possibly state sequence
      
      Note: label != NULL iff (model_type & kLabeledStates) != 0  */
    int * label;
    alphabet_s * labelAlphabet;

    alphabet_s * alphabet;
  };
  typedef struct ghmm_dsmodel ghmm_dsmodel;


#ifdef __cplusplus
}
#endif
/*
  Important: The inclusion of sequence.h ist not done before this point in
  order to avoid error by compiling.
*/
#include "sequence.h"
#include "scanner.h"
#ifdef __cplusplus
extern "C" {
#endif


  /** Frees the memory of a model.
      @return 0 for succes; -1 for error
      @param mo:  pointer to a ghmm_dsmodel */
  int ghmm_ds_free (ghmm_dsmodel ** mo);

  int ghmm_ds_init_silent_states (ghmm_dsmodel * mo);

  /** 
      Produces sequences to a given model. All memory that is needed for the 
      sequences is allocated inside the function. It is possible to define
      the length of the sequences global (global_len > 0) or it can be set 
      inside the function, when a final state in the model is reach (a state
      with no output). If the model has no final state, the sequences will
      have length MAX_SEQ_LEN.
      @return             pointer to an array of sequences
      @param mo:          model
      @param seed:        initial parameter for the random value generator
      (an integer). If seed == 0, then the random value
      generator is not initialized.
      @param global_len:  length of sequences (=0: automatically via final states)
      @param seq_number:  number of sequences
      @param T_max:  maximal number of consecutive silent states in model (used to
	  identify silent circles).
  */
  ghmm_dseq * ghmm_ds_generate_sequences(ghmm_dsmodel * mo, int seed,
					 int global_len, long seq_number,
					 int Tmax);

  /**
     Copies a given model. Allocates the necessary memory.
     @returns a copy of the model
     @param mo:  dmodel to copy */
  ghmm_dsmodel * ghmm_ds_copy(const ghmm_dsmodel * mo);

  /** Utility for converting between single discrete model and switching model */
  ghmm_dmodel * ghmm_ds_to_dmodel(const ghmm_dsmodel * mo, int kclass);

  /** */
  void ghmm_ds_from_dmodel (const ghmm_dmodel * mo, ghmm_dsmodel * smo, int klass);

  /**
     Writes a model in matrix format.
     @param file: output file
     @param mo:   model
  */
  void sdmodel_print (FILE * file, ghmm_dsmodel * mo);


  /**
     Writes transition matrix of a model.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void ghmm_ds_Ak_print (FILE * file, ghmm_dsmodel * mo, int k, char *tab,
                         char *separator, char *ending);
  /**
     Writes output matrix of a model.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void ghmm_ds_B_print (FILE * file, ghmm_dsmodel * mo, char *tab, char *separator,
                        char *ending);

  /**
     Writes initial allocation vector of a matrix.
     @param file: output file
     @param mo:   model
     @param tab:  format: leading Tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void ghmm_ds_Pi_print (FILE * file, ghmm_dsmodel * mo, char *tab,
                         char *separator, char *ending);

  /*============================================================================*/
  /** ghmm_ds_viterbi is working for switching discrete model
   *  ghmm_ds_topo_order -- need to be implemented with DFS (as in model_util.c)
   *============================================================================
   **/

  void ghmm_ds_topo_order (ghmm_dsmodel * mo);

  int *ghmm_ds_viterbi (ghmm_dsmodel * mo, int *o, int len, double *log_p);

  /** Forward-Algorithm.
      Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
      a given double sequence and a given model.
      @param mo:      model
      @param O        sequence
      @param length: length of sequence
      @param alpha:  alpha[t][i]
      @param scale:   a reference for double type, scale factors
      @param log\_p:  a reference for double type, log likelihood log( P(O|lambda) )
      @return 0 for success, -1 for error
  */
  int ghmm_ds_forward (ghmm_dsmodel * mo, const int *O, int len, double **alpha,
                      double *scale, double *log_p);


  /** Descale
      descales the alpha matrix from the forward algorithm
      @param alpha: alpha matrix from forward
      @param scale: scale vector from forward
      @param t:     number of timesteps
      @param n:     number of states
      @param newalpha: unscaled alpha matrix
      @return 0 for success, -1 for error
  */
  int ghmm_ds_forward_descale (double **alpha, double *scale, int t, int n,
                      double **newalpha);


/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
  double ghmm_ds_likelihood (ghmm_dsmodel * mo, ghmm_dseq * sq);


/** 
    Writes the parameters of a model sorted by states. 
    Is not very concise.   
    @param file: output file
    @param mo:   model
*/
  void ghmm_ds_states_print (FILE * file, ghmm_dsmodel * mo);


#ifdef __cplusplus
}
#endif
#endif
/*@} (Doc++-Group: model) */
