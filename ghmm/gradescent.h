/*******************************************************************************
  author       : Janne Grunau, Alexander Riemer
  filename     : ghmm/ghmm/gradescent.h
  created      : TIME: 12:21:05     DATE: Mon 07. June 2004
  $Id$

__copyright__

*******************************************************************************/


#ifndef GRADESCENT_H
#define GRADESCENT_H

#ifdef __cplusplus
extern "C" {
#endif

/*----------------------------------------------------------------------------*/
double compute_performance(model* mo, sequence_t *sq);

/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences until the training
   converges using gradient descent.
   Model must not have silent states. (checked in Python wrapper)
   @return            0/-1 success/error
   @param mo:         reference to a pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:
 */
int gradient_descent(model** mo, sequence_t *sq);


/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences using gradient descent.
   Model must not have silent states. (checked in Python wrapper)
   @return            0/-1 success/error
   @param mo:         pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:        learning factor
 */
int gradient_descent_onestep(model* mo, sequence_t *sq, double eta);


/*----------------------------------------------------------------------------*/
/**
   computes matrices of n and m variables (expected values for how often a
   certain parameter from A or B is used)
   computes Baum-Welch variables implicit 
   @return                 0/-1 success/error
   @param mo:              pointer to a model
   @param alpha:           matrix of forward variables
   @param backward:        matrix of backward variables
   @param scale:           scaling vector from forward-backward-algorithm
   @param seq:             sequence in internal representation
   @param seq_len:         length of sequence
   @param matrix_b:        matrix for parameters from B (n_b or m_b)
   @param matrix_a:        matrix for parameters from A (n_a or m_a)
   @param vec_pi:          vector for parameters in PI (n_pi or m_pi)
 */
int gradescent_compute_expectations(model* mo, double** alpha, double** beta,
				     double* scale,
				     int* seq, int seq_len, double** matrix_b,
				     double* matrix_a, double* vec_pi);



#ifdef __cplusplus
}
#endif

#endif
