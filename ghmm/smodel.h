/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/smodel.h
  created      : TIME: 21:44:45     DATE: Sun 14. November 1999
  $Id$

__copyright__

------------------------------------------------------------------------------*/

#ifndef SMODEL_H
#define SMODEL_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/const.h>
#include <ghmm/scanner.h>

/**@name SHMM-Modell */
/*@{ (Doc++-Group: smodel) */

/** Continuous HMM. Structures and function. 
    smodel includes continuous model with one transition matrix 
    (COS  is set to 1) and an extension for
    models with several matrices
    (COS is set to a value greater than 1). In the latter case
    a suitable (depending on the spezific application) function 
    sequence\_get\_class has to be defined */

/**
   typedef density\_t fuer cmodel u. smodel.
*/
typedef enum {
  normal, 
  normal_pos, 
  normal_approx,
  density_number
} density_t;


/** @name sstate
    Structure for one state.
*/
struct sstate {
  /** initial prob. */ 
  double pi;
  /** IDs of successor states */ 
  int *out_id;  
  /** IDs of predecessor states */
  int *in_id;
  /** transition probs to successor states. It is a
   matrix in case of mult. transition matrices (COS > 1)*/
  double **out_a; 
  /** transition probs from predecessor states. It is a
   matrix in case of mult. transition matrices (COS > 1) */ 
  double **in_a;
  /** number of  successor states */     
  int out_states; 
  /** number of  predecessor states */  
  int in_states;    
  /** weight vector for output function components */
  double *c;
  /** mean vector for output functions (normal density and truncated normal
      density */
  double *mue;
  /** variance vector for output functions */
  double *u;
  /** flag for fixation of parameter. If fix = 1 do not change parameters of
      output functions, if fix = 0 do normal training. Default is 0. */
  int fix;
};
typedef struct sstate sstate;

/** @name smodel
    continous HMM    
*/
struct smodel{
  /** Number of states */
  int N;
  /** Number of output densities per state */
  int M;
  /** smodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
  int cos;
  /** Flag for density function. 0: normal density, 1: truncated normal 
      density, 2: approximated normal density */
  density_t density;
  /** prior for a priori prob. of the model. -1 means no prior specified (all
      models have equal prob. a priori. */
  double prior;
  /** All states of the model. Transition probs are part of the states. */
  sstate *s; 
};
typedef struct smodel smodel;

/* don't include this earlier: in sequence.h smodel has to be known */
#include <ghmm/sequence.h>

/** Free memory smodel 
    @return 0: success, -1: error
    @param smo  pointer pointer of smodel */
int     smodel_free(smodel **smo);

/** Reads an ascii file with specifications for one or more smodels.
    All parameters in matrix or vector form.
    This is necessary whenever an initial model is needed (e.g. 
    training) or sequences have to be generated from trained models.
    For each smodel block smodel\_read\_block() is called.
   @return vector of read smodels
   @param filename   input ascii file
   @param smo_number  number of read smodels */
smodel** smodel_read(const char *filename, int *smo_number);

/** Reads one smodel block. It is possible to generate multiple
    identical copies of the model read. Memory allocation is here.
   @return pointer of smode read
   @param s        scanner for reading
   @param multip   number ob identical copies
*/
smodel*  smodel_read_block(scanner_t *s, int *multip);

/**
   Copies one smodel. Memory alloc is here.
   @return pointer to smodel copy
   @param smo   smodel to be copied  */
smodel*  smodel_copy(const smodel *smo);

/**
   Checks if smodel is well definded. E.g. sum pi = 1, only positive values 
   etc.
   @return 0 if smodel is ok, -1 for error
   @param smo   smodel for  checking
*/
int     smodel_check(const smodel* smo);

/**
   For a vector of smodels: check that the number of states and the number
   of output function components are the same in each smodel.
   @return 0 if smodels are  ok, -1 for error
   @param smo    vector of smodels for checking
   @param smodel_number  number of smodels
 */
int     smodel_check_compatibility(smodel **smo, int smodel_number);

/**
   Generates random symbol.
   Generates one random number for a specified state and specified
   output component of the given smodel.
   @return               random number
   @param smo     smodel
   @param state    state
   @param m         index of output component
*/
double smodel_get_random_var(smodel *smo, int state, int m);


/** 
    Produces sequences to a given model. All memory that is needed for the 
    sequences is allocated inside the function. It is possible to define
    the length of the sequences global (global_len > 0) or it can be set 
    inside the function, when a final state in the model is reach (a state
    with no output). If the model has no final state, the sequences will
    have length MAX_SEQ_LEN.
    @return             pointer to an array of sequences
    @param smo:         model
    @param seed:        initial parameter for the random value generator
                        (an integer). If seed == 0, then the random value
			generator is not initialized.
    @param global_len:  length of sequences (=0: automatically via final states)
    @param seq_number:  number of sequences
    @param label:       label tag
    @param Tmax:        maximal sequence length, set to MAX_SEQ_LEN if -1 
*/

sequence_d_t *smodel_generate_sequences(smodel* smo, int seed, int global_len,
					long seq_number, long label, int Tmax);

/** 
    Computes sum over all sequence of
    seq_w * log( P ( O|lambda ) ). If a sequence can't be generated by smo
    error cost of seq_w * PRENALTY_LOGP are imposed.
   @return       n: number of evaluated sequences, -1: error
   @param smo   smodel
   @param sqd    sequence struct
   @param log\_p  evaluated log likelihood
*/
int smodel_likelihood(smodel *smo, sequence_d_t *sqd, double *log_p);

/**
   Prints one smodel in matrix form.
   @param file     output file
   @param smo   smodel
*/
void smodel_print(FILE *file, smodel *smo); 

/**
   Prints one smodel with only one transition Matrix A (=Ak\_0).
   @param file     output file
   @param smo   smodel
*/
void smodel_print_oneA(FILE *file, smodel *smo);

/**
   Prints transition matrix of specified class.
   @param file       output file
   @param smo     smodel
   @param k          transition class
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_Ak_print(FILE *file, smodel *smo, int k, char *tab,
		     char *separator, char *ending);

/**
   Prints weight matrix of output functions of an smodel.
   @param file       output file
   @param smo     smodel
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_C_print(FILE *file, smodel *smo, char *tab, char *separator, 
		    char *ending);

/**
   Prints mean matrix of output functions of an smodel.
   @param file       output file
   @param smo     smodel
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_Mue_print(FILE *file, smodel *smo, char *tab, char *separator, 
		      char *ending);
/**
   Prints variance matrix of output functions of an smodel.
   @param file       output file
   @param smo     smodel
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_U_print(FILE *file, smodel *smo, char *tab, char *separator, 
			char *ending);
/**
   Prints initial prob vector of an smodel.
   @param file       output file
   @param smo     smodel
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_Pi_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending);
/**
   Prints vector of fix\_states.
   @param file       output file
   @param smo     smodel
   @param tab      format: leading tab
   @param separator  format: seperator
   @param ending     format: end of data in line
*/
void smodel_fix_print(FILE *file, smodel *smo, char *tab, char *separator, 
		     char *ending);

/** Computes the density of one symbol (omega) in a given state and a 
    given output component
    @return calculated density
    @param smo smodel
    @param state state 
    @param m output component
    @param omega given symbol
*/
double smodel_calc_cmbm(smodel *smo, int state, int m, double omega);

/** Computes the density of one symbol (omega) in a given state (sums over
    all output components
    @return calculated density
    @param smo smodel
    @param state state 
    @param omega given symbol
*/
double smodel_calc_b(smodel *smo, int state, double omega);

/** Computes probabilistic distance of two models
    @return the distance
    @param cm0  smodel used for generating random output
    @param cm   smodel to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal langth of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
double smodel_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, 
			    int verbose);

/** 
    Computes value of distribution function for a given symbol omega, a given
    state and a given output component.
    @return   value of distribution function
    @param smo   smodel
    @param state  state
    @param m      component
    @param omega symbol
*/
double smodel_calc_cmBm(smodel *smo, int state, int m, double omega);

/** 
    Computes value of distribution function for a given symbol omega and
    a given  state. Sums over all components.
    @return   value of distribution function
    @param smo   smodel
    @param state  state
    @param omega symbol
*/
double smodel_calc_B(smodel *smo, int state, double omega);

/** Computes the number of free parameters in an array of
   smodels. E.g. if the number of parameter from pi is N - 1.
   Counts only those parameters, that can be changed during  
   training. If pi[i] = 0 it is not counted, since it can't be changed.
   @return number of free parameters
   @param smo smodel
   @param smo\_number number of smodels
*/
int smodel_count_free_parameter(smodel **smo, int smo_number);


/*============================================================================*/

/* keep the following functions for first distribution???
   --> BK ? 
*/


/** Generates interval(a,b) with  B(a) < 0.01, B(b) > 0.99
    @param smo    continous HMM
    @param state  given state
    @param a      return-value: left side
    @param b      return-value: right side
*/
void smodel_get_interval_B(smodel *smo, int state, double *a, double *b);

#ifdef __cplusplus
}
#endif


#endif

/*@} (Doc++-Group: smodel) */


