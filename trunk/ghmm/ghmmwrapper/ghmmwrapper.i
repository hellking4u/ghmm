/* /*******************************************************************************
*
*       This file is part of the General Hidden Markov Model Library,
*       GHMM version __VERSION__, see http://ghmm.org
*
*       Filename: ghmmwrapper.i
*       Authors:   Wasinee Rungsarityotin, Benjamin Georgi, Janne Grunau
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


%module ghmmwrapper

/* this will not work if we build ghmmwrapper out of the ghmm tree
   include ../config.h twice, once for swig and once for CC */
%include "../config.h"
%{
#include "../config.h"
#include <stdio.h>
#include <ghmm/ghmm.h>
#include <ghmm/sequence.h>
#include <ghmm/smodel.h>
#include <ghmm/sreestimate.h>
#include <ghmm/rng.h>
#include <ghmm/reestimate.h>
#include <ghmm/foba.h>
#include <ghmm/scluster.h>
#include <ghmm/gradescent.h>
#include <ghmm/kbest.h>
#include "sclass_change.h"
#include <ghmm/psequence.h>
#include <ghmm/pmodel.h>
#include <ghmm/pviterbi.h>
#include <ghmm/pviterbi_propagate.h>
#include "pclasschange.h"
#include <ghmm/obsolete.h>
%}

%include carrays.i
%include cmalloc.i
%include cpointer.i
%include cstring.i
// %pointer_functions(int, intp)

%include constraints.i
%include exception.i    
%include typemaps.i

// Constraints on GHMM date types - no NULL pointers as function arguments
%apply Pointer NONNULL { ghmm_dmodel * };
%apply Pointer NONNULL { ghmm_dmodel ** };
%apply Pointer NONNULL { ghmm_cmodel * };
%apply Pointer NONNULL { ghmm_cmodel ** };
%apply Pointer NONNULL { ghmm_dstate * };
%apply Pointer NONNULL { ghmm_cstate * };
%apply Pointer NONNULL { ghmm_dseq * };
%apply Pointer NONNULL { ghmm_dseq ** };
%apply Pointer NONNULL { ghmm_cseq * };
%apply Pointer NONNULL { ghmm_cseq ** };
%apply Pointer NONNULL { scluster * };
 
// Constraints on general C data types - no NULL pointers as function arguments
// XXX certain arguments are supposed to be NULL
// %apply Pointer NONNULL { int *, double *, int **, double **, void * };

/*==========================================================================
  ========================== test for obsolete features ==================== */
#ifdef GHMM_OBSOLETE

#define SMO_FILE_SUPPORT 1
#define ASCI_SEQ_FILE    1
#define SEQ_LABEL_FIELD  1

#else

#define SMO_FILE_SUPPORT 0
#define ASCI_SEQ_FILE    0
#define SEQ_LABEL_FIELD  0

#endif /* GHMM_OBSOLETE */


/*=============================================================================================
  =============================== Random Number Generator (RNG) ================================= */

/* The global RNG */
//extern gsl_rng * RNG; 
extern GHMM_RNG * RNG;

/* Important! initialise rng  */
extern void ghmm_rng_init(void);

/* Initialise random timeseed */
extern void ghmm_rng_timeseed(GHMM_RNG * r);

%inline %{
	void time_seed(){
		ghmm_rng_timeseed(RNG);
	}
%}
/*=============================================================================================
  =============================== Logging Function Callbacks ================================= */

extern void ghmm_set_logfunc(void (* fptr)(int, const char *, void *), void * clientdata);

// Grab a Python function object as a Python object.
%typemap(python,in) PyObject *pyfunc {
  if (!PyCallable_Check($input)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }
  $1 = $input;
}

%{
  /* This function matches the prototype of the normal C callback
     function for our widget. However, we use the clientdata pointer
     for holding a reference to a Python callable object. */
  
  static void PythonCallBack(int level, const char * message, void *clientdata)
    {
      PyObject *func, *arglist;
      
      // Get Python function
      func = (PyObject *) clientdata;
      // Build Python arguments
      arglist = Py_BuildValue("(is)", level, message);
      // Call Python
      PyEval_CallObject(func, arglist);
      // Trash arguments
      Py_DECREF(arglist);
    }
%}

%inline %{
// Set a Python function object as a callback function
// Note : PyObject *pyfunc is remapped with a typempap
static void set_pylogging(PyObject *pyfunc) {
  ghmm_set_logfunc(PythonCallBack, (void *) pyfunc);
  Py_INCREF(pyfunc);
}
%}


				
/*=============================================================================================
  =============================== sequence.c  ============================================== */
		
		
/**@name sequences  (double and int) */
/*@{ (Doc++-Group: sequence) */
/** @name struct ghmm_dseq
    Sequence structure for integer sequences. 
    Contains an array of sequences and corresponding
    data like sequence label, sequence weight, etc. Sequences may have different
    length.    
 */

struct ghmm_dseq {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq.
   */
  int **seq;
  /** matrix of state ids, can be used to save the viterbi path during sequence generation.
   ATTENTION: is NOT allocated by ghmm_dseq_calloc  */
  int **states;
  /** array of sequence length */
  int *seq_len;
#ifdef GHMM_OBSOLETE
  /**  array of sequence labels */
  long *seq_label;
#endif /* GHMM_OBSOLETE */
  /**  array of sequence IDs*/
  double *seq_id;
  /** positiv! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;
  
  /* matrix of state labels corresponding to seq */
  int **state_labels; 
  /* number of labels for each sequence */  
  int *state_labels_len;        
};
typedef struct ghmm_dseq ghmm_dseq;

%pointer_functions(ghmm_dseq **, ghmm_dseq_setPtr)

/** @name struct ghmm_cseq
    Sequence structure for double sequences. 
    Contains an array of sequences and corresponding
    data like sequnce label, sequence weight, etc. Sequences may have different
    length.    
 */
struct ghmm_cseq {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq. */
  double **seq;
  /** array of sequence length */
  int *seq_len;
#ifdef GHMM_OBSOLETE
  /**  array of sequence labels */
  long *seq_label;
#endif /* GHMM_OBSOLETE */
  /**  array of sequence IDs*/
  double *seq_id;
  /** positive! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;  
};
typedef struct ghmm_cseq ghmm_cseq;

%pointer_functions(ghmm_cseq **, ghmm_cseq_setPtr)

/**
   Memory allocation for an integer sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
extern ghmm_dseq *ghmm_dseq_calloc(long seq_number);


/**
   Memory allocation for a double  sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
extern ghmm_cseq *ghmm_cseq_calloc(long seq_number);

/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
extern void ghmm_dseq_clean(ghmm_dseq *sq);

/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
extern void ghmm_cseq_clean(ghmm_cseq *sq);

/**
  Prints one array of integer sequences in a file.
  @param file       output file
  @param sequence    array of sequences
  */
void ghmm_dseq_print(FILE *file, ghmm_dseq *sequence);


/**
  copy one integer sequence. Memory for target has to be allocated outside.
  @param target  target sequence
  @param source source sequence
  @param len     length of source sequence
  */
void ghmm_dseq_copy(int *target, int *source, int len);


/**
  copy one double sequence. Memory for target has to be allocated outside.
  @param target  target sequence
  @param source source sequence
  @param len     length of source sequence
  */
void ghmm_cseq_copy(double *target, double *source, int len);


#ifdef GHMM_OBSOLETE
/**
   Reads one or several arrays of double sequences. 
   Calls sequence\_read\_alloc, where reading
   and memory allocation is done. 
   @return pointer to sequence array
   @param filename    input filename
*/
ghmm_cseq **ghmm_cseq_read(const char *filename, int *sqd_number);
#endif /* GHMM_OBSOLETE */

/**
  Adds all integer sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
int ghmm_dseq_add(ghmm_dseq *target, ghmm_dseq *source);

/**
  Adds all double sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
extern int ghmm_cseq_add(ghmm_cseq *target, ghmm_cseq *source);

/**
  Frees all memory in a given array of integer sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
int ghmm_dseq_free(ghmm_dseq **sq);

/**
  Frees all memory in a given array of double sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
int ghmm_cseq_free(ghmm_cseq **sq);

extern void ghmm_cseq_print(FILE *file, ghmm_cseq *sqd, int discrete);


/**
  Extract a single sequence from a larger ghmm_dseq into a new struct.
  
  @return ghmm_dseq struct containing a single sequence
  @param sq   source ghmm_dseq
  @param index   index of sequence to extract
*/
extern ghmm_dseq *ghmm_dseq_get_singlesequence(ghmm_dseq *sq, int index);

/**
  Extract a single ghmm_cseq from a larger ghmm_cseq into a new struct.
  
  @return ghmm_cseq struct containing a single sequence
  @param sq   source ghmm_cseq
  @param index   index of sequence to extract
*/
extern ghmm_cseq *ghmm_cseq_get_singlesequence(ghmm_cseq *sq, int index);

/**
  Free a ghmm_dseq struct which holds as sequence a reference to a sequence in a different
  sequence_t. The function deallocates everything but the reference.
*/
extern int ghmm_dseq_subseq_free (ghmm_dseq ** sq);

/**
  Free a ghmm_cseq struct which holds as sequence a reference to a sequence in a different
  sequence_d_t. The function deallocates everything but the reference.
*/
extern int ghmm_cseq_subseq_free (ghmm_cseq ** sqd);





/*** !!!!!!!! TO DO: Free functions for all types of pointers !!!!!!!! Why? free() does not depend on the data type ***/

%inline%{

  /* return a C-pointer to an integer sequence */
  int *get_onesequence(ghmm_dseq *seqpt, int seqnumber) { 
    return (int *) seqpt->seq[seqnumber];
  }


#ifdef GHMM_OBSOLETE
  ghmm_dseq *seq_read(char* filename ){
	  int i;
	  ghmm_dseq** s;
	  //s = (ghmm_cseq **) malloc(1*sizeof(ghmm_cseq*));
	  s = ghmm_dseq_read(filename, &i);
	  return s[0];
  }
#endif /* GHMM_OBSOLETE */
  
  ghmm_dseq* get_seq_ptr(ghmm_dseq** array, int index){
	  return (ghmm_dseq*) array[index];
  }
  
  /*** create and manipulate an array of pointers of pointers to ghmm_cseq structs ***/
  ghmm_cseq **sequence_d_t_array(int size) {
     int i;
	 ghmm_cseq **s;
	 s = (ghmm_cseq **) malloc(size*sizeof(ghmm_cseq*));
	 for(i =0;i<size;i++){	
		 s[i] = NULL;
	 }
	 return s;	 
  }

  ghmm_cseq* get_seq_d_ptr(ghmm_cseq** array, int index){
	  return (ghmm_cseq*) array[index];
  }	  
  
  void set_seq_d_array(ghmm_cseq** array, int index,ghmm_cseq* seq){
	array[index] = seq;
  }	

#ifdef GHMM_OBSOLETE
  void set_sequence_d_label(ghmm_cseq* seq, int seq_num, long label){
	  seq->seq_label[seq_num] = label;
  }
  
  long get_sequence_d_label(ghmm_cseq* seq, int seq_num ){
	  return seq->seq_label[seq_num];
  }	  
   
  ghmm_cseq *seq_d_read(char* filename ){
	  int i;
      ghmm_cseq** s;
	  //s = (ghmm_cseq **) malloc(1*sizeof(ghmm_cseq*));
	  s = ghmm_cseq_read(filename, &i);
	  return *s;
  }
#endif /* GHMM_OBSOLETE */

  void call_ghmm_dseq_print (char* ch, ghmm_dseq* seq){
    FILE* file_name;
    file_name = fopen (ch, "at");
    ghmm_dseq_print (file_name, seq);
    fclose (file_name);
  }	  

  void call_ghmm_cseq_print (char* ch,ghmm_cseq* seq, int disc){
    FILE* file_name;
    file_name = fopen (ch, "at");
    ghmm_cseq_print (file_name, seq, disc);
    fclose (file_name);
  }	  


  void call_ghmm_dseq_free(ghmm_dseq *sq) {ghmm_dseq_free(&sq);}
  void call_ghmm_cseq_free(ghmm_cseq *sq) {ghmm_cseq_free(&sq);}

  void call_ghmm_dseq_subseq_free(ghmm_dseq *sq) {ghmm_dseq_subseq_free (&sq);}
  void call_ghmm_cseq_subseq_free(ghmm_cseq *sq) {ghmm_cseq_subseq_free (&sq);}

   
%}

/*=============================================================================================
  =============================== model.c  ============================================== */

/** @name ghmm_d_background_distributions
    A container for background distributions to be used in the reestimation. Model
    has an ID (== index) to be used for the arrays background_distributions.order
    and background_distributions.b
*/
struct ghmm_d_background_distributions {
  /** Number of distributions */
  int n;
  /** Number of symbols in alphabet */
  int m;
  /** Order of the respective distribution */
  int* order;
  /** The probabilities */ 
  double **b;
};
typedef struct ghmm_d_background_distributions ghmm_d_background_distributions;


/** @name ghmm_dstate
    The basic structure, keeps all parameters that belong to a state. 
*/
struct ghmm_dstate {
  /** Initial probability */ 
  double pi;
  /** Output probability */
  double *b;
  int order;
  
  /** IDs of the following states */ 
  int *out_id;  
  /** IDs of the previous states */    
  int *in_id;

  /** transition probs to successor states. */
  double *out_a; 
  /** transition probs from predecessor states. */ 
  double *in_a;

  /** Transition probability to a successor 
      double *out_a; */
  /** Transition probablity to a precursor 
      double *in_a;*/

  /** Number of successor states */     
  int out_states; 
  /** Number of precursor states */
  int in_states;  
  /** if fix == 1 --> b stays fix during the training */
  int fix;

  int label;  
};
typedef struct ghmm_dstate ghmm_dstate;


struct ghmm_coord_t {
  double x;
  double y;
};
typedef struct ghmm_coord_t ghmm_coord_t;


/** @name ghmm_dmodel
    The complete HMM. Contains all parameters, that define a HMM.
*/
  struct ghmm_dmodel {
  /** Number of states */
    int N;
  /** Number of outputs */
    int M;
  /** Vector of the states */
    ghmm_dstate *s;
  /** The a priori probability for the model.
      A value of -1 indicates that no prior is defined. 
      Note: this is not to be confused with priors on emission
      distributions*/
    double prior;

    /* contains a arbitrary name for the model */
    char *name;

   /** Contains bit flags for varios model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
    int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
    int *silent;
      /*AS*/
  /** Int variable for the maximum level of higher order emissions */
    int maxorder;
  /** saves the history of emissions as int, 
      the nth-last emission is (emission_history * |alphabet|^n+1) % |alphabet|
      see ...*/
    int emission_history;

  /** Flag variables for each state indicating whether the states emissions
      are tied to another state. Groups of tied states are represented
      by their tie group leader (the lowest numbered member of the group).
      
      tied_to[s] == kUntied  : s is not a tied state
      
      tied_to[s] == s        : s is a tie group leader

      tied_to[t] == s        : t is tied to state s

      Note: tied_to != NULL iff (model_type & kTiedEmissions) == 1  */
    int *tied_to;

  /** Note: State store order information of the emissions.
      Classic HMMS have emission order 0, that is the emission probability
      is conditioned only on the state emitting the symbol.

      For higher order emissions, the emission are conditioned on the state s
      as well as the previous emission_order[s] observed symbols.

      The emissions are stored in the state's usual double* b. The order is
      set state.order.

      Note: state.order != NULL iff (model_type & kHigherOrderEmissions) == 1  */

  /** bp is a pointer to a ghmm_d_background_distributions structure,
      which holds (essentially) an array of background distributions
      (which are just vectors of floating point numbers like state.b).

      For each state the array background_id indicates which of the background
      distributions to use in parameter estimation. A value of kNoBackgroundDistribution
      indicates that none should be used.


      Note: background_id != NULL iff (model_type & kHasBackgroundDistributions) == 1  */
    int *background_id;
    ghmm_d_background_distributions *bp;

  /** (WR) added these variables for topological ordering of silent states 
      Condition: topo_order != NULL iff (model_type & kSilentStates) == 1
   */
    int *topo_order;
    int topo_order_length;

  /** pow_lookup is a array of precomputed powers

      It contains in the i-th entry M (alphabet size) to the power of i
      The last entry is maxorder+1
  */
    int *pow_lookup;
    
    
    
    /*  storage of model representation information ( read in from XML ) */
    
    /* emission alphabet  */ 
    char*** alphabet;
    
    /* sizes of the different alphabets (for pair HMMs there might be more than one)  */
    int*  alphabet_size;
    
    /* number of entries in alphabet_size */
    int S;
    
    /* an arry of positions of states for graphical representation */ 
    ghmm_coord_t *position;
    
    /* state label alphabet (only for labelled HMMs)  */ 
    char** label_alphabet;
 
    /* size of the label_alphabet */
    int label_size;

  };
  typedef struct ghmm_dmodel ghmm_dmodel;
  

/** Frees the memory of a model.
    @return 0 for succes; -1 for error
    @param mo:  address of a pointer to a ghmm_dmodel 
*/
extern int     ghmm_d_free(ghmm_dmodel **mo);

#ifdef GHMM_OBSOLETE
/**
   Reads in ASCII data to initialize an array of models. Memory allocation for
   the models is done here.
   @return array of pointers to the models
   @param filename:   the ASCII input file
   @param mo_number:  filled with number of models read */
extern ghmm_dmodel** ghmm_d_read(char *filename, int *mo_number);
#endif /* GHMM_OBSOLETE */

/**
   Writes a model in matrix format.
   @param file: output file
   @param mo:   model
*/
extern void ghmm_d_print(FILE *file, ghmm_dmodel *mo); 


#ifdef GHMM_OBSOLETE
/**
   Produces simple left-right models given sequences. 
   The function "ghmm_d_generate_from_sequence" is called for each 
   model that should be made. The sequences are read in from the
   ASCII file and thrown away again when leaving the function.
   @return vector of models
   @param s:          scanner
   @param new_models: number of models to produce */
extern ghmm_dmodel **ghmm_d_from_sequence_ascii(scanner_t *s, long *mo_number);
#endif  /* GHMM_OBSOLETE */

/** 
    Produces simple left-right models given sequences. The sequences
    are not read in from file, but exists already as a structur.
    @return vector of models
    @param s:          scanner
    @param new_models: number of models to produce */
extern ghmm_dmodel **ghmm_d_from_sequence(ghmm_dseq *sq, long *mo_number);

/**
   Copies a given model. Allocates the necessary memory.
   @return copy of the model
   @param mo:  model to copy */
extern ghmm_dmodel*  ghmm_d_copy(const ghmm_dmodel *mo);

/**
   Tests if all standardization requirements of model are fulfilled. 
   (That is, if the sum of the probabilities is 1).
   @return 0 for succes; -1 for error
   @param mo:  model to test */
extern int     ghmm_d_check(const ghmm_dmodel* mo);

/**
   Tests if number of states and number of outputs in the models match.
   @return 0 for succes; -1 for error
   @param mo:           vector of models
   @param model_number: numbr of models */
extern int     ghmm_d_check_compatibility(ghmm_dmodel **mo, int model_number);

/**
   Produces a model, which generates the given sequence with probability 1.
   The model is a strict left-right model with one state for each element 
   in the sequence and the output in state i is the i-th value in the sequence 
   with probability 1. The model also has a final state, a state with no output.
   @return         pointer to the produced model 
   @param seq:      sequence
   @param seq_len:  length of the sequence
   @param anz_symb: number of symbols in the sequence
*/
extern ghmm_dmodel*  ghmm_d_generate_from_sequence(const int *seq, int seq_len, 
				     int anz_symb);

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
*/
extern ghmm_dseq *ghmm_d_generate_sequences(ghmm_dmodel* mo, int seed, int global_len,
				     long seq_number, int Tmax);

/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
extern double ghmm_d_likelihood(ghmm_dmodel *mo, ghmm_dseq *sq);


/** Computes probabilistic distance of two models
    @return the distance
    @param m0  model used for generating random output
    @param m  model to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal langth of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
double ghmm_d_prob_distance(ghmm_dmodel *m0, ghmm_dmodel *m, int maxT, int symmetric, int verbose);


/** 
   Allocates a new ghmm_d_background_distributions struct and assigns the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return     :               0 on success, -1 on error
   @param mo   :               one model
   @param cur  :               a id of a state
   @param times:               number of times the state cur is at least evaluated
*/
extern int ghmm_d_duration_apply(ghmm_dmodel* mo, int cur, int times);


/******** Reestimate Baum-Welch (reestimate.c) *******/
extern int ghmm_d_baum_welch(ghmm_dmodel *mo, ghmm_dseq *sq);
extern int ghmm_d_baum_welch_nstep(ghmm_dmodel *mo, ghmm_dseq *sq, int max_step, double likelihood_delta);

extern void ghmm_d_update_tied_groups(ghmm_dmodel *mo);

/**
  Baum Welch training for StateLabelHMMs
*/
extern int ghmm_dl_baum_welch(ghmm_dmodel *mo, ghmm_dseq *sq);

/**
  Just like reestimate_baum_welch_label, but you can limit
  the maximum number of steps
  */
extern int ghmm_dl_baum_welch_nstep(ghmm_dmodel *mo, ghmm_dseq *sq, int max_step, double likelihood_delta);


/*----------------------------------------------------------------------------*/
/**
   Trains the model with a set of annotated sequences till convergence using
   gradient descent.
   Model must not have silent states. (checked in Python wrapper)
   @return            0/-1 success/error
   @param mo:         pointer to a model
   @param sq:         struct of annotated sequences
   @param eta:        intial parameter eta (learning rate)
   @param no_steps    number of training steps
 */
extern int ghmm_dl_gradient_descent(ghmm_dmodel** mo, ghmm_dseq* sq, double eta, int no_steps);


/********  K-Best decoding (kbest.c) ********/
/**
   Calculates the most probable labeling for the given sequence in the given
   model using k-best decoding.
   @return array of labels (internal representation)
   @param mo:         pointer to a model
   @param o_seq:      output sequence (array of internal representation chars)
   @param seq_len:    length of output sequence
   */
extern int* ghmm_dl_kbest(ghmm_dmodel* mo, int* o_seq, int seq_len, int k, double* log_p);

/******* Viterbi (viterbi.c)*******/
/**
  Viterbi algorithm. Calculates the Viterbi path (the optimal path trough
  the model) and the Viterbi probability to a given model and a given 
  sequence.
  @return Viterbi path
  @param mo:    model
  @param o:     sequence
  @param len:   length of the sequence
  @param log_p: probability of the sequence in the Viterbi path
  */

// XXX use swig OUTPUT typemap
extern int * ghmm_d_viterbi(ghmm_dmodel *mo, int *o, int len, double *log_p);

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
extern double ghmm_d_viterbi_logp(ghmm_dmodel *mo, int *o, int len, int *state_seq);


/********  Discriminative Training (discrime.c) ********/
/**
   Trains two or more models to opimise the discrimination between the
   classes in the trainingset.
   The models must have the same topology. (checked)
   @return                 0/-1 success/error
   @param mo:              array of pointers to some models
   @param sqs:             array of annotated sequence sets
   @param noC:             number of classes
   @param max_steps:       maximum number of training steps for a class
   @param gradient:        if gradient == 0 try a closed form solution
                           otherwise a gradient descent
 */
extern int ghmm_d_discriminative(ghmm_dmodel** mo, ghmm_dseq** sqs, int noC, int max_steps,
			  int gradient);

/**
   Returns the value of teh in this discriminative training algorithm optimised
   function for a tupel of HMMs and sequencesets.
   @return                 value of funcion
   @param mo:              array of pointers to some models
   @param sqs:             array of annotated sequence sets
   @param noC:             number of classes
*/
extern double ghmm_d_discrim_performance(ghmm_dmodel** mo, ghmm_dseq** sqs, int noC);


/******* Forward , backward (foba.c) ******/

/** Forward-Algorithm.
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param mo:      model
  @param O:       sequence
  @param len:     length of sequence
  @param alpha:   alpha[t][i]
  @param scale:   scale factors
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int ghmm_d_forward (ghmm_dmodel * mo, const int *O, int length, double **alpha,
                    double *scale, double *log_p);

/** 
  Backward-Algorithm. 
  Calculates beta[t][i] given an integer sequence and a model. Scale factors 
  given as parameter (come from foba\_forward).
  @param mo:      model
  @param O:       sequence
  @param len:     length of sequence
  @param beta:    empty beta matrix
  @param scale:   scale factors
  @return 0 for success, -1 for error
  */
extern int ghmm_d_backward (ghmm_dmodel * mo, const int *O, int len, double **beta,
                     const double *scale);

/** 
  Termination of Backward-Algorithm. 
  Calculates Backward-probability given an integer sequence, a model and
  the beta matrix. Scale factors given as parameter (come from foba\_forward).
  @param mo:      pointer to the model
  @param O:       sequence
  @param len:     length of sequence
  @param beta:    beta matrix
  @param scale    scale factors
  @param log\_p:  log probability
  @return 0 for success, -1 for error
  */
extern int ghmm_d_backward_termination (ghmm_dmodel *mo, const int *O, int len,
				 double **beta, double *scale, double *log_p);


/**
  Calculation of  log( P(O|lambda) ). 
  Done by calling foba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed.
  @param  mo      model
  @param O        sequence
  @param len       length of sequence
  @param log\_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int ghmm_d_logp(ghmm_dmodel *mo, const int *O, int len, double *log_p);

/**
    Set transition from state 'i' to state 'j' to value 'prob'.
    NOTE: No internal checks - model might get broken if used carelessly.
    @param mo model
    @param i state index
    @param j state index
    @param prob probabilitys
    
*/
extern void ghmm_d_transition_set(ghmm_dmodel *mo, int i, int j, double prob);

/** Forward-Algorithm (lean version).
  Calculates log( P(O|lambda) ) for a given double sequence and a given model.
  @param mo       model
  @param O        sequence
  @param length: length of sequence
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
int ghmm_d_forward_lean(ghmm_dmodel *mo, const int *O, int len, double *log_p); 


/* Labeled HMMs */
extern int ghmm_dl_forward(ghmm_dmodel *mo, const int *O, const int *label, int len, double **alpha, double *scale, double *log_p);
extern int ghmm_dl_logp(ghmm_dmodel *mo, const int *O, const int *label, int len, double *log_p);
extern int ghmm_dl_backward(ghmm_dmodel *mo, const int *O, const int *label, int len, double **alpha, double *scale, double *log_p);


/** 
   Allocates a new ghmm_d_background_distributions struct and assigns the arguments to
   the respective fields. Note: The arguments need allocation outside of this
   function.
   
   @return    :               new pointer to a ghmm_d_background_distributions struct
   @param n   :               number of distributions
   @param order:              orders of the distribtions
   @param B:                  matrix of distribution parameters
*/
ghmm_d_background_distributions *ghmm_d_background_alloc(int n,int m, int *orders, double **B);

extern ghmm_d_background_distributions *ghmm_d_background_copy(ghmm_d_background_distributions *bg);
extern int ghmm_d_background_free(ghmm_d_background_distributions *bg);


/**
    Uses ighmm_cvector_normalize in vector.h
    Normalizes the transition and output probs for each state
    in the given model
    @return 0 if normalization went through
    @param mo: model to be normalized

*/
extern int ghmm_d_normalize(ghmm_dmodel* mo);

/**
   Add a specific level of noise to the model parameters
   @return     :        -1 on error, 0 else
   @param mo   :        a pointer to the model
   @param level:        amount of noise to use,
                        a noise level of 0.0 doesn't change the model
   @param seed :        seed for ramdom number generator
*/
extern int ghmm_d_add_noise(ghmm_dmodel* mo, double level, int seed);

/**
   Apply the background distributions to the emission probabilities of states of
   the model which have one specified (background_id[state_id] != kNoBackgroundDistribution).

   @return    :                -1 on error, 0 else
   @param mo  :                a pointer to the model
   @param background_weight:   a parameter controlling the weight given to the
                               background. Note, should be between 0 and 1.
*/
extern int ghmm_d_background_apply(ghmm_dmodel *mo, double* background_weight);




%inline%{
	
  /* allocation of an empty ghmm_dmodel struct */
  ghmm_dmodel *new_model() {
     return (struct ghmm_dmodel *) calloc(1, sizeof(struct ghmm_dmodel));    
  }

   
  /* allocation of an array of ghmm_dstate structs*/
  ghmm_dstate *arraystate(int size) {
    return (ghmm_dstate *) calloc(size*sizeof(ghmm_dstate), 1);
  }
  
  /* extract pointer to a ghmm_dstate */
  ghmm_dstate *get_stateptr(ghmm_dstate *ary, int index) { return ary + index; }
  

  void call_model_print(char *filename, ghmm_dmodel *mo) {
    FILE *fp=fopen(filename, "a");
    if (fp == NULL) {
      fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);    
    } 
    else {
      ghmm_d_print(fp, mo);
      fclose(fp);
    }
  }
  
  void call_ghmm_d_free (ghmm_dmodel *mo) {ghmm_d_free(&mo);}

  ghmm_dmodel *get_model_ptr(ghmm_dmodel **mo, int index) {return mo[index];}
    
  ghmm_dmodel * * cast_model_ptr(ghmm_dmodel *mo){
    ghmm_dmodel * * result = &mo;
    return result;
  }   

  
  
%}

/*=============================================================================================
  =============================== labeled models (model.c)  ============================================== */

extern ghmm_dseq *ghmm_dl_generate_sequences(ghmm_dmodel* mo, int seed, int global_len, long seq_number, int Tmax);

/*=============================================================================================
  =============================== sdmodel.c  ============================================== */

/** @name ghmm_dsstate
    The basic structure, keeps all parameters that belong to a ghmm_dsstate. 
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
 /** ghmm_dsmodel includes discrete model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
  int cos;
  /** Vector of the states */
  ghmm_dsstate *s; 
  /** Prior for the a priori probability for the model.
      A value of -1 indicates that no prior is defined. */
  double prior;
  /** pointer to class function
   */
  int (*get_class)(const int*,int);

 /** Contains bit flags for various model extensions such as
     kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
 */
  int model_type;
  
  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
  int* silent; /*AS*/

};
typedef struct ghmm_dsmodel ghmm_dsmodel;

/** Frees the memory of a model.
    @return 0 for succes; -1 for error
    @param mo:  pointer to a ghmm_dmodel */
int     ghmm_ds_free (ghmm_dsmodel **mo);

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
*/
ghmm_dseq *ghmm_ds_generate_sequences (ghmm_dsmodel* mo, int seed, int global_len,
				     long seq_number, int Tmax);


/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
extern double ghmm_ds_likelihood (ghmm_dsmodel *mo, ghmm_dseq *sq);

/** 
    Computes log likelihood for all sequence of
    seq_w * log( P ( O|lambda ) ). If a sequence can't be generated by smo
    error cost of seq_w * PRENALTY_LOGP are imposed.
   @return       n: number of evaluated sequences, -1: error
   @param smo    ghmm_cmodel
   @param sqd    sequence struct
   @param log\_p array of evaluated likelihoods
*/
extern int ghmm_c_individual_likelihoods(ghmm_cmodel *smo, ghmm_cseq *sqd, double *log_ps);

/******* Viterbi for switching discrete ghmm_dmodel (sdviterbi.c) *******/
int *ghmm_ds_viterbi (ghmm_dsmodel *mo, int *o, int len, double *log_p);

/******* Forward-Algorithm. (sdfoba.c) *******
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param mo       model
  @param O        sequence
  @param length: length of sequence
  @param alpha:  alpha[t][i]
  @param scale:   a reference for double type, scale factors
  @param log\_p:  a reference for double type, log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int ghmm_ds_forward (ghmm_dsmodel *mo, const int *O, int len, double **alpha, 
		 double *scale, double *log_p); 


%inline%{
  
  void call_sdmodel_free (ghmm_dsmodel * sdm)	{ghmm_ds_free(&sdm);}
	
  /* allocation of an array of ghmm_dsstate structs*/	
  ghmm_dsstate *arraysdstate(int size) {
    return (ghmm_dsstate *) malloc(size*sizeof(ghmm_dsstate));
  }	
  
  /* extract pointer to a ghmm_dsstate  */
  ghmm_dsstate *get_sdstateptr(ghmm_dsstate *ary, int index) { return ary + index; }
  
%}


/*=============================================================================================
  =============================== smodel.c  ============================================== */

   typedef enum {
    normal,
    normal_right,
    normal_approx,
    normal_left,
    uniform,
    density_number
  } ghmm_density_t;

/** @name ghmm_cstate
    Structure for one state.
*/
struct ghmm_cstate {
  /** Maximun number of output densities per state */
  int M;
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
  /** value where the normal is truncated (only for truncated distributions) */
    double *a;
  /** flag for fixation of parameter. If fix = 1 do not change parameters of
      output functions, if fix = 0 do normal training. Default is 0. */
  int fix;
  /** Flag for density function for each component of the mixture
      0: normal density, 1: truncated normal (right side) 
      density, 2: approximated normal density, 3: truncated normal (left side)
      4: uniform distribution */
  ghmm_density_t *density;  
  /**  array of flags for fixing mixture components in the reestimation
        fix[i] = 1 means mu and sigma of component i are fixed.  **/
  int *mixture_fix;
};
typedef struct ghmm_cstate ghmm_cstate;

struct ghmm_cmodel;

struct ghmm_c_class_change_context{

    /* Names of class change module/function (for python callback) */
    char* python_module;
    char* python_function;
       
    /* index of current sequence */ 
    int k;

    /** pointer to class function */
    int (*get_class)(struct ghmm_cmodel*,double*,int,int);
    
    /* space for any data necessary for class switch, USER is RESPONSIBLE */
    void* user_data;
};
typedef struct ghmm_c_class_change_context ghmm_c_class_change_context;

/** @name ghmm_cmodel
    continous HMM    
*/
struct ghmm_cmodel{
  /** Number of states */
  int N;
  /** Maximun number of output densities per state */
  int M;
  /** ghmm_cmodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
  int cos;
  /** prior for a priori prob. of the model. -1 means no prior specified (all
      models have equal prob. a priori. */
  double prior;
  /** All states of the model. Transition probs are part of the states. */
  ghmm_cstate *s; 
 
  /* pointer to a ghmm_c_class_change_context struct necessary for multiple transition
   classes */
  ghmm_c_class_change_context *class_change;  
  
 };
typedef struct ghmm_cmodel ghmm_cmodel;


extern int ghmm_c_class_change_alloc(ghmm_cmodel *smo);

/** Free memory ghmm_cmodel 
    @return 0: success, -1: error
    @param smo  address of a pointer to a ghmm_cmodel */
extern int     ghmm_c_free(ghmm_cmodel **smo);

#ifdef GHMM_OBSOLETE
/** Reads an ascii file with specifications for one or more models.
    All parameters in matrix or vector form.
    This is necessary whenever an initial model is needed (e.g. 
    training) or sequences have to be generated from trained models.
    For each ghmm_cmodel block smodel\_read\_block() is called.
   @return vector of read ghmm_cmodels
   @param filename   input ascii file
   @param smo_number  number of read ghmm_cmodels */
extern ghmm_cmodel** ghmm_c_read(const char *filename, int *smo_number);
#endif /* GHMM_OBSOLETE */

/**
   Copies one ghmm_cmodel. Memory alloc is here.
   @return pointer to ghmm_cmodel copy
   @param smo   model to be copied  */
extern ghmm_cmodel*  ghmm_c_copy(const ghmm_cmodel *smo);

/**
   Produces sequences to a given model. All memory that is needed for the 
    sequences is allocated inside the function. It is possible to define
    the length of the sequences global (global_len > 0) or it can be set 
    inside the function, when a final state in the model is reach (a state
    with no output). If the model has no final state, the sequences will
    have length MAX_SEQ_LEN.
*/
extern ghmm_cseq *ghmm_c_generate_sequences(ghmm_cmodel* smo, int seed, int global_len,
					       long seq_number, long label, int Tmax);

/** Computes the density of one symbol (omega) in a given state (sums over
    all output components
    @return calculated density
    @param smo   model
    @param state state 
    @param omega given symbol
*/
extern double ghmm_c_calc_b (ghmm_cmodel * smo, int ghmm_dstate, double omega);

/** Computes probabilistic distance of two models
    @return the distance
    @param cm0  ghmm_cmodel used for generating random output
    @param cm   ghmm_cmodel to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal length of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
extern double ghmm_c_prob_distance(ghmm_cmodel *cm0, ghmm_cmodel *cm, int maxT, int symmetric, int verbose);


/** Forward-Algorithm.  (sfoba.c)

  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param T        length of sequence
  @param b        matrix with precalculated output probabilities. May be NULL
  @param alpha    alpha[t][i]
  @param scale    scale factors
  @param log_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int ghmm_c_forward(ghmm_cmodel *smo, double *O, int T, double ***b, 
		  double **alpha, double *scale, double *log_p);

/** 
  Backward-Algorithm.  (sfoba.c)
  Calculates beta[t][i] given a double sequence and a model. Scale factors 
  given as parameter (come from sfoba\_forward).
  @param smo      model
  @param O          sequence
  @param T        length of sequence
  @param b        matrix with precalculated output probabilities. May be NULL
  @param beta     beta[t][i]
  @param scale    scale factors
  @return 0 for success, -1 for error
  */
extern int ghmm_c_backward(ghmm_cmodel *smo, double *O, int T, double ***b,
		   double **beta, const double *scale);

/**
  Calculation of  log( P(O|lambda) ).  (sfoba.c)s
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed, alpha matrix is allocated with
  ighmm_cmatrix_stat_alloc
  @param smo      model
  @param O        sequence
  @param T         length of sequence
  @param log_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int ghmm_c_logp(ghmm_cmodel *smo, double *O, int T, double *log_p);


/** 
    Computes sum over all sequence of
    seq_w * log( P ( O|lambda ) ). If a sequence can't be generated by smo
    error cost of seq_w * PRENALTY_LOGP are imposed.
   @return       n: number of evaluated sequences, -1: error
   @param smo    model
   @param sqd    sequence struct
   @param log\_p  evaluated log likelihood
*/

extern ghmm_cmodel *smodel_alloc_fill(int N, int M, int cos, double prior, int density);

//extern void smodel_set_densityvector(ghmm_cmodel *smo, int i, int density);

extern void smodel_set_pivector(ghmm_cmodel *smo, int i, double piv);

extern void smodel_set_fixvector(ghmm_cmodel *smo, int i, double fixv);

extern void smodel_set_transition(ghmm_cmodel *smo, int i, int j, int cos, double prob);

extern double smodel_get_transition(ghmm_cmodel *smo, int i, int j, int cos);

extern void smodel_set_mean(ghmm_cmodel *smo, int i, double *mu);

extern void smodel_set_variance(ghmm_cmodel *smo, int i, double *variance);

// extern void call_smodel_print(char *filename, ghmm_cmodel *smo);

extern int ghmm_c_likelihood(ghmm_cmodel *smo, ghmm_cseq *sqd, double *log_p);

extern int smodel_sorted_individual_likelihoods(ghmm_cmodel *smo, ghmm_cseq *sqd, double *log_ps, int *seq_rank);

/** Viterbi for ghmm_cmodel (sviterbi.c)

   Viterbi algorithm: calculation of the viterbi path (best possible
   state sequenz for a given sequenz and a given model (smo)). Also
   calculates logp according to this path, the matrices in the local_store
   struct are allocated using stat_matrix_d_alloc.
  @return        Viterbi-path 
  @param smo     model
  @param o       double-sequence
  @param T       sequence length
  @param log_p   log(p) of the sequence using the vitberbi path
  */
extern int *ghmm_c_viterbi(ghmm_cmodel *smo, double *o, int T, double *log_p);


extern int cp_class_change(ghmm_cmodel *smo, double *seq, int k, int t);
extern void setSwitchingFunction( ghmm_cmodel *smd );

extern int python_class_change(ghmm_cmodel *smo, int *seq, int k, int t);
extern void setPythonSwitching( ghmm_cmodel *smd, char* python_module, char* python_function);


/* TEST */
/* static PyObject *pyCallBack = NULL;*/

extern void setPythonCallback(ghmm_cmodel *smo, PyObject *py_cb);
extern int executePythonCallback(ghmm_cmodel* smo, double *seq, int k, int t);



%inline %{
 
  void set_density(ghmm_cstate *state, int m, int value){
    state->density[m] = (ghmm_density_t)value;
  }

  int blatestbla(ghmm_cmodel* smo,double *seq ,int k, int t){
      return smo->class_change->get_class(smo,seq,k,t);
  }

  /* allocation of an empty ghmm_cmodel struct */
  ghmm_cmodel *new_smodel() {
     return (struct ghmm_cmodel *)(struct ghmm_cmodel *) calloc(1, sizeof(ghmm_cmodel));
  }

  /** allocation of an empty ghmm_cmodel struct */
  ghmm_density_t *arraydensity(int size) {
     return (ghmm_density_t *) malloc(size*sizeof(ghmm_density_t));
  }

  int get_density(ghmm_cstate *state, int m){
    return (int)(state->density[m]);
  }

  /* array of ghmm_cstate structs */
  ghmm_cstate *arraysstate(int size) {
    return (ghmm_cstate *) malloc(size*sizeof(ghmm_cstate));
  }

  ghmm_cstate *get_sstate_ptr(ghmm_cstate *states, int k) {
    return &(states[k]);
  }

  void call_smodel_free (ghmm_cmodel *smo) {ghmm_c_free(&smo);}

  void free_smodel_array (ghmm_cmodel **smo) {if (smo) {free(smo);} }

  void smodel_print_stdout(ghmm_cmodel *smo) {
    ghmm_c_print(stdout, smo);
  }

  /* extract pointer to ghmm_cstate  */
  ghmm_cstate *get_sstate(ghmm_cmodel *smo, int k) {
    return &(smo->s[k]);
  }

  /* creation and assessor functions for an array of ghmm_cmodels */
  ghmm_cmodel **smodel_array(int size){
     return (ghmm_cmodel**) malloc(size * sizeof(ghmm_cmodel*));
  }

  ghmm_cmodel *get_smodel_ptr(ghmm_cmodel **smo, int index) { return smo[index]; }

  void set_smodel_ptr(ghmm_cmodel **smo_array ,ghmm_cmodel *smo, int index) { smo_array[index] = smo; }

  ghmm_cmodel **cast_smodel_ptr(ghmm_cmodel *smo){
     ghmm_cmodel** res = (ghmm_cmodel**) malloc(sizeof(ghmm_cmodel*));
     res[0] = smo;
     return res;
  }

  // write a ghmm_cmodel to a file
  void call_smodel_print(char *filename, ghmm_cmodel *smo) {
  FILE *fp=fopen(filename, "a");
  if (fp == NULL) {
    fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);
  } else {
    ghmm_c_print(fp, smo);
    fclose(fp);
  }
}

%}



/* =============================================================================================
   ============================== scluster.c  ================================================== */
#ifdef GHMM_OBSOLETE
/**
   Cluster structure: All models and sequences. */
struct scluster_t{
  /** 
  Vector of SHMMs pointers */
  ghmm_cmodel **smo;
  /** 
    Vector of ghmm_cseq pointers; to store the sequences, that  belong to the models */
  ghmm_cseq **smo_seq;
  /** 
    Number of models to read in */
  int smo_number;
  /** 
    Number of sequences for each model */
  long *seq_counter;
  /** 
    log(P) for the model */
  double *smo_Z_MD;
  /** a posteriori probability for the models to calculate the objective
      fucntion in case of a MAW classificator. Is calculated using ghmm_smap_bayes */
  double *smo_Z_MAW;
}; typedef struct scluster_t scluster_t;


/**
   Frees the memory allocated for a scluster_t struct.
   @return 0 for success; -1 for error
   @param scl pointer to scl struct
*/
extern int ghmm_scluster_t_free(scluster_t *scl);


/**
   Creates random labels for a vector of sequences
   @return 0 for success; -1 for error
   @param sqd vector of sequences
   @param smo_number number of models (needed to determine the interval
   for the random numbers)
*/
extern int ghmm_scluster_random_labels(ghmm_cseq *sqd, int smo_number);

/**
   Calculates the logarithmic probability of sequences for a model.
   @param cs sequences and model 
 */
extern void ghmm_scluster_prob(ghmm_c_baum_welch_context *cs);

/** 
    Determines form an already calculated probability matrix, which model 
    fits best to a certain sequence. 
    @return index of best model if success, otherwize -1
    @param cl cluster 
    @param seq_id ID of the sequence in question
    @param all_log_p matrix containing the probability of each sequence
                     for each model
    @param log_p the probability of the sequence in question for the 
                 best fitting model
*/
extern int ghmm_scluster_best_model(scluster_t *cl, long seq_id, double **all_log_p,
			double *log_p);

/**
   Updates the cluster with additional sequences.
   @return 0 for success; -1 for error
   @param cl cluster to update
   @param sqd sequences to update the cluster with
 */
extern int ghmm_scluster_update(scluster_t *cl, ghmm_cseq *sqd);

/**
   Prints the input vector for ghmm_scluster_hmm
 */
void ghmm_scluster_print_likelihood(FILE *outfile, scluster_t *cl);
	

/**
   Writes out the final models.
   @return 0 for success; -1 for error
   @param cl cluster of models to write
   @param sqd
   @param outfile output file
   @param out_filename name of the output file
 */
int ghmm_scluster_out(scluster_t *cl, ghmm_cseq *sqd, FILE *outfile, char *argv[]);

/** Calculates the aposteriori prob. $\log(p(\lambda_best | O[seq\_id]))$, 
    where $\lambda_best$ is the model with highest apost. prob.
    @return 0 for success; -1 for error
    @param cl cluster
    @param sqd the sequence in question
    @param seq_id the ID of the sequence
    @param lob_apo the results
*/
extern int  ghmm_scluster_log_aposteriori(scluster_t *cl, ghmm_cseq *sqd, int seq_id, double *log_apo);

/**
   Avoids empty models going out as outputs by assigning them a random 
   sequence. This may lead to a produce of a new empty model - therefore
   change out sequences until a non empty model is found. (Quit after 100 
   iterations to avoid an infinite loop). 
   @return 0 for success; -1 for error
   @param sqd sequences to generate the model from
   @param cl cluster for the models
 */
extern int ghmm_scluster_avoid_empty_smodel(ghmm_cseq *sqd, scluster_t *cl);

/**
   Makes a cluster and reestimates the HMMs.
   @return 0 for success; -1 for error
   @param argv vector of input files, one with sequences, one with models, 
   one for output and one with labels for the sequences - in this order.
 */
extern int ghmm_scluster_hmm(char *argv[]);


%inline %{
  void scluster_printl_stdout(scluster_t *scl) {
    /*ghmm_scluster_print_likelihood(stdout, scl);*/
  }
  
%}

#endif /* GHMM_OBSOLETE */


/*=============================================================================================
  =============================== sreestimate.c  ============================================== */

/** @name struct ghmm_c_baum_welch_context
    structure that combines a continuous model (smo) and an integer
    sequence struct. Is used by sreestimate\_baum\_welch for 
    parameter reestimation.
 */
struct ghmm_c_baum_welch_context {
  /** pointer of continuous model*/
  ghmm_cmodel *smo;
  /** sequence\_d\__t pointer */
  ghmm_cseq *sqd;
  /** calculated log likelihood */
  double* logp;
  /** leave reestimation loop if diff. between successive logp values 
      is smaller than eps */
  double eps;
  /** max. no of iterations */
  int max_iter;
}; 
typedef struct ghmm_c_baum_welch_context ghmm_c_baum_welch_contex;


/**
  Baum-Welch Algorithm for SHMMs.
  Training of model parameter with multiple double sequences (incl. scaling).
  New parameters set directly in hmm (no storage of previous values!).
  @return            0/-1 success/error
  @param cs         initial model and train sequences
  */
extern int ghmm_c_baum_welch(ghmm_c_baum_welch_context *cs);

/********** Here comes all the Pair HMM stuff  **********/

/********** Pair HMM ghmm_dpseq (psequence.c) **********/

struct ghmm_dpseq {
  /** for each alphabet in model->number_of_alphabets there is one int seq **/
  int ** seq;
  /** number of alphabets (same as in model) **/
  int number_of_alphabets;
  /** for each sequence position there are also double values (e.g) Ka **/
  double ** d_value;
  /** number of continous sequences **/
  int number_of_d_seqs;
  /** length of the sequence **/
  int length;
};

typedef struct ghmm_dpseq ghmm_dpseq;

extern ghmm_dpseq * ghmm_dpseq_init(int length, int number_of_alphabets, int number_of_d_seqs);

extern int ghmm_dpseq_free(ghmm_dpseq * seq, int number_of_alphabets, int number_of_d_seqs);

extern void ghmm_dpseq_set_discrete(ghmm_dpseq * seq_pointer, int index, int * sequence);

extern void ghmm_dpseq_set_continuous(ghmm_dpseq * seq_pointer, int index, double * sequence);

extern int * ghmm_dpseq_get_discrete(ghmm_dpseq * seq_pointer, int index);

extern double * ghmm_dpseq_get_continuous(ghmm_dpseq * seq_pointer, int index);

/********** End of Pair HMM ghmm_dpseq (psequence.c) **********/

/********** Pair HMM ghmm_dmodel (pmodel.c) **********/
struct ghmm_dp_class_change_context{

    /* Names of class change module/function (for python callback) */
    char* python_module;
    char* python_function;
    
    /** pointer to class function called with seq X, Y and resp indices */
    int (*get_class)(struct ghmm_dpmodel*, ghmm_dpseq*, ghmm_dpseq*, int, int,void*);
    
    /* space for any data necessary for class switch, USER is RESPONSIBLE */
    void* user_data;
};
typedef struct ghmm_dp_class_change_context ghmm_dp_class_change_context;

struct ghmm_dpstate {
  /** Initial probability */ 
  double pi;
  /** Log of the initial probability */
  double log_pi;
  /** Output probability */
  double *b;
  int order;
  
  /** IDs of the following states */ 
  int *out_id;  
  /** IDs of the previous states */    
  int *in_id;

  /** transition probs to successor states. (for each transition class)*/
  double **out_a; 
  /** transition probs from predecessor states. (for each transition class)*/ 
  double **in_a;
  /** number of transition classes in this state **/
  int kclasses;
  /** pointer to class function   */
  ghmm_dp_class_change_context *class_change; 
  /** int (*get_class)(int*,int); */
  
  /** Transition probability to a successor 
      double *out_a; */
  /** Transition probablity to a precursor 
      double *in_a;*/

  /** Number of successor states */     
  int out_states; 
  /** Number of precursor states */
  int in_states;  
  /** if fix == 1 --> b stays fix during the training */
  int fix;

  int label;
  
  /** EXTENSIONS for Pair HMMs **/
  /** read offset_x characters from sequence X **/
  int offset_x;
  /** read offset_y characters from sequence Y **/
  int offset_y;
  /** which emission alphabet **/
  int alphabet;
};
typedef struct ghmm_dpstate ghmm_dpstate;

/** @name ghmm_dpmodel
    The complete HMM. Contains all parameters, that define a HMM.
*/
struct ghmm_dpmodel {
  /** Number of states */
  int N;
  /** Number of outputs */   
  int M;   
  /** Vector of the states */
  ghmm_dpstate *s; 
  /** The a priori probability for the model.
      A value of -1 indicates that no prior is defined. 
      Note: this is not to be confused with priors on emission
      distributions*/
  double prior;

  /* contains a arbitrary name for the model */
  char* name;
  
   /** Contains bit flags for varios model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
  int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
  int* silent; /*AS*/

  /** Int variable for the maximum level of higher order emissions */
  int maxorder;
  /** saves the history of emissions as int, 
      the nth-last emission is (emission_history * |alphabet|^n+1) % |alphabet|
      see ...*/
  int emission_history;

  /** Flag variables for each state indicating whether the states emissions
      are tied to another state. Groups of tied states are represented
      by their tie group leader (the lowest numbered member of the group).
      
      tied_to[s] == kUntied  : s is not a tied state
      
      tied_to[s] == s        : s is a tie group leader

      tied_to[t] == s        : t is tied to state s

      Note: tied_to != NULL iff (model_type & kTiedEmissions) == 1  */
  int* tied_to; 
  
  /** Note: State store order information of the emissions.
      Classic HMMS have emission order 0, that is the emission probability
      is conditioned only on the state emitting the symbol.

      For higher order emissions, the emission are conditioned on the state s
      as well as the previous emission_order[s] observed symbols.

      The emissions are stored in the state's usual double* b. The order is
      set state.order.

      Note: state.order != NULL iff (model_type & kHigherOrderEmissions) == 1  */
  
  /** bp is a pointer to a ghmm_d_background_distributions structure,
      which holds (essentially) an array of background distributions 
      (which are just vectors of floating point numbers like state.b).

      For each state the array background_id indicates which of the background
      distributions to use in parameter estimation. A value of kNoBackgroundDistribution
      indicates that none should be used.


      Note: background_id != NULL iff (model_type & kHasBackgroundDistributions) == 1  */
  int *background_id;
  ghmm_d_background_distributions* bp; 

  /** (WR) added these variables for topological ordering of silent states 
      Condition: topo_order != NULL iff (model_type & kSilentStates) == 1
   */
  int* topo_order; 
  int  topo_order_length;

  /** EXTENSIONS for Pair HMMs **/
  /** total number of alphabets **/
  int number_of_alphabets;
  /** list of sizes of the alphabets **/
  int * size_of_alphabet;
  /** number of double sequences to modify the transition classes */
  int number_of_d_seqs;
  /** maximal offset in sequence X (for the viterbi lookback matrix) **/
  int max_offset_x;
  /** maximal offset in sequence Y (for the viterbi lookback matrix) **/
  int max_offset_y;

};
typedef struct ghmm_dpmodel ghmm_dpmodel;

extern void ghmm_dp_state_clean(ghmm_dpstate *my_state);

extern ghmm_dpmodel * ghmm_dp_init();

ghmm_dp_class_change_context * ghmm_dp_init_class_change();

extern int ghmm_dp_free(ghmm_dpmodel *mo);

%inline %{ ghmm_dpstate * get_pstateptr(ghmm_dpstate * ary, int index) {return &(ary[index]);} %}

extern int ghmm_dp_pair(int symbol_x, int symbol_y, int alphabet_size, int off_x, int off_y);

extern int ghmm_dp_default_transition_class(ghmm_dpmodel * mo, ghmm_dpseq * X, ghmm_dpseq * Y, int index_x, int index_y, void * user_data);

extern void ghmm_dp_set_to_default_transition_class(ghmm_dp_class_change_context * pccc);

/********** End of Pair HMM ghmm_dmodel (pmodel.c) **********/

/********  Pair HMM viterbi (pviterbi.c) *********/

extern int *ghmm_dp_viterbi(ghmm_dpmodel *mo, ghmm_dpseq * X, ghmm_dpseq * Y, double *log_p, int * length);

int *ghmm_dp_viterbi_variable_tb(ghmm_dpmodel *mo, ghmm_dpseq * X, ghmm_dpseq * Y, double *log_p, int *path_length, int start_traceback_with);

double ghmm_dp_viterbi_logp(ghmm_dpmodel *mo, ghmm_dpseq * X, ghmm_dpseq * Y, int *state_seq, int state_seq_len);

/********  End of Pair HMM viterbi (pviterbi.c) *********/

/********  Pair HMM viterbi linear space (pviterbi_propagate.c) *********/

extern int * ghmm_dp_viterbi_propagate (ghmm_dpmodel *mo, ghmm_dpseq * X, ghmm_dpseq * Y,
				 double *log_p, int *path_length, double max_size);

extern int * ghmm_dp_viterbi_propagate_segment (ghmm_dpmodel *mo, ghmm_dpseq * X, ghmm_dpseq * Y,
					 double *log_p, int *path_length, double max_size,
					 int start_x, int start_y, int stop_x, int stop_y,
					 int start_state, int stop_state, double start_log_p,
					 double stop_log_p);


/********  End of Pair HMM viterbi linear space (pviterbi_propagate.c) *******/

/********  Pair HMM class change (pclasschange.c) *********/

struct threshold_user_data {
  /** which double value in myseq **/
  int seq_index;
  /** cut off value **/
  double threshold;
  /** add this to the index in sequence X **/
  int offset_x;
  /** add this to the index in sequence Y **/
  int offset_y;
};
typedef struct threshold_user_data threshold_user_data;

extern int gt_sum(ghmm_dpmodel * mo, ghmm_dpseq * X, ghmm_dpseq * Y, int index_x, int index_y, void * user_data);

extern int lt_sum(ghmm_dpmodel * mo, ghmm_dpseq * X, ghmm_dpseq * Y, int index_x, int index_y, void * user_data);

extern int boolean_and(ghmm_dpmodel * mo, ghmm_dpseq * X, ghmm_dpseq * Y, int index_x, int index_y, void * user_data);

extern int boolean_or(ghmm_dpmodel * mo, ghmm_dpseq * X, ghmm_dpseq * Y, int index_x, int index_y, void * user_data);

extern void set_to_lt_sum(ghmm_dp_class_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y);

extern void set_to_gt_sum(ghmm_dp_class_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y);

void set_to_boolean_and(ghmm_dp_class_change_context * pccc, int seq_index, int offset_x, int offset_y);

void set_to_boolean_or(ghmm_dp_class_change_context * pccc, int seq_index, int offset_x, int offset_y);

/********  End of Pair HMM class change (pclasschange.c) *********/

/********** utility functions for Pair HMMs **********/
%inline%{
  
  /* allocation of an array of ghmm_dpstate structs*/
  ghmm_dpstate *arraypstate(int size) {
    return (ghmm_dpstate *) malloc(size*sizeof(ghmm_dpstate));
  }
  
%}

%inline%{
  
  /************ create and manipulate an array of pointers to ghmm_c_baum_welch_context structs **************/
  ghmm_c_baum_welch_context *smosqd_t_array(int size) {
     return (ghmm_c_baum_welch_context *) malloc(size*sizeof(ghmm_c_baum_welch_context));
  }

  void set_smosq_t_smo(ghmm_c_baum_welch_context *cs, ghmm_cmodel *smo, int index){
	  cs[index].smo = smo;
  }	  
  
  ghmm_c_baum_welch_context get_smosqd_t_ptr(ghmm_c_baum_welch_context *cs, int i) {return cs[i];}
  
  void free_smosqd_t(ghmm_c_baum_welch_context *s){
	  if(s) {free(s);}	
  }	  
  		
%}	


/*=============================================================================================
  =============================== miscellanous functions ====================================== */

%inline %{
   
  // Some array helpers	
  /* Create any sort of int[size] array */
  int *int_array(int size) {
     return (int *) malloc(size*sizeof(int));
  }

  void set_arrayint(int  *ary, int index, int value) {
    ary[index] = value;
  }

  int  get_arrayint(int  *ary, int index) { return ary[index]; }

  void free_arrayi(int *pt ) { 
          free (pt); 
          (pt) = NULL;
  }

  /************ Create and access double[size] arrays ************/
 
   double *double_array(int size) {
     return (double *) malloc(size*sizeof(double));
  }

  void set_arrayd(double *ary, int index, double value) {
    ary[index] = value;
  }

  double get_arrayd(double *ary, int index) { return ary[index]; }
  
  void free_arrayd(double *pt) {free (pt); (pt) = NULL;}
  
  /************  Create and access sort of long[size] arrays **************/
    
  long *long_array(int size) {
     return (long *) malloc(size*sizeof(long));
  }

  void set_arrayl(long *ary, int index, long value) {
    ary[index] = value;
  }
  
  double get_arrayl(long *ary, int index) { return ary[index]; }
  
  void free_arrayl(long *ary) {free (ary); (ary) = NULL;}
  
  /*********** Create and access char** arrays ***********/
  char **char_array(int size) {
     return (char **) malloc(size*sizeof(char*));
  }

  void set_arraychar(char** ary, int index, char* value) {
    ary[index] = value;
  }

  char *get_arraychar(char** ary, int index) { return ary[index]; }

  void free_arraychar (char * * ary) {free (ary);}

  /********** Create and access ghmm_dmodel arrays  ****************************/
ghmm_dmodel * * modelarray_alloc (int size) {
  ghmm_dmodel * * retval;
  retval = calloc (size, sizeof(ghmm_dmodel*));
  return retval;
}

void modelarray_free (ghmm_dmodel ** mos) {
  free (mos);
}

void modelarray_setptr (ghmm_dmodel ** mos, ghmm_dmodel * mo, int pos) {mos[pos] = mo;}

ghmm_dmodel * modelarray_getptr (ghmm_dmodel ** mos, int pos) {return mos[pos];}


  /********** Create and access sequence* arrays  *************************/
ghmm_dseq * * seqarray_alloc (int size) {
  ghmm_dseq * * retval;
  retval = calloc (size, sizeof(ghmm_dseq*));
  return retval;
}

void seqarray_free (ghmm_dseq ** seqs) {
  free (seqs);
}

void seqarray_setptr (ghmm_dseq ** seqs, ghmm_dseq * seq, int pos)
{seqs[pos] = seq;}

ghmm_dseq * seqarray_getptr (ghmm_dseq ** seqs, int pos)
{return seqs[pos];}
   

  /************  Create and access double[size1][size2] arrays ************/

  double * * double_2d_array(int rows, int cols) {

    double **matrix;
    int i;

    matrix = calloc(rows, sizeof(double*));
    if (matrix)
      for (i=0; i<rows; i++) {
        matrix[i] = calloc(cols, sizeof(double));
        /* clean up and abort if an allocation fails */
        if (!matrix[i]) {
          for (i=i-1; i>=0; --i)
            free(matrix[i]);
          free(matrix);
          matrix = NULL;
          break;
        }
      }
    return matrix;
  }

  double * * double_2d_array_nocols(int rows) {
    return (double **) malloc(rows*sizeof(double*));
  }

  void set_2d_arrayd_col(double **ary, int index, double *col) {
    ary[index] = col;
  }

  void set_2d_arrayd(double **ary, int index1,int index2, double value) {
    ary[index1][index2] = value;
  }

  double get_2d_arrayd(double **ary, int index1, int index2) {
    return ary[index1][index2];
  }

  double *get_col_pointer_d(double **ary, int index) {
    return ary[index];
  }

  void double_2d_print(double **ary, int row, int col) {
    int i,j;
    printf("(%dx%d) Matrix", row, col);
    for (i =0;i<row;i++) {
      printf("\n");
      for (j =0;j<col;j++){
        printf("%f ",ary[i][j]);
      }
    }
    printf("\n");
  }

  void free_2darrayd(double **pt, int row) {
    unsigned int i;

    if (pt)
      for (i=0; i<row; i++)
        free(pt[i]);
    free (pt);
    return;
  }

  /************  Create and access int[size1][size2] arrays ************/

  int * * int_2d_array_nocols(int rows) {
    int * * matrix = malloc(rows*sizeof(int*));
    return matrix;
  }

  void set_2d_arrayint_col(int **ary, int index, int *col) {
    ary[index] = col;
  }

  int *get_col_pointer_int(int **ary, int index) {
    return ary[index];
  }

  void set_2d_arrayint(int **ary, int index1,int index2, int value) {
    ary[index1][index2] = value;
  }

  /* Get two dimensional array entry */
  int  get_2d_arrayint (int **ary, int index1, int index2) {
    return ary[index1][index2];
  }


  /**************** generalized deallocation *******************/
  void freearray(void *pt) {
    free(pt);
  }

%}


/******************* typemap ********************************/

%typemap(memberin) ghmm_cseq * {
  ghmm_cseq *src_pt;
  ghmm_cseq *dest_pt;

  dest_pt = (ghmm_cseq *) malloc(sizeof(ghmm_cseq));
  src_pt  = $input;

  dest_pt.seq        = src_pt.seq;
  dest_pt.seq_len    = src_pt.seq_len;
#ifdef GHMM_OBSOLETE
  dest_pt.seq_label  = src_pt.seq_label;
#endif /* GHMM_OBSOLETE */
  dest_pt.seq_id     = src_pt.seq_id;
  dest_pt.seq_w      = src_pt.seq_w;
  dest_pt.seq_number = src_pt.seq_number;
  dest_pt.total_w    = src_pt.total_w;

  $1 = dest_pt;
}

%typemap(check) double **{
    if ($1 == 0) {
      PyErr_SetString(PyExc_TypeError,"NULL Pointer not allowed");
      return NULL;
    }
}


%typemap(check) ghmm_cseq *{
    if ($1 == 0) {
      PyErr_SetString(PyExc_TypeError,"NULL Pointer not allowed");
      return NULL;
    }
}


%typemap(check) ghmm_dseq *{
    if ($1 == 0) {
      PyErr_SetString(PyExc_TypeError,"NULL Pointer not allowed");
      return NULL;
    }
}

