/* author       : Wasinee Rungsarityotin and Benjamin Georgi
 *  filename    : ghmmwrapper/gql.c
 *  created      : DATE: September, 2003
 *
 * $Id$
 */

/*
  __copyright__
*/

%module ghmmwrapper

%{
#include <stdio.h>
#include <ghmm/ghmm.h>
#include <ghmm/vector.h>
#include <ghmm/sequence.h>
#include <ghmm/scanner.h>
#include <ghmm/smodel.h>
#include <ghmm/rng.h>
#include <ghmm/reestimate.h>
%}

%include carrays.i
%include cmalloc.i
%include cpointer.i
%include cstring.i
%pointer_functions(int, intp)

struct sequence_t {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq.      
   */
  int **seq;
  /** array of sequence length */
  int *seq_len;
  /**  array of sequence labels */
  long *seq_label;
  /**  array of sequence IDs*/
  double *seq_id;
  /** positiv! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;
};
typedef struct sequence_t sequence_t;

/** @name struct sequence_d_t
    Sequence structure for double sequences. 
    Contains an array of sequences and corresponding
    data like sequnce label, sequence weight, etc. Sequences may have different
    length.    
 */
struct sequence_d_t {
  /** sequence array. sequence[i] [j] = j-th symbol of i-th seq. */
  double **seq;
  /** array of sequence length */
  int *seq_len;
  /**  array of sequence labels */
  long *seq_label;
  /**  array of sequence IDs*/
  double *seq_id;
  /** positive! sequence weights.  default is 1 = no weight */
  double *seq_w;
  /** total number of sequences */
  long seq_number;
  /** sum of sequence weights */
  double total_w;  
};
typedef struct sequence_d_t sequence_d_t;


/** @name state
    The basic structure, keeps all parameters that belong to a state. 
*/
struct state {
  /** Initial probability */ 
  double pi;
  /** Output probability */
  double *b;
  /** ID of the following state */ 
  int *out_id;  
  /** ID of the previous state */    
  int *in_id;
  /** Transition probability to a successor */
  double *out_a; 
  /** Transition probablity to a precursor */
  double *in_a;
  /** Number of successor states */     
  int out_states; 
  /** Number of precursor states */
  int in_states;  
  /** if fix == 1 --> b stays fix during the training */
  int fix;
};
typedef struct state state;

/** @name model
    The complete HMM. Contains all parameters, that define a HMM.
*/
struct model {
  /** Number of states */
  int N;
  /** Number of outputs */   
  int M;   
  /** Vector of the states */
  state *s; 
  /** Prior for the a priori probability for the model.
      A value of -1 indicates that no prior is defined. */
  double prior;
};
typedef struct model model;


/** @name model_direct
    The complete HMM. Keeps the model parameters in a matrix form. 
*/
struct model_direct {
  /** Number of states */
  int N;
  /** Number of outputs */  
  int M;
  /** Prior for the a priori probability for the model.
      Gets the value -1 if no prior defined. */
  double prior;
  /** Transition matrix  */
  double **A;
  /** Output matrix */
  double **B;
  /** Initial matrix */
  double* Pi;
  /** A vector to know the states where the output should not be trained.
      Default value is 0 for all states. */
  int *fix_state;
};

typedef struct model_direct model_direct;




/**
   Memory allocation for an integer sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
extern sequence_t *sequence_calloc(long seq_number);


/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
void sequence_clean(sequence_t *sq);

/**
  Prints one array of integer sequences in a file.
  @param file       output file
  @param sequence    array of sequences
  */
void sequence_print(FILE *file, sequence_t *sequence);


/**
  copy one integer sequence. Memory for target has to be allocated outside.
  @param target  target sequence
  @param source source sequence
  @param len     length of source sequence
  */
void sequence_copy(int *target, int *source, int len);

/**
   Reads one or several arrays of double sequences. 
   Calls sequence\_read\_alloc, where reading
   and memory allocation is done. 
   @return pointer to sequence array
   @param filename    input filename
*/
sequence_d_t **sequence_d_read(const char *filename, int *sqd_number);

/**
  Frees all memory in a given array of integer sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
int sequence_free(sequence_t **sq);

/**
  Frees all memory in a given array of double sequences.
  @param sq sequence  structure
  @return 0 for succes, -1 for error
  */
int sequence_d_free(sequence_d_t **sq);

/*** !!!!!!!! TO DO: Free functions for all types of pointers !!!!!!!!***/

/* Some array helpers */
%inline %{
/* Create any sort of [size] array */
  int *int_array(int size) {
    return (int *) malloc(size*sizeof(int));
  }

  /* Create any sort of [size] array */
  double *double_array(int size) {
    return (double *) malloc(size*sizeof(double));
  }
  
  void set_arrayd(double *ary, int index, double value) {
    ary[index] = value;
  }

  double get_arrayd(double *ary, int index) { return ary[index]; }

  void set_arrayint(int  *ary, int index, int value) {
    ary[index] = value;
  }

  int  get_arrayint(int  *ary, int index) { return ary[index]; }

  /* Get two dimensional array entry */
  int  get_2d_arrayint (int **ary, int index1, int index2) {return ary[index1][index2]; }
		  
  /* Get two dimensional array entry */
  double get_2d_arrayd (double **ary, int index1, int index2) {return ary[index1][index2]; }

  /* Set two dimensional array entry */
  void  set_2d_arrayd (double **ary, int index1, int index2, double val) 
    {ary[index1][index2] = val; }
		  

  /* Create a two-dimension array [size][10] */
  int (*int_array_10(int size))[10] {
     return (int (*)[10]) malloc(size*10*sizeof(int));
  }

  state *get_stateptr(state *ary, int index) { return ary + index; }

  state *arraystate(int size) {
    return (state *) malloc(size*sizeof(state));
  }

  /* return a C-pointer to an integer sequence */
  int *get_onesequence(sequence_t *seqpt, int seqnumber) { 
    return (int *) seqpt->seq[seqnumber];
  }

  sequence_d_t *get_onesequence_d(sequence_d_t **seqpt, int seqnumber) { 
    return seqpt[seqnumber];
  }

  void freearray(void *pt)  { free(pt); }

  void free_arrayi(int *pt) { free(pt); }

  void free_arrayd(double *pt) { free(pt); }

  void free_2darrayd(double **pt, int rows) 
    {
      matrix_d_free(&pt, rows);
    }

  void free_2darrayint(int **pt, int rows,int cols) 
    {
      matrix_i_free(&pt, rows);
    }
%}

/** Frees the memory of a model.
    @return 0 for succes; -1 for error
    @param mo:  pointer to a model 
*/
extern int     model_free(model **mo);

/**
   Reads in ASCII data to initialize an array of models. Memory allocation for
   the models is done here.
   @return array of pointers to the models
   @param filename:   the ASCII input file
   @param mo_number:  filled with number of models read */
extern model** model_read(char *filename, int *mo_number);

/**
   Reads in a model, where the model parameters are explicit given in
   matrix form. Memory allocation for the model is also done here.
   @return pointer to the model
   @param s:       scanner
   @param multip:  multiplicity; gives how many copies should 
   be made of the model */
extern model*  model_direct_read(scanner_t *s, int *multip);

/**
   Produces simple left-right models given sequences. 
   The function "model_generate_from_sequence" is called for each 
   model that should be made. The sequences are read in from the
   ASCII file and thrown away again when leaving the function.
   @return vector of models
   @param s:          scanner
   @param new_models: number of models to produce */
extern model **model_from_sequence_ascii(scanner_t *s, long *mo_number);

/** 
    Produces simple left-right models given sequences. The sequences
    are not read in from file, but exists already as a structur.
    @return vector of models
    @param s:          scanner
    @param new_models: number of models to produce */
extern model **model_from_sequence(sequence_t *sq, long *mo_number);

/**
   Copies a given model. Allocates the necessary memory.
   @return copy of the model
   @param mo:  model to copy */
extern model*  model_copy(const model *mo);

/**
   Tests if all standardization requirements of model are fulfilled. 
   (That is, if the sum of the probabilities is 1).
   @return 0 for succes; -1 for error
   @param mo:  model to test */
extern int     model_check(const model* mo);

/**
   Tests if number of states and number of outputs in the models match.
   @return 0 for succes; -1 for error
   @param mo:           vector of models
   @param model_number: numbr of models */
extern int     model_check_compatibility(model **mo, int model_number);

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
extern model*  model_generate_from_sequence(const int *seq, int seq_len, 
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
extern sequence_t *model_generate_sequences(model* mo, int seed, int global_len,
				     long seq_number, int Tmax);


/* Important! initialise rng  */
extern void gsl_rng_init(void);

/* Reestimate Baum-Welch */
extern int reestimate_baum_welch(model *mo, sequence_t *sq);

/******  Viterbi ********/
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
extern int *viterbi(model *mo, int *o, int len, double *log_p);

/*
// Some array helpers
%inline %{
  viterbi_t *call_viterbi(model *mo, int *o, int len) {
    double log_p;
    int   *viterbi_path;
    viterbi_t *answer;
    viterbi_path = viterbi(model *mo, int *o, int len, &log_p);
    answer->log_p = log_p;
    answer->viterbi_path = viterbi_path;
    return answer;
  }
  %}*/

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
extern double viterbi_logp(model *mo, int *o, int len, int *state_seq);


/******* Forward , backward ******/

/** Forward-Algorithm.
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param length: length of sequence
  @param alpha:  alpha[t][i]
  @param scale:  scale factors
  @param log\_p:  log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int foba_forward(model *mo, const int *O, int length, double **alpha, 
		 double *scale, double *log_p);

/** 
  Backward-Algorithm. 
  Calculates beta[t][i] given an integer sequence and a model. Scale factors 
  given as parameter (come from sfoba\_forward).
  @param  mo      model
  @param O          sequence
  @param length   length of sequence
  @param beta     beta[t][i]
  @param scale    scale factors
  @return 0 for success, -1 for error
  */
extern int foba_backward(model *mo, const int *O, int length, double **beta, 
		  const double *scale);

/**
  Calculation of  log( P(O|lambda) ). 
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed.
  @param  mo      model
  @param O        sequence
  @param len       length of sequence
  @param log\_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int foba_logp(model *mo, const int *O, int len, double *log_p);


/******* Utility function: Matrix allocation and destruction ******/
/**
  Allocation of a double matrix. 
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
extern double** matrix_d_alloc(int rows, int columns);

/**
  Copying and allocation of a double matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  @param copymatrix: matrix to copy 
  */
extern double** matrix_d_alloc_copy(int rows, int columns, double **copymatrix);

/**
  Free the memory of a double matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
extern int matrix_d_free(double ***matrix, long rows);

/**
  Allocation of a integer matrix.
  @return pointer to a matrix
  @param rows: number of rows
  @param columns: number of columns
  */
extern int** matrix_i_alloc(int rows, int columns);

/**
  Free the memory of a integer matrix.
  @return 0 for succes; -1 for error
  @param  matrix: matrix to free
  @param  rows: number of rows
  */
extern int matrix_i_free(int ***matrix, long rows); 

/**
  Writes a double matrix (without parenthesis).
  @param file:       output file
  @param matrix:     matrix to write
  @param rows:       number of rows
  @param columns:    number of columns
  @param tab:        format: leading tabs
  @param separator:  format: separator for columns
  @param ending:     format: end of a row  
  */
extern void matrix_d_print(FILE *file, double **matrix, int rows, int columns, 
		    char *tab, char *separator, char *ending);




/*============================================================================*/
/** Wrappers for Continuous model
 *  Structure: sstate, smodel
 */

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


/** Free memory smodel 
    @return 0: success, -1: error
    @param smo  pointer pointer of smodel */
extern int     smodel_free(smodel **smo);

/** Reads an ascii file with specifications for one or more smodels.
    All parameters in matrix or vector form.
    This is necessary whenever an initial model is needed (e.g. 
    training) or sequences have to be generated from trained models.
    For each smodel block smodel\_read\_block() is called.
   @return vector of read smodels
   @param filename   input ascii file
   @param smo_number  number of read smodels */
extern smodel** smodel_read(const char *filename, int *smo_number);

%inline %{

  int call_sequence_free(sequence_t *sq) {
    return sequence_free(&sq);
  }

  int call_sequenced_free(sequence_d_t *sq) {
    return sequence_d_free(&sq);
  }

  smodel* get_onesmodel(smodel** smos, int snum) {
    return smos[snum];
  }
  void call_smodel_free(smodel *smo) {
    smodel_free(&smo);
  }
  void smodel_print_stdout(smodel *smo) {
    smodel_print(stdout, smo);
  }
  sstate *get_sstate(smodel *smo, int k) {
    return &(smo->s[k]);
  }
%}

/**
   Copies one smodel. Memory alloc is here.
   @return pointer to smodel copy
   @param smo   smodel to be copied  */
extern smodel*  smodel_copy(const smodel *smo);

/**
   Produces sequences to a given model. All memory that is needed for the 
    sequences is allocated inside the function. It is possible to define
    the length of the sequences global (global_len > 0) or it can be set 
    inside the function, when a final state in the model is reach (a state
    with no output). If the model has no final state, the sequences will
    have length MAX_SEQ_LEN.
*/
extern sequence_d_t *smodel_generate_sequences(smodel* smo, int seed, int global_len,
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

extern smodel *smodel_alloc_fill(int N, int M, int cos, double prior, int density);

extern void smodel_set_pivector(smodel *smo, int i, double piv);

extern void smodel_set_fixvector(smodel *smo, int i, double fixv);

extern void smodel_set_transition(smodel *smo, int i, int j, int cos, double prob);

extern double smodel_get_transition(smodel *smo, int i, int j, int cos);

extern void smodel_set_mean(smodel *smo, int i, double *mu);

extern void smodel_set_variance(smodel *smo, int i, double *variance);

extern void call_smodel_print(char *filename, smodel *smo);

extern int smodel_likelihood(smodel *smo, sequence_d_t *sqd, double *log_p);

extern int smodel_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps);

extern int smodel_sorted_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps, int *seq_rank);
