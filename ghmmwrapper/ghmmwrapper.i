/* author       : Wasinee Rungsarityotin and Benjamin Georgi
 *  filename    : ghmmwrapper/ghmmwrapper.i
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
#include <ghmm/foba.h>
#include <ghmm/scluster.h>
#include "sdclass_change.h"
#include "read_cxml.h"
%}

%include carrays.i
%include cmalloc.i
%include cpointer.i
%include cstring.i
%pointer_functions(int, intp)


/*=============================================================================================
  =============================== Random Number Generator (RNG ================================= */

/* The global RNG */
extern gsl_rng * RNG; 

/* Important! initialise rng  */
extern void gsl_rng_init(void);

/* Initialise random timeseed */
extern void gsl_rng_timeseed(gsl_rng * r);

%inline %{
	void time_seed(){
		gsl_rng_timeseed(RNG);
	}
%}
		
/*=============================================================================================
  ================= Utility functions: Matrix allocation and destruction  =====================*/

/*
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
extern int matrix_d_free(double ***matrix,int row);

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



				
/*=============================================================================================
  =============================== sequence.c  ============================================== */
		
		
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

%pointer_functions(sequence_t **, sequence_setPtr)

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

%pointer_functions(sequence_d_t **, sequence_d_setPtr)

/**
   Memory allocation for an integer sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
extern sequence_t *sequence_calloc(long seq_number);


/**
   Memory allocation for a double  sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
extern sequence_d_t *sequence_d_calloc(long seq_number);

/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
extern void sequence_clean(sequence_t *sq);

/**
   Cleans integer sequence pointers in sequence struct. sets 
   seq\_number to zero.
   Differs from sequence\_free since memory is not freed here. 
   @param sq sequence structure
  */
extern void sequence_d_clean(sequence_d_t *sq);

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
  Adds all integer sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
int sequence_add(sequence_t *target, sequence_t *source);

/**
  Adds all double sequences, sequence lengths etc 
  from source to target. Memory allocation is done here.
  @param target target sequence structure
  @param source  source sequence structure
  @return -1 for error, 0 for success
  */
extern int sequence_d_add(sequence_d_t *target, sequence_d_t *source);

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

extern void sequence_d_print(FILE *file, sequence_d_t *sqd, int discrete);

/*** !!!!!!!! TO DO: Free functions for all types of pointers !!!!!!!!***/

%inline%{

  /* return a C-pointer to an integer sequence */
  int *get_onesequence(sequence_t *seqpt, int seqnumber) { 
    return (int *) seqpt->seq[seqnumber];
  }

  
  sequence_t *seq_read(char* filename ){
	  int i;
	  sequence_t** s;
	  //s = (sequence_d_t **) malloc(1*sizeof(sequence_d_t*));
	  s = sequence_read(filename, &i);
	  return s[0];
  }
  
  sequence_t* get_seq_ptr(sequence_t** array, int index){
	  return (sequence_t*) array[index];
  }
  
  /*** create and manipulate an array of pointers of pointers to sequence_d_t structs ***/
  sequence_d_t **sequence_d_t_array(int size) {
     int i;
	 sequence_d_t **s;
	 s = (sequence_d_t **) malloc(size*sizeof(sequence_d_t*));
	 for(i =0;i<size;i++){	
		 s[i] = NULL;
	 }
	 return s;	 
  }

  sequence_d_t* get_seq_d_ptr(sequence_d_t** array, int index){
	  return (sequence_d_t*) array[index];
  }	  
  
  void set_seq_d_array(sequence_d_t** array, int index,sequence_d_t* seq){
	array[index] = seq;
  }	
  
  void set_sequence_d_label(sequence_d_t* seq, int seq_num, long label){
	  seq->seq_label[seq_num] = label;
  }
  
  long get_sequence_d_label(sequence_d_t* seq, int seq_num ){
	  return seq->seq_label[seq_num];
  }	  
  
 /* to free a sequence* from python use free_sequence_d */
   void free_sequence_d(sequence_d_t *seq) { 
	  sequence_d_free(&seq);}
  
  sequence_d_t *seq_d_read(char* filename ){
	  int i;
	  sequence_d_t** s;
	  //s = (sequence_d_t **) malloc(1*sizeof(sequence_d_t*));
	  s = sequence_d_read(filename, &i);
	  return s[0];
  }

  void call_sequence_print(char* ch, sequence_t* seq){
    FILE* file_name;
    file_name = fopen(ch,"at");
    sequence_print(file_name,seq);
    fclose(file_name);
  }	  

  void call_sequence_d_print(char* ch,sequence_d_t* seq, int disc){
    FILE* file_name;
    file_name = fopen(ch,"at");
    sequence_d_print(file_name,seq, disc);
    fclose(file_name);
  }	  

   
%}

/*=============================================================================================
  =============================== model.c  ============================================== */


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

  /** Flag variables for each state indicating whether the states emissions
      are tied to another state. Groups of tied states are represented
      by their tie group leader (the lowest numbered member of the group).
      
      tied_to[s] == kUntied  : s is not a tied state
      
      tied_to[s] == s        : s is a tie group leader

      tied_to[t] == s        : t is tied to state s

      Note: tied_to != NULL iff (model_type & kTiedEmissions) == 1  */
  int* tied_to; 
  
  /** Flag variables for each state giving the order of the emissions
      Classic HMMS have emission order 0, that is the emission probability
      is conditioned only on the state emitting the symbol.

      For higher order emissions, the emission are conditioned on the state s
      as well as the previous emission_order[s] observed symbols.

      Note: emission_order != NULL iff (model_type & kHigherOrderEmissions) == 1  */
  int* emission_order; 
  int  topo_order_length; /*WR*/
  int* topo_order;        /*WR*/
  	
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
   Writes a model in matrix format.
   @param file: output file
   @param mo:   model
*/
void model_print(FILE *file, model *mo); 

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

/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
extern double model_likelihood(model *mo, sequence_t *sq);


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
double model_prob_distance(model *m0, model *m, int maxT, int symmetric, int verbose);


/******** Reestimate Baum-Welch (reestimate.c) *******/
extern int reestimate_baum_welch(model *mo, sequence_t *sq);
// XXX need source code XXX 
/// extern int reestimate_baum_welch_nstep(model *mo, sequence_t *sq, int max_step, double likelihood_delta);

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
extern int *viterbi(model *mo, int *o, int len, double *log_p);

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


/******* Forward , backward (foba.c) ******/

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
  given as parameter (come from foba\_forward).
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
  Done by calling foba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed.
  @param  mo      model
  @param O        sequence
  @param len       length of sequence
  @param log\_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int foba_logp(model *mo, const int *O, int len, double *log_p);


%inline%{
	
  /* allocation of an empty model struct */
  model *new_model() {
     return (struct model *)(struct model *) calloc(1, sizeof(struct model));    
  }
    
    
  /* allocation of an array of state structs*/
  state *arraystate(int size) {
    return (state *) malloc(size*sizeof(state));
  }
  
  /* extract pointer to a state */
  state *get_stateptr(state *ary, int index) { return ary + index; }
  

  void call_model_free(model *mo) {
    model_free(&mo);
  }

  void call_model_print(char *filename, model *mo) {
    FILE *fp=fopen(filename, "a");
    if (fp == NULL) {
      fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);    
    } 
    else {
      model_print(fp, mo);
      fclose(fp);
    }
  }
  

  void call_model_free(model * mo)	{ model_free(&mo); }
  
  model *get_model_ptr(model **mo, int index) { return mo[index]; }
    
  model **cast_model_ptr(model *mo){
    model** result = &mo;
    return result;
  }   

  
  
%}



/*=============================================================================================
  =============================== sdmodel.c  ============================================== */

/** @name state
    The basic structure, keeps all parameters that belong to a sdstate. 
*/
struct sdstate {
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
typedef struct sdstate sdstate;

/** @name model
    The complete HMM. Contains all parameters, that define a HMM.
*/
struct sdmodel {
  /** Number of states */
  int N;
  /** Number of outputs */   
  int M;   
 /** smodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
  int cos;
  /** Vector of the states */
  sdstate *s; 
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
typedef struct sdmodel sdmodel;

/** Frees the memory of a model.
    @return 0 for succes; -1 for error
    @param mo:  pointer to a model */
int     sdmodel_free(sdmodel **mo);

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
sequence_t *sdmodel_generate_sequences(sdmodel* mo, int seed, int global_len,
				     long seq_number, int Tmax);


/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
extern double sdmodel_likelihood(sdmodel *mo, sequence_t *sq);

/** 
    Computes log likelihood for all sequence of
    seq_w * log( P ( O|lambda ) ). If a sequence can't be generated by smo
    error cost of seq_w * PRENALTY_LOGP are imposed.
   @return       n: number of evaluated sequences, -1: error
   @param smo   smodel
   @param sqd    sequence struct
   @param log\_p array of evaluated likelihoods
*/
extern int smodel_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps);

/******* Viterbi for switching discrete model (sdviterbi.c) *******/
int *sdviterbi(sdmodel *mo, int *o, int len, double *log_p);

/******* Forward-Algorithm. (sdfoba.c) *******
  Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
  a given double sequence and a given model.
  @param smo      model
  @param O        sequence
  @param length: length of sequence
  @param alpha:  alpha[t][i]
  @param scale:   a reference for double type, scale factors
  @param log\_p:  a reference for double type, log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int sdfoba_forward(sdmodel *mo, const int *O, int len, double **alpha, 
		 double *scale, double *log_p); 

// default class change function (dummy)
extern int cp_class_change(int *seq, int len);

extern void setSwitchingFunction( sdmodel *smd );		  

%inline%{
  
  void call_sdmodel_free(sdmodel * sdm)	{ sdmodel_free(&sdm); }
	
  /* allocation of an array of sdstate structs*/	
  sdstate *arraysdstate(int size) {
    return (sdstate *) malloc(size*sizeof(sdstate));
  }	
  
  /* extract pointer to a sdstate  */
  sdstate *get_sdstateptr(sdstate *ary, int index) { return ary + index; }
  
%}


/*=============================================================================================
  =============================== smodel.c  ============================================== */

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

/** Computes probabilistic distance of two models
    @return the distance
    @param cm0  smodel used for generating random output
    @param cm   smodel to compare with
    @param maxT  maximum output length (for HMMs with absorbing states multiple
                 sequences with a toal length of at least maxT will be 
		 generated)
    @param symmetric  flag, whether to symmetrize distance (not implemented yet)
    @param verbose  flag, whether to monitor distance in 40 steps. 
                    Prints to stdout (yuk!)
*/
extern double smodel_prob_distance(smodel *cm0, smodel *cm, int maxT, int symmetric, int verbose);


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
extern int sfoba_forward(smodel *smo, const double *O, int T, double ***b, 
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
extern int sfoba_backward(smodel *smo, const double *O, int T, double ***b,
		   double **beta, const double *scale);

/**
  Calculation of  log( P(O|lambda) ).  (sfoba.c)s
  Done by calling sfoba\_forward. Use this function if only the
  log likelihood and not alpha[t][i] is needed, alpha matrix is allocated with
  stat_matrix_d_alloc
  @param smo      model
  @param O        sequence
  @param T         length of sequence
  @param log_p    log likelihood log( P(O|lambda) )
  @return 0 for success, -1 for error
  */
extern int sfoba_logp(smodel *smo, const double *O, int T, double *log_p);





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

// extern void call_smodel_print(char *filename, smodel *smo);

extern int smodel_likelihood(smodel *smo, sequence_d_t *sqd, double *log_p);

extern int smodel_sorted_individual_likelihoods(smodel *smo, sequence_d_t *sqd, double *log_ps, int *seq_rank);

/** Viterbi for smodel (sviterbi.c)

   Viterbi algorithm: calculation of the viterbi path (best possible
   state sequenz for a given sequenz and a given model (smo)). Also 
   calculates logp according to this path, the matrices in the local_store
   struct are allocated using stat_matrix_d_alloc.
  @return        Viterbi-path 
  @param smo    model
  @param o       double-sequence
  @param T       sequence length
  @param log_p   log(p) of the sequence using the vitberbi path
  */
extern int *sviterbi(smodel *smo, double *o, int T, double *log_p);



%inline %{

  /* allocation of an empty smodel struct */
  smodel *new_smodel() {
     return (struct smodel *)(struct smodel *) calloc(1, sizeof(struct smodel));    
  }

    
  /* array of sstate structs */
  sstate *arraysstate(int size) { 
    return (sstate *) malloc(size*sizeof(sstate));
  }	
  
   sstate *get_sstate_ptr(sstate *states, int k) {
    return &(states[k]);
  }
		  	
  void call_smodel_free(smodel *smo ) {smodel_free(&smo);}

  void free_smodel_array(smodel **smo) { if (smo){ free(smo);} }
  		
  void smodel_print_stdout(smodel *smo) {
    smodel_print(stdout, smo);
  }
 
  /* extract pointer to sstate  */	
  sstate *get_sstate(smodel *smo, int k) {
    return &(smo->s[k]);
  }


  /* creation and assessor functions for an array of smodels */
  smodel **smodel_array(int size){
  	return (smodel**) malloc(size * sizeof(smodel*));
  } 

  smodel *get_smodel_ptr(smodel **smo, int index) { return smo[index]; }
 
  void set_smodel_ptr(smodel **smo_array ,smodel *smo, int index) { smo_array[index] = smo; }
  
  smodel **cast_smodel_ptr(smodel *smo){
     smodel** res = (smodel**) malloc(sizeof(smodel*));
     res[0] = smo;
     return res;
  }   
  
  // write a smodel to a file
  void call_smodel_print(char *filename, smodel *smo) {
  FILE *fp=fopen(filename, "a");
  if (fp == NULL) {
    fprintf(stderr, "call_smodel_print(0): cannot open file %s\n", filename);    
  } else {
    smodel_print(fp, smo);
    fclose(fp);
  }
}

%}



/* =============================================================================================
   ============================== scluster.c  ================================================== */

/**
   Cluster structure: All models and sequences. */
struct scluster_t{
  /** 
  Vector of SHMMs pointers */
  smodel **smo;
  /** 
    Vector of sequence_d_t pointers; to store the sequences, that  belong to the models */
  sequence_d_t **smo_seq;
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
      fucntion in case of a MAW classificator. Is calculated using smap_bayes */
  double *smo_Z_MAW;
};
/**
 */
typedef struct scluster_t scluster_t;


/**
   Frees the memory allocated for a scluster_t struct.
   @return 0 for success; -1 for error
   @param scl pointer to scl struct
*/
extern int scluster_t_free(scluster_t *scl);


/**
   Creates random labels for a vector of sequences
   @return 0 for success; -1 for error
   @param sqd vector of sequences
   @param smo_number number of models (needed to determine the interval
   for the random numbers)
*/
extern int scluster_random_labels(sequence_d_t *sqd, int smo_number);

/**
   Calculates the logarithmic probability of sequences for a model.
   @param cs sequences and model 
 */
extern void scluster_prob(smosqd_t *cs);

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
extern int scluster_best_model(scluster_t *cl, long seq_id, double **all_log_p,
			double *log_p);

/**
   Updates the cluster with additional sequences.
   @return 0 for success; -1 for error
   @param cl cluster to update
   @param sqd sequences to update the cluster with
 */
extern int scluster_update(scluster_t *cl, sequence_d_t *sqd);

/**
   Prints the input vector for scluster_hmm
 */
void scluster_print_likelihood(FILE *outfile, scluster_t *cl);
	

/**
   Writes out the final models.
   @return 0 for success; -1 for error
   @param cl cluster of models to write
   @param sqd
   @param outfile output file
   @param out_filename name of the output file
 */
int scluster_out(scluster_t *cl, sequence_d_t *sqd, FILE *outfile, char *argv[]);

/** Calculates the aposteriori prob. $\log(p(\lambda_best | O[seq\_id]))$, 
    where $\lambda_best$ is the model with highest apost. prob.
    @return 0 for success; -1 for error
    @param cl cluster
    @param sqd the sequence in question
    @param seq_id the ID of the sequence
    @param lob_apo the results
*/
extern int  scluster_log_aposteriori(scluster_t *cl, sequence_d_t *sqd, int seq_id, double *log_apo);

/**
   Avoids empty models going out as outputs by assigning them a random 
   sequence. This may lead to a produce of a new empty model - therefore
   change out sequences until a non empty model is found. (Quit after 100 
   iterations to avoid an infinite loop). 
   @return 0 for success; -1 for error
   @param sqd sequences to generate the model from
   @param cl cluster for the models
 */
extern int scluster_avoid_empty_smodel(sequence_d_t *sqd, scluster_t *cl);

/**
   Makes a cluster and reestimates the HMMs.
   @return 0 for success; -1 for error
   @param argv vector of input files, one with sequences, one with models, 
   one for output and one with labels for the sequences - in this order.
 */
extern int scluster_hmm(char *argv[]);


%inline %{
  void scluster_printl_stdout(scluster_t *scl) {
    scluster_print_likelihood(stdout, scl);
  }
  
%}


/*=============================================================================================
  =============================== sreestimate.c  ============================================== */

/** @name struct smosqd_t
    structure that combines a continuous model (smo) and an integer
    sequence struct. Is used by sreestimate\_baum\_welch for 
    parameter reestimation.
 */
struct smosqd_t {
  /** pointer of continuous model*/
  smodel *smo;
  /** sequence\_d\__t pointer */
  sequence_d_t *sqd;
  /** calculated log likelihood */
  double* logp;
  /** leave reestimation loop if diff. between successive logp values 
      is smaller than eps */
  double eps;
  /** max. no of iterations */
  int max_iter;
}; 
typedef struct smosqd_t smosqd_t;


/**
  Baum-Welch Algorithm for SHMMs.
  Training of model parameter with multiple double sequences (incl. scaling).
  New parameters set directly in hmm (no storage of previous values!).
  @return            0/-1 success/error
  @param cs         initial model and train sequences
  */
extern int sreestimate_baum_welch(smosqd_t *cs);


%inline%{
  
  /************ create and manipulate an array of pointers to smosqd_t structs **************/
  smosqd_t *smosqd_t_array(int size) {
     return (smosqd_t *) malloc(size*sizeof(smosqd_t));
  }

  void set_smosq_t_smo(smosqd_t *cs, smodel *smo, int index){
	  cs[index].smo = smo;
  }	  
  
  smosqd_t *get_smosqd_t_ptr(smosqd_t *cs, int i){ return &cs[i];}
  
  void free_smosqd_t(smosqd_t *s){
	  if(s){
		free(s);
	 }	
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

  void free_arrayi(int *pt) { free(pt); }

  /************ Create and access double[size] arrays ************/
 
   double *double_array(int size) {
     return (double *) malloc(size*sizeof(double));
  }

  void set_arrayd(double *ary, int index, double value) {
    ary[index] = value;
  }

  double get_arrayd(double *ary, int index) { return ary[index]; }
  
  void free_arrayd(double *pt) { free(pt); }
  
  /************  Create and access sort of long[size] arrays **************/
    
  long *long_array(int size) {
     return (long *) malloc(size*sizeof(long));
  }

  void set_arrayl(long *ary, int index, long value) {
    ary[index] = value;
  }
  
  double get_arrayl(long *ary, int index) { return ary[index]; }
  
  void free_arrayl(long *ary) {free(ary);}
  
  /*********** Create and access char** arrays ***********/
  char **char_array(int size) {
     return (char **) malloc(size*sizeof(char*));
  }

  void set_arraychar(char** ary, int index, char* value) {
    ary[index] = value;
  }

  char *get_arraychar(char** ary, int index) { return ary[index]; }

   
  /************  Create and access double[size1][size2] arrays ************/
  
 
  double **double_2d_array(int rows, int cols) {
    //int i;
    //double **array = (double **) malloc(rows*sizeof(double*));
    //for(i=0; i < rows; i++)
    //  array[i] = double_array(cols);
    return matrix_d_alloc(rows,cols);
  }
  
  double **double_2d_array_nocols(int rows){
	return (double **) malloc(rows*sizeof(double*));
  }	  
  
  void set_2d_arrayd_col(double **ary, int index, double *col){
	  ary[index] = col;
  }	  
  
  void set_2d_arrayd(double **ary, int index1,int index2, double value) {
    ary[index1][index2] = value;
  }
  
  double get_2d_arrayd(double **ary, int index1, int index2) { return ary[index1][index2]; }

  int *get_col_pointer_int(int **ary, int index) {
	  return ary[index];
  }
  
  double *get_col_pointer_d(double **ary, int index) {
	  return ary[index];
  }
  
  void double_2d_print(double **ary, int row, int col){
	int i,j;
	printf("(%dx%d) Matrix", row, col);
	for(i =0;i<row;i++){		  
      printf("\n");
      for(j =0;j<col;j++){		  		 
        printf("%f ",ary[i][j]);
	  }
	}
	printf("\n");  	
  }

  double** cast_ptr_d(double* array){
     double ** res = (double **) malloc(sizeof(double*));
	 res[0] = array;
	 return res;
	  
	//return (double**) array;
  }	  
  
  void free_2darrayd(double **pt,int row) { matrix_d_free(&pt,row); }
   
  /************  Create and access int[size1][size2] arrays ************/
  
   int **int_2d_array_nocols(int rows){
	return (int **) malloc(rows*sizeof(int*));
  }	  
  
  void set_2d_arrayint_col(int **ary, int index, int *col){
	  ary[index] = col;
  }	  
  
  void set_2d_arrayint(int **ary, int index1,int index2, int value) {
    ary[index1][index2] = value;
  }
  
  /* Get two dimensional array entry */
  int  get_2d_arrayint (int **ary, int index1, int index2) {return ary[index1][index2]; }
   
   int** cast_ptr_int(int* array){

	 int ** res = (int **) malloc(sizeof(int*));
	 res[0] = array;
	 return res; 
	 
	 //return (int**) array;
  }	 
  
  void freearray(void *pt)  { free(pt); }

  void free_2darrayint(int **pt, int rows,int cols) {  matrix_i_free(&pt, rows); }

 
%}


/******************* typemap ********************************/

%typemap(memberin) sequence_d_t * {
  sequence_d_t *src_pt;
  sequence_d_t *dest_pt;

  dest_pt = (sequence_d_t *) malloc(sizeof(sequence_d_t));
  src_pt  = $input;

  dest_pt.seq        = src_pt.seq;
  dest_pt.seq_len    = src_pt.seq_len;
  dest_pt.seq_label  = src_pt.seq_label;
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


%typemap(check) sequence_d_t *{
    if ($1 == 0) {
      PyErr_SetString(PyExc_TypeError,"NULL Pointer not allowed");
      return NULL;
    }
}


%typemap(check) sequence_t *{
    if ($1 == 0) {
      PyErr_SetString(PyExc_TypeError,"NULL Pointer not allowed");
      return NULL;
    }
}

