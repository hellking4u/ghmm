/*
 * created: 21 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include <xmlio/XMLIO_Definitions.h>
#include "ppghmm++/GHMM_DiscreteModel.h"
#include "ppghmm++/GHMM_Sequences.h"
#include "ghmm/viterbi.h"
#include "ghmm/mes.h"
#include "ghmm/foba.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DiscreteModel::GHMM_DiscreteModel(model* my_model) {
  //  int i;

  c_model = my_model;

  /* false means: C style states are owned by the model 
                  rather than by the GHMM_State object. */
  //  for (i = 0; i < c_model->N; ++i)
  //    states.push_back(new GHMM_State(this,&c_model->s[i],false));
}


GHMM_DiscreteModel::GHMM_DiscreteModel(int number_of_states, int my_M, double my_prior) {
  //  int i;

  c_model = (model*) calloc(1,sizeof(model));
  if (!c_model) {
    fprintf(stderr,"GHMM_DiscreteModel::GHMM_DiscreteModel() could not allocate c_model\n");
    exit(1);
  }

  c_model->N       = number_of_states;
  c_model->M       = my_M;
  c_model->prior   = my_prior;
  c_model->s       = (state*) malloc(sizeof(state) * c_model->N);

  /* false means: C style states are owned by the model 
                  rather than by the GHMM_State object. */
  //  for (i = 0; i < c_model->N; ++i)
  //    states.push_back(new GHMM_State(this,&c_model->s[i],false));
}


GHMM_DiscreteModel::~GHMM_DiscreteModel() {
  //  unsigned int i;

  /* frees model. */  
  model_free(&c_model);

  /* delete all states. */
  //  for (i = 0; i < states.size(); ++i)
  //    delete states[i];
}


const char* GHMM_DiscreteModel::toString() const {
  return "GHMM_DiscreteModel";
}


int* GHMM_DiscreteModel::Viterbi(GHMM_Sequences* sequences, int index, double *log_p) const {
  return viterbi(c_model,sequences->getIntSequence(index),sequences->getLength(index),log_p);
}


GHMM_DiscreteModel* GHMM_DiscreteModel::copy() {
  return new GHMM_DiscreteModel(model_copy(c_model));
}


int GHMM_DiscreteModel::check() {
  return model_check(c_model);
}


/**
   Produces a model, which generates the given sequence with probability 1.
   The model is a strict left-right model with one state for each element 
   in the sequence and the output in state i is the i-th value in the sequence 
   with probability 1. The model also has a final state, a state with no output.
   @param seq:      sequence
   @param seq_len:  length of the sequence
   @param anz_symb: number of symbols in the sequence
*/
//GHMM_DiscreteModel::GHMM_DiscreteModel(const int *seq, int seq_len, int anz_symb) {
//  c_model = model_generate_from_sequence(seq,seq_len,anz_symb);
//}


GHMM_Sequences* GHMM_DiscreteModel::generate_sequences(int seed, int global_len, long seq_number) {
  return new GHMM_Sequences(model_generate_sequences(c_model,seed,global_len,seq_number));
}


/**
   Calculates the sum log( P( O | lambda ) ).
   Sequences, that can not be generated from the given model, are neglected.
   @return    log(P)
   @param mo model
   @param sq sequences       
*/
//double GHMM_DiscreteModel::likelihood(sequence_t *sq) {
//  return model_likelihood(c_model,sq);
//}


void GHMM_DiscreteModel::print(FILE *file) {
  model_print(file,c_model);
}


/**
   Writes transition matrix of a model.
   @param file: output file
   @param tab:  format: leading tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
void GHMM_DiscreteModel::A_print(FILE *file, char *tab, char *separator, char *ending) {
  model_A_print(file,c_model,tab,separator,ending);
}


/**
   Writes output matrix of a model.
   @param file: output file
   @param tab:  format: leading tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
void GHMM_DiscreteModel::B_print(FILE *file, char *tab, char *separator, char *ending) {
  model_B_print(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::Pi_print(FILE *file, char *tab, char *separator, char *ending) {
  model_Pi_print(file,c_model,tab,separator,ending);
}


/**
   Writes fix vector of a matrix.
   @param file: output file
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
void GHMM_DiscreteModel::fix_print(FILE *file, char *tab, char *separator, char *ending) {
  model_fix_print(file,c_model,tab,separator,ending);
}
 

void GHMM_DiscreteModel::A_print_transp(FILE *file, char *tab, char *separator, char *ending) {
  model_A_print_transp(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::B_print_transp(FILE *file, char *tab, char *separator, char *ending) {
  model_B_print_transp(file,c_model,tab,separator,ending);
}

/**
   Writes transposed initial allocation vector of a matrix.
   @param file: output file
   @param mo:   model
   @param tab:  format: leading Tabs
   @param separator: format: seperator for columns
   @param ending:    format: end of a row  
*/
//void GHMM_DiscreteModel::model_Pi_print_transp(FILE *file, model *mo, char *tab, char *ending) {
//}


/** 
    Writes the parameters of a model sorted by states. 
    Is not very concise.   
    @param file: output file
    @param mo:   model
*/
//void GHMM_DiscreteModel::model_states_print(FILE *file, model *mo) {
//}


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
//double GHMM_DiscreteModel::model_prob_distance(model *m0, model *m, int maxT, int symmetric, int verbose) {
//}


/** Generates all possible integer sequence of lenght n from an alphabet with
    M letters. Use lexicographical ordering. Memory allocation here.
    @param n      length of sequences
    @param M     size of alphabet
    @return array of generated integer sequences
*/
//sequence_t *GHMM_DiscreteModel::sequence_lexWords(int n, int M) {
//}


int GHMM_DiscreteModel::fobaForward(GHMM_Sequences* seq, int index, double **alpha, 
			     double *scale, double *log_p) {
  return foba_forward(c_model,seq->getIntSequence(index),seq->getLength(index),alpha,scale,log_p);
}


int GHMM_DiscreteModel::fobaBackward(GHMM_Sequences* seq, int index, double **beta, 
			      const double *scale) {
  return foba_backward(c_model,seq->getIntSequence(index),seq->getLength(index),beta,scale);
}


int GHMM_DiscreteModel::fobaLogp(GHMM_Sequences* seq, int index, double *log_p) {
  return foba_logp(c_model,seq->getIntSequence(index),seq->getLength(index),log_p);
}


state* GHMM_DiscreteModel::getState(int index) const {
  if (index >= c_model->N) {
    fprintf(stderr,"GHMM_DiscreteModel::getState(int):\n");
    fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,c_model->N);
    exit(1);
  }

  return &c_model->s[index];
}


int GHMM_DiscreteModel::getNumberOfTransitionMatrices() const {
  return 1;
}
