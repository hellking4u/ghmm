/*
  created: 21 Jan 2002 by Peter Pipenbacher
  authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
  file   : $Source$
  $Id$

  __copyright__

*/

#include <xmlio/XMLIO_Definitions.h>
#include "ghmm++/GHMM_IntVector.h"
#include "ghmm++/GHMM_DoubleVector.h"
#include "ghmm++/GHMM_DiscreteModel.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_DoubleMatrix.h"
#include "ghmm/viterbi.h"
#include "ghmm/mes.h"
#include "ghmm/foba.h"
#include "ghmm/matrix.h"
#include "ghmm/reestimate.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DiscreteModel::GHMM_DiscreteModel(model* my_model) {
  c_model = my_model;

  buildCppData();
}


GHMM_DiscreteModel::GHMM_DiscreteModel(int number_of_states, int my_M, double my_prior) {
  int i;
  int j;

  c_model = (model*) calloc(1,sizeof(model));
  if (!c_model) {
    fprintf(stderr,"GHMM_DiscreteModel::GHMM_DiscreteModel() could not allocate c_model\n");
    exit(1);
  }

  c_model->N       = number_of_states;
  c_model->M       = my_M;
  c_model->prior   = my_prior;
  c_model->s       = (state*) malloc(sizeof(state) * c_model->N);
  /* initialize all states. */
  for (i = 0; i < number_of_states; ++i) {
    c_model->s[i].pi         = 0;
    c_model->s[i].b          = (double*) malloc(sizeof(double) * my_M);
    /* output probabilities are initialized with 0. */
    for (j = 0; j < my_M; ++j)
      c_model->s[i].b[j] = 0;
    c_model->s[i].out_id     = NULL;
    c_model->s[i].in_id      = NULL;
    c_model->s[i].out_a      = NULL;
    c_model->s[i].in_a       = NULL;
    c_model->s[i].out_states = 0;
    c_model->s[i].in_states  = 0;
    c_model->s[i].fix        = 0;
  }

  buildCppData();
}


GHMM_DiscreteModel::~GHMM_DiscreteModel() {
  /* frees model. */  
  model_free(&c_model);
}


const char* GHMM_DiscreteModel::toString() const {
  return "GHMM_DiscreteModel";
}


GHMM_IntVector* GHMM_DiscreteModel::viterbi(GHMM_Sequences* sequences, 
					    int index, double *log_p) const {
  int len = sequences->getLength(index);

  return new GHMM_IntVector(::viterbi(c_model,sequences->getIntSequence(index),len,log_p),len);
}


GHMM_DiscreteModel* GHMM_DiscreteModel::copy() const {
  return new GHMM_DiscreteModel(model_copy(c_model));
}


int GHMM_DiscreteModel::check() const {
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


GHMM_Sequences* GHMM_DiscreteModel::generate_sequences(int seed, int global_len, long seq_number) const {
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


void GHMM_DiscreteModel::print(FILE *file) const {
  model_print(file,c_model);
}


void GHMM_DiscreteModel::A_print(FILE *file, char *tab, char *separator, char *ending) const {
  model_A_print(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::B_print(FILE *file, char *tab, char *separator, char *ending) const {
  model_B_print(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::Pi_print(FILE *file, char *tab, char *separator, char *ending) const {
  model_Pi_print(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::fix_print(FILE *file, char *tab, char *separator, char *ending) const {
  model_fix_print(file,c_model,tab,separator,ending);
}
 

void GHMM_DiscreteModel::A_print_transp(FILE *file, char *tab, char *separator, char *ending) const {
  model_A_print_transp(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::B_print_transp(FILE *file, char *tab, char *separator, char *ending) const {
  model_B_print_transp(file,c_model,tab,separator,ending);
}


void GHMM_DiscreteModel::Pi_print_transp(FILE *file, char *tab, char *ending) {
  model_Pi_print_transp(file,c_model,tab,ending);
}


void GHMM_DiscreteModel::states_print(FILE *file) {
  model_states_print(file,c_model);
}


double GHMM_DiscreteModel::prob_distance(GHMM_DiscreteModel* m, int maxT, int symmetric, int verbose) {
  return model_prob_distance(c_model,m->c_model,maxT,symmetric,verbose);
}


GHMM_DoubleMatrix* GHMM_DiscreteModel::fobaForward(GHMM_Sequences* seq, int index, 
						   GHMM_DoubleVector *scale, double *log_p) const {
  int len = seq->getLength(index);
  GHMM_DoubleMatrix *alpha = new GHMM_DoubleMatrix(len,c_model->N);

  bool delete_scale = false;
  if (! scale) {
    scale        = new GHMM_DoubleVector();
    delete_scale = true;
  }

  scale->resize(len);

  int result = foba_forward(c_model,seq->getIntSequence(index),len,alpha->c_matrix,scale->c_vector,log_p);

  if (result == -1)
    SAFE_DELETE(alpha);

  if (delete_scale)
    SAFE_DELETE(scale);

  return alpha;
}


int GHMM_DiscreteModel::fobaBackward(GHMM_Sequences* seq, int index, double **beta, 
			      const double *scale) const {
  return foba_backward(c_model,seq->getIntSequence(index),seq->getLength(index),beta,scale);
}


int GHMM_DiscreteModel::fobaLogp(GHMM_Sequences* seq, int index, double *log_p) const {
  return foba_logp(c_model,seq->getIntSequence(index),seq->getLength(index),log_p);
}


state* GHMM_DiscreteModel::getCState(int index) const {
  if (index >= c_model->N) {
    fprintf(stderr,"GHMM_DiscreteModel::getCState(int):\n");
    fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,c_model->N);
    exit(1);
  }

  return &c_model->s[index];
}


int GHMM_DiscreteModel::getNumberOfTransitionMatrices() const {
  return 1;
}


void GHMM_DiscreteModel::buildCppData() {
  /* Create C++ wrapper for all states and fill C states with usefull 
     information. */
  int i;
  for (i = 0; i < c_model->N; ++i) {
    GHMM_State* state = new GHMM_State(this,i,&c_model->s[i]);
    states.push_back(state);
  }
}


int GHMM_DiscreteModel::reestimate_baum_welch(GHMM_Sequences* seq) {
  return ::reestimate_baum_welch(c_model,seq->c_i_sequences);
}
