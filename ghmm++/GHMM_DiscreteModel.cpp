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
#include "ghmm++/GHMM_Alphabet.h"
#include "ghmm++/GHMM_Emission.h"

#include "ghmm/viterbi.h"
#include "ghmm/mes.h"
#include "ghmm/foba.h"
#include "ghmm/matrix.h"
#include "ghmm/reestimate.h"

#include <iostream>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_DiscreteModel::GHMM_DiscreteModel() {
  init();
}


GHMM_DiscreteModel::GHMM_DiscreteModel(model* my_model) {
  init();
  c_model = my_model;
  buildCppData();
}


GHMM_DiscreteModel::GHMM_DiscreteModel(int number_of_states, int my_M, double my_prior) {
  init(number_of_states,my_M,my_prior);
}


GHMM_DiscreteModel::GHMM_DiscreteModel(GHMM_Alphabet* my_alphabet) 
{  
  init(my_alphabet);
}


GHMM_DiscreteModel::~GHMM_DiscreteModel() 
{
  cerr << " Delete " << toString() << endl;
  /* frees model. */  
  model_free(&c_model);
  if (own_alphabet)
    SAFE_DELETE(alphabet);
  cleanCPP();
}

void GHMM_DiscreteModel::setNodeTag(const string& tag)       
{ GHMM_AbstractModel::setNodeTag( tag ); }

void GHMM_DiscreteModel::setTransitionTag(const string& tag) 
{ GHMM_AbstractModel::setTransitionTag( tag ); }

const char* GHMM_DiscreteModel::toString() const 
{
  return "GHMM_DiscreteModel";
}


GHMM_IntVector* GHMM_DiscreteModel::viterbi(GHMM_Sequences* sequences, 
					    int index, double *log_p) const {
  double my_logp;

  if (!log_p)
    log_p = &my_logp;

  int len = sequences->getLength(index);

  return new GHMM_IntVector(::viterbi(c_model,sequences->getIntSequence(index),len,log_p),len);
}


double GHMM_DiscreteModel::viterbi_logp(GHMM_Sequences* seq, int index, int* state_seq) {
  return ::viterbi_logp(c_model,seq->getIntSequence(index),seq->getLength(index),state_seq);
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


GHMM_Sequences* GHMM_DiscreteModel::generate_sequences(int seed, int global_len, long seq_number, int maxT) const {

  return new GHMM_Sequences(model_generate_sequences(c_model,seed,global_len,seq_number, maxT));
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


GHMM_DoubleMatrix* GHMM_DiscreteModel::foba_forward(GHMM_Sequences* seq, int index, 
						    GHMM_DoubleVector *scale, double *log_p) const {
  int len = seq->getLength(index);
  GHMM_DoubleMatrix *alpha = new GHMM_DoubleMatrix(len,c_model->N);

  bool delete_scale = false;
  if (! scale) {
    scale        = new GHMM_DoubleVector();
    delete_scale = true;
  }

  scale->resize(len);

  int result = ::foba_forward(c_model,seq->getIntSequence(index),len,alpha->c_matrix,scale->c_vector,log_p);

  if (result == -1)
    SAFE_DELETE(alpha);

  if (delete_scale)
    SAFE_DELETE(scale);

  return alpha;
}


int GHMM_DiscreteModel::foba_backward(GHMM_Sequences* seq, int index, double **beta, 
				      const double *scale) const {
  return ::foba_backward(c_model,seq->getIntSequence(index),seq->getLength(index),beta,scale);
}


int GHMM_DiscreteModel::foba_logp(GHMM_Sequences* seq, int index, double *log_p) const {
  return ::foba_logp(c_model,seq->getIntSequence(index),seq->getLength(index),log_p);
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


void GHMM_DiscreteModel::addState(const string& my_id) {
  int i;

  c_model->N += 1; /* increase number of states. */
  c_model->s  = (state*) realloc(c_model->s,sizeof(state) * c_model->N);

  int index = c_model->N - 1;
  /* initialize new state. */
  c_model->s[index].pi = 0;
  c_model->s[index].b  = (double*) malloc(sizeof(double) * c_model->M);
  /* output probabilities are initialized with 0. */
  for (i = 0; i < c_model->M; ++i)
    c_model->s[index].b[i] = 0;
  c_model->s[index].out_id     = NULL;
  c_model->s[index].in_id      = NULL;
  c_model->s[index].out_a      = NULL;
  c_model->s[index].in_a       = NULL;
  c_model->s[index].out_states = 0;
  c_model->s[index].in_states  = 0;
  c_model->s[index].fix        = 0;

  GHMM_State* state = new GHMM_State(this,index,&c_model->s[index]);
  state->setID(my_id);
  states.push_back(state);

  /* reallocation of c_model->c might change position of c_states
     in memory. */
  for (i = 0; i < index; ++i)
    states[i]->c_state = &c_model->s[i];
}


void GHMM_DiscreteModel::init() {
  attributes["type"] = "discrete";
  alphabet           = NULL;
  c_model            = NULL;
  own_alphabet       = false;
}

void GHMM_DiscreteModel::init(GHMM_Alphabet *my_alphabet) {
  attributes["type"] = "discrete";
  alphabet           = my_alphabet;
  c_model            = NULL;
  own_alphabet       = true; 

  c_model = (model*) calloc(1,sizeof(model));
  if (!c_model) {
    fprintf(stderr,"GHMM_DiscreteModel::GHMM_DiscreteModel() could not allocate c_model\n");
    exit(1);
  }

  c_model->N        = 0;
  c_model->M        = alphabet->size();
  buildCppData();
}

void GHMM_DiscreteModel::init(int number_of_states, int my_M, double my_prior) {
  init();

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
  c_model->s       = (state*) malloc(sizeof(state) * max(c_model->N,1));
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


void GHMM_DiscreteModel::cleanCPP() {
  unsigned int i;

  for (i = 0; i < states.size(); ++i)
    SAFE_DELETE(states[i]);
  states.clear();
}


GHMM_Alphabet* GHMM_DiscreteModel::getAlphabet() const 
{
  return alphabet;
}


void GHMM_DiscreteModel::XMLIO_finishedReading() {
  unsigned int i;

  if (c_model == NULL)
    {
      c_model = (model*) calloc(1,sizeof(model));

      if (!c_model) {
	fprintf(stderr,"GHMM_ContinuousModel::XMLIO_finishedReading() could not allocate c_model\n");
	exit(1);
      }
    }

  int alphabet_size = 0;
  if (getAlphabet())
    alphabet_size = getAlphabet()->size();
  else
    alphabet_size = states[0]->demission->weights.size();

  c_model->N       = states.size();
  c_model->M       = alphabet_size;
  c_model->prior   = -1;
  c_model->s       = (state*) malloc(sizeof(state) * c_model->N);

  /* fill states. */
  for (i = 0; i < states.size(); ++i)
    states[i]->fillState(&c_model->s[i]);

  /* check whether sum of initial probs is 1. */
  //if (check() == -1)  to allow silent states
  //  exit(1);

  /* we dont need the transitions any more. */
  cleanTransitions();
}


XMLIO_Element* GHMM_DiscreteModel::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) 
{

  /* For the simple XML format */
  if (tag == "alphabet") { 
    alphabet     = new GHMM_Alphabet();
    own_alphabet = true;
    return alphabet;
  }

  return GHMM_AbstractModel::XMLIO_startTag(tag,attrs);
}


GHMM_ModelType GHMM_DiscreteModel::getModelType() const {
  return GHMM_DISCRETE;
}
