
#include <ghmm/sdmodel.h>
#include "ghmm++/GHMM_GMLState.h"
#include "ghmm++/GHMM_DoubleMatrix.h"
#include "ghmm++/GHMM_IntVector.h"
#include "ghmm++/GHMM_DoubleVector.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_DiscreteModelT.hh"



#ifdef HAVE_NAMESPACES
using namespace std;
#endif

void GHMM_DiscreteModelT::addState(const string& my_id)
{
      int i;
      
      c_model->N += 1; /* increase number of states. */
      c_model->s  = (sdstate*) realloc(c_model->s,sizeof(sdstate) * c_model->N);
      
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

      GHMM_GMLState* state = new GHMM_GMLState(this,index,&c_model->s[index]);
      state->setID(my_id);
      states.push_back(state);
      
      /* reallocation of c_model->c might change position of c_states
	 in memory. */
      for (i = 0; i < index; ++i)
	states[i]->c_sdstate = &c_model->s[i];
      
}


GHMM_DoubleMatrix* GHMM_DiscreteModelT::foba_forward(GHMM_Sequences* seq, int index, GHMM_DoubleVector* scale, 				double *log_p) const
{
  int len = seq->getLength(index);
  GHMM_DoubleMatrix *alpha = new GHMM_DoubleMatrix(len,c_model->N);
  
  bool delete_scale = false;
  if (! scale) {
    scale        = new GHMM_DoubleVector();
    delete_scale = true;
  }
  
  scale->resize(len);
  
  int result = ::sdfoba_forward(c_model,seq->getIntSequence(index),len,alpha->c_matrix,scale->c_vector,log_p);
  
  if (result == -1)
    SAFE_DELETE(alpha);
  
  if (delete_scale)
    SAFE_DELETE(scale);
  
  return alpha;
}


/**
   Viterbi algorithm. Calculates the Viterbi path (the optimal path trough
   the model) and the Viterbi probability to a given model and a given 
   sequence.
     @return Viterbi path
     @param sequences  sequences structure
     @param index      index of sequence to take
     @param log_p      probability of the sequence in the Viterbi path 
     (return value).
  */
GHMM_IntVector* GHMM_DiscreteModelT::viterbi(GHMM_Sequences* sequences, int index, double* log_p) const
{
  double my_logp;
  
  if (!log_p)
    log_p = &my_logp;
  
  int len = sequences->getLength(index);
  
  return new GHMM_IntVector(sdviterbi(this->c_model,sequences->getIntSequence(index),len,log_p),len);
}


 /** */
void GHMM_DiscreteModelT::cleanCPP()
{
  unsigned int i;
  for (i = 0; i < states.size(); ++i)
    SAFE_DELETE(states[i]);
  states.clear();
}

