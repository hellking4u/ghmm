/*
  created: 05 Feb 2002 by Peter Pipenbacher
  authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
  file   : $Source$
  $Id$

  __copyright__
  
 */

#include "ghmm/sreestimate.h"
#include "ghmm/sviterbi.h"
#include <xmlio/XMLIO_Document.h>
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_State.h"
#include "ghmm++/GHMM_Transition.h"
#include "ghmm++/GHMM_Emission.h"
#include "ghmm++/GHMM_IntVector.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_ContinuousModel::GHMM_ContinuousModel() {
  /* Do not initialize anything. It will be read from xml file. */
  c_model            = NULL;
  attributes["type"] = "continuous";
}


GHMM_ContinuousModel::GHMM_ContinuousModel(int N, int M, int cos, density_t density,
					   double prior) {
  int i;
  int j;

  attributes["type"] = "continuous";

  c_model = (smodel*) calloc(1,sizeof(smodel));
  if (!c_model) {
    fprintf(stderr,"GHMM_ContinuousModel::GHMM_ContinuousModel() could not allocate c_model.\n");
    exit(1);
  }

  c_model->N       = N;
  c_model->M       = M;
  c_model->cos     = cos;
  c_model->density = density;
  c_model->prior   = prior;
  c_model->s       = (sstate*) malloc(sizeof(sstate) * c_model->N);
  /* initialize all states. */
  for (i = 0; i < N; ++i) {
    c_model->s[i].pi         = 0;
    c_model->s[i].c          = (double*) malloc(sizeof(double) * M);
    c_model->s[i].mue        = (double*) malloc(sizeof(double) * M);
    c_model->s[i].u          = (double*) malloc(sizeof(double) * M);
    /* output probabilities are initialized with 0. */
    for (j = 0; j < M; ++j) {
      c_model->s[i].c[j]   = 0;
      c_model->s[i].mue[j] = 0;
      c_model->s[i].u[j]   = 0;
    }
    c_model->s[i].out_id     = NULL;
    c_model->s[i].in_id      = NULL;
    c_model->s[i].out_a      = (double**) malloc(sizeof(double) * cos);
    c_model->s[i].in_a       = (double**) malloc(sizeof(double) * cos);
    for (j = 0; j < cos; ++j) {
      c_model->s[i].out_a[j] = NULL;
      c_model->s[i].in_a[j]  = NULL;
    }
    c_model->s[i].out_states = 0;
    c_model->s[i].in_states  = 0;
    c_model->s[i].fix        = 0;
  }

  /* create C++ data structure. */
  buildCppData();
}


GHMM_ContinuousModel::~GHMM_ContinuousModel() {
  clean();
}


const char* GHMM_ContinuousModel::toString() const {
  return "GHMM_ContinuousModel";
}


sstate* GHMM_ContinuousModel::getCState(int index) const {
  if (index >= c_model->N) {
    fprintf(stderr,"GHMM_ContinuousModel::getCState(int):\n");
    fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,c_model->N);
    exit(1);
  }

  return &c_model->s[index];
}


GHMM_Sequences* GHMM_ContinuousModel::generate_sequences(int seed, int global_len, long seq_number, 
							 long label, int Tmax) const {
  return new GHMM_Sequences(smodel_generate_sequences(c_model,seed,global_len,seq_number,label,Tmax));
}


void GHMM_ContinuousModel::print(FILE *file) const {
  smodel_print(file,c_model);
}


int GHMM_ContinuousModel::reestimate_baum_welch(GHMM_Sequences* seq, double* logp, double eps, int max_iter) {
  smosqd_t cs;

  cs.smo      = c_model;
  cs.sqd      = seq->c_d_sequences;
  cs.logp     = logp;
  cs.eps      = eps;
  cs.max_iter = max_iter;

  return sreestimate_baum_welch(&cs);
}


XMLIO_Element* GHMM_ContinuousModel::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) 
{
  return GHMM_AbstractModel::XMLIO_startTag(tag,attrs);
}


void GHMM_ContinuousModel::XMLIO_finishedReading() {
  unsigned int i;

  c_model = (smodel*) calloc(1,sizeof(smodel));
  if (!c_model) {
    fprintf(stderr,"GHMM_ContinuousModel::XMLIO_finishedReading() could not allocate c_model\n");
    exit(1);
  }

  c_model->N       = states.size();
  c_model->M       = 1;
  c_model->cos     = 1;
  c_model->prior   = -1;
  c_model->s       = (sstate*) malloc(sizeof(sstate) * c_model->N);

  c_model->density = normal;
  if (states.size() > 0)
    c_model->density = states[0]->cemission->density;
  for (i = 1; i < states.size(); ++i)
    if (c_model->density != states[i]->cemission->density) {
      fprintf(stderr,"Not all gaussian functions are equal.\nThis is not yet supported by the library.\n");
      exit(1);
    }

  /* fill up states. */
  for (i = 0; i < states.size(); ++i)
    states[i]->fillState(&c_model->s[i]);

  /* check whether sum of initial probs is 1. */
  if (check() == -1) 
    exit(1);

  /* we dont need the transitions any more. */
  cleanTransitions();
}


int GHMM_ContinuousModel::check() const {
  return smodel_check(c_model);
}


int GHMM_ContinuousModel::getNumberOfTransitionMatrices() const {
  return c_model->cos;
}


void GHMM_ContinuousModel::read(const string& filename) {
  int smo_number;

  /* first clean up old model. */
  clean();

  smodel** models = smodel_read(filename.c_str(),&smo_number);

  if (smo_number > 0)
    copyFromModel(models[0]);

  int i;
  for (i = 0; i < smo_number; ++i)
    smodel_free(&models[i]);

  free(models);
}


void GHMM_ContinuousModel::clean() {
  /* frees model. */  
  if (c_model)
    smodel_free(&c_model);

  GHMM_AbstractModel::clean();
}


void GHMM_ContinuousModel::copyFromModel(smodel* smo) {
  clean();
  c_model = smodel_copy(smo);

  buildCppData();
}


void GHMM_ContinuousModel::buildCppData() {
  /* Create C++ wrapper for all states and fill C states with usefull 
     information. */
  int i;
  for (i = 0; i < c_model->N; ++i) {
    GHMM_State* state = new GHMM_State(this,i,&c_model->s[i]);
    states.push_back(state);
  }
}


GHMM_ModelType GHMM_ContinuousModel::getModelType() const {
  return GHMM_CONTINOUS;
}


GHMM_IntVector* GHMM_ContinuousModel::viterbi(GHMM_Sequences* seq, int index, double *log_p) {
  double my_logp;

  if (!log_p)
    log_p = &my_logp;

  int len = seq->getLength(index);

  return new GHMM_IntVector(sviterbi(c_model,seq->getDoubleSequence(index),len,log_p),len);
}
