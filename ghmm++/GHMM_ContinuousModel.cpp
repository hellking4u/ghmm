/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ghmm/mes.h"
#include "ghmm/sreestimate.h"
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_State.h"
#include "ghmm++/GHMM_Transition.h"
#include "ghmm++/GHMM_Emission.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_ContinuousModel::GHMM_ContinuousModel() {
  /* Do not initialize anything. It will be read from xml file. */
  c_model = NULL;
}


GHMM_ContinuousModel::GHMM_ContinuousModel(int N, int M, int cos, density_t density,
					   double prior) {
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


XMLIO_Element* GHMM_ContinuousModel::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "state") {
    GHMM_State* ghmm_state = new GHMM_State(this,states.size(),attrs);
    states.push_back(ghmm_state);

    state_by_id[ghmm_state->id] = states.size() - 1;

    /* Pass all nested elements to this state. */
    return ghmm_state;
  }

  if (tag == "transition") {
    GHMM_Transition* transition = new GHMM_Transition(attrs);
    transitions.push_back(transition);

    return transition;
  }

  fprintf(stderr,"tag '%s' not recognized in hmm element\n",tag.c_str());
  exit(1);
  
  return NULL;
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
    c_model->density = states[0]->emission->density;
  for (i = 1; i < states.size(); ++i)
    if (c_model->density != states[i]->emission->density) {
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
  for (i = 0; i < transitions.size(); ++i)
    SAFE_DELETE(transitions[i]);
  transitions.clear();
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
