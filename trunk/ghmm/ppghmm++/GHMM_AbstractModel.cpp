/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ppghmm++/GHMM_AbstractModel.h"
#include "ppghmm++/GHMM_Transition.h"
#include "ppghmm++/GHMM_State.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_AbstractModel::GHMM_AbstractModel() {
}


GHMM_AbstractModel::~GHMM_AbstractModel() {
}


const char* GHMM_AbstractModel::toString() const {
  return "GHMM_AbstractModel";
}


void GHMM_AbstractModel::print(FILE *file) {
}


int GHMM_AbstractModel::check() {
  return 0;
}


void GHMM_AbstractModel::setTransition(GHMM_Transition* transition) {
  GHMM_State* source = getState(transition->source);
  GHMM_State* target = getState(transition->target);

  source->changeOutEdge(target->index,transition->prob);
  target->changeInEdge(source->index,transition->prob);
}


void GHMM_AbstractModel::setTransition(int start, int end, double prob) {
  GHMM_State* source = getState(start);
  GHMM_State* target = getState(end);

  source->changeOutEdge(target->index,prob);
  target->changeInEdge(source->index,prob);
}


GHMM_State* GHMM_AbstractModel::getState(const string& id) const {
  map<string,int>::const_iterator iter;
  iter = state_by_id.find(id);
  if (iter == state_by_id.end()) {
    fprintf(stderr,"State '%s' not found in model.\n",id.c_str());
    return NULL;
  }
  
  return states[iter->second];
}


GHMM_State* GHMM_AbstractModel::getState(int index) const {
  if (index >= (int) states.size()) {
    fprintf(stderr,"GHMM_AbstractModel::getState(int):\n");
    fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,(int) states.size());
    exit(1);
  }

  return states[index];
}


int GHMM_AbstractModel::getNumberOfTransitionMatrices() const {
  return 0;
}


void GHMM_AbstractModel::clean() {
  unsigned int i;
  for (i = 0; i < states.size(); ++i)
    SAFE_DELETE(states[i]);
  states.clear();

  for (i = 0; i < transitions.size(); ++i)
    SAFE_DELETE(transitions[i]);
  transitions.clear();
}


int GHMM_AbstractModel::getStateID(const string& id) {
  return state_by_id[id];
}
