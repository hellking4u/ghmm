/*
 * created: 21 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include <xmlio/XMLIO_Definitions.h>
#include "ghmm/matrix.h"
#include "ppghmm++/GHMM_State.h"
#include "ppghmm++/GHMM_Transition.h"
#include "ppghmm++/GHMM_ContinuousModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_State::GHMM_State(XMLIO_Attributes& attrs) {
  c_state = NULL;
  reading = GHMM_STATE_NONE;
  id      = attrs["id"];
}


GHMM_State::~GHMM_State() {
  if (c_state) {
    state_clean(c_state);
    free(c_state);
  }
}


const char* GHMM_State::toString() const {
  return "GHMM_State";
}


XMLIO_Element* GHMM_State::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "initial") {
    reading = GHMM_STATE_INITIAL;
    
    return this;
  }

  if (tag == "emission") {
    emission_weight   = atoi(attrs["weight"].c_str());
    emission_mue      = atof(attrs["mue"].c_str());
    emission_variance = atof(attrs["variance"].c_str());

    return this;
  }

  fprintf(stderr,"tag '%s' not recognized in state element\n",tag.c_str());
  exit(1);
  
  return NULL;
}


void GHMM_State::XMLIO_endTag(const string& tag) {
  reading = GHMM_STATE_NONE;
}


void GHMM_State::XMLIO_getCharacters(const string& characters) {
  switch (reading) {

  case GHMM_STATE_INITIAL:
    initial = atof(characters.c_str());
    break;
    
  case GHMM_STATE_NONE:
    break;
  }
}


void GHMM_State::fillState(GHMM_ContinuousModel* model, sstate* s) {
  smodel* c_model = model->c_model;

  s->c   = (double*) malloc(sizeof(double) * c_model->M);
  s->mue = (double*) malloc(sizeof(double) * c_model->M);
  s->u   = (double*) malloc(sizeof(double) * c_model->M);

  vector<GHMM_Transition*> out_edges;
  vector<GHMM_Transition*> in_edges;
  unsigned int i;
  for (i = 0; i < model->transitions.size(); ++i) {
    if (model->transitions[i]->source == id)
      out_edges.push_back(model->transitions[i]);
    if (model->transitions[i]->target == id)
      in_edges.push_back(model->transitions[i]);
  }

  if (out_edges.size() > 0) {
    s->out_id = (int*) malloc(sizeof(int) * out_edges.size());
    s->out_a  = matrix_d_alloc(c_model->cos, out_edges.size());
  }

  if (in_edges.size() > 0) {
    s->in_id = (int*) malloc(sizeof(int) * in_edges.size());
    s->in_a  = matrix_d_alloc(c_model->cos, in_edges.size());
  }

  /* now fill with useful data. */
  s->pi         = initial;

  for (i = 0; i < out_edges.size(); ++i)
    s->out_id[i] = model->getStateID(out_edges[i]->target);

  for (i = 0; i < in_edges.size(); ++i)
    s->in_id[i] = model->getStateID(in_edges[i]->source);

  int cos;
  for (cos = 0; cos < c_model->cos; ++cos) {
    for (i = 0; i < out_edges.size(); ++i)
      s->out_a[cos][i] = out_edges[i]->prob;

    for (i = 0; i < in_edges.size(); ++i)
      s->in_a[cos][i] = in_edges[i]->prob;
  }

  s->out_states = out_edges.size();
  s->in_states  = in_edges.size();

  if (c_model->M != 1) {
    fprintf(stderr,"M != 1 not yet supported in GHMM_State.cpp\n");
    exit(1);
  }

  for (i = 0; (int) i < c_model->M; ++i) {
    s->c[i]   = emission_weight;
    s->mue[i] = emission_mue;
    s->u[i]   = emission_variance;
  }

  s->fix        = 0;
}
