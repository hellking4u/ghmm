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
#include "ppghmm++/GHMM_Emission.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_State::GHMM_State(GHMM_AbstractModel* my_model, int my_index, XMLIO_Attributes& attrs) {
  index        = my_index;
  c_sstate     = NULL;
  reading      = GHMM_STATE_NONE;
  emission     = NULL;
  parent_model = my_model;

  /* by default take index as id. */
  id = attrs["id"];
  if (id == "") {
    char mem[100];
    snprintf(mem,sizeof(mem),"%d",my_index);
    id = string(mem);
  }
}


GHMM_State::GHMM_State(GHMM_AbstractModel* my_model, int my_index, sstate* my_state) {
  index        = my_index;
  c_sstate     = my_state;
  reading      = GHMM_STATE_NONE;
  emission     = NULL;
  parent_model = my_model;
  
  /* take index as id. */
  char mem[100];
  snprintf(mem,sizeof(mem),"%d",my_index);
  id = string(mem);
}


GHMM_State::~GHMM_State() {
  SAFE_DELETE(emission);
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
    SAFE_DELETE(emission);
    emission = new GHMM_Emission(this);

    return emission;
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


void GHMM_State::fillState(sstate* s) {
  GHMM_ContinuousModel* model = dynamic_cast<GHMM_ContinuousModel*>(parent_model);
  smodel* c_model             = model->c_model;

  /* store current c representation of state. */
  c_sstate = s;

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
  s->pi = initial;

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
    s->c[i]   = emission->weight;
    s->mue[i] = emission->mue;
    s->u[i]   = emission->variance;
  }

  s->fix = 0;
}


void GHMM_State::XMLIO_finishedReading() {
  if (!emission) {
    fprintf(stderr,"state with id='%s' lacks emission element.\n",id.c_str());
    exit(1);
  }
}


void GHMM_State::changeOutEdge(int target, double prob) {
  changeOutEdge(0,target,prob);
}


void GHMM_State::changeInEdge(int target, double prob) {
  changeInEdge(0,target,prob);
}


void GHMM_State::changeOutEdge(int matrix_index, int target, double prob) {
  if (prob == 0) {
    removeOutEdge(target);
    return;
  }
  
  int i;
  for (i = 0; i < c_sstate->out_states; ++i)
    if (c_sstate->out_id[i] == target) {
      c_sstate->out_a[matrix_index][i] = prob;
      return;
    }

  /* create new state. */
  c_sstate->out_states += 1;
  c_sstate->out_id = (int*) realloc(c_sstate->out_a,sizeof(int) * c_sstate->out_states);
  c_sstate->out_id[c_sstate->out_states - 1] = target;

  for (i = 0; i < parent_model->getNumberOfTransitionMatrices(); ++i) {
    c_sstate->out_a[i] = (double*) realloc(c_sstate->out_a,sizeof(double) * c_sstate->out_states);
    c_sstate->out_a[i][c_sstate->out_states - 1]  = 0;
  }
  c_sstate->out_a[matrix_index][c_sstate->out_states - 1]  = prob;
}


void GHMM_State::removeOutEdge(int target) {
  int i;
  for (i = 0; i < c_sstate->out_states; ++i)
    if (c_sstate->out_id[i] == target) {
      c_sstate->out_id[i] = c_sstate->out_id[c_sstate->out_states - 1];
      for (int m = 0; m < parent_model->getNumberOfTransitionMatrices(); ++m)
	c_sstate->out_a[m][i]  = c_sstate->out_a[m][c_sstate->out_states - 1];
      c_sstate->out_states -= 1;
      return;
    }
}


void GHMM_State::removeInEdge(int source) {
  int i;
  for (i = 0; i < c_sstate->in_states; ++i)
    if (c_sstate->in_id[i] == source) {
      c_sstate->in_id[i] = c_sstate->in_id[c_sstate->in_states - 1];
      for (int m = 0; m < parent_model->getNumberOfTransitionMatrices(); ++m)
	c_sstate->in_a[m][i]  = c_sstate->in_a[m][c_sstate->in_states - 1];
      c_sstate->in_states -= 1;
      return;
    }
}


void GHMM_State::changeInEdge(int matrix_index, int source, double prob) {
  if (prob == 0) {
    removeInEdge(source);
    return;
  }
  
  int i;
  for (i = 0; i < c_sstate->in_states; ++i)
    if (c_sstate->in_id[i] == source) {
      c_sstate->in_a[matrix_index][i] = prob;
      return;
    }

  /* create new state. */
  c_sstate->in_states += 1;
  c_sstate->in_id = (int*) realloc(c_sstate->in_a,sizeof(int) * c_sstate->in_states);
  c_sstate->in_id[c_sstate->in_states - 1] = source;

  for (i = 0; i < parent_model->getNumberOfTransitionMatrices(); ++i) {
    c_sstate->in_a[i] = (double*) realloc(c_sstate->in_a,sizeof(double) * c_sstate->in_states);
    c_sstate->in_a[i][c_sstate->in_states - 1]  = 0;
  }
  c_sstate->in_a[matrix_index][c_sstate->in_states - 1]  = prob;
}
