/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <iostream>
#include <assert.h>

#include <xmlio/XMLIO_Definitions.h>
#include <xmlio/XMLIO_Document.h>
#include "ghmm/ghmm.h"
#include "ghmm/matrix.h"
#include <ghmm++/GHMM_GMLTransition.h>
#include "ghmm++/GHMM_GMLEmission.h"
#include "ghmm++/GHMM_StateT.hh"
#include "ghmm++/GHMM_GMLDataNode.h"
#include "ghmm++/GHMM_SWDiscreteModel.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLState::GHMM_GMLState(GHMM_SWDiscreteModel* my_model, int my_index, XMLIO_Attributes& attrs) :
  GHMM_StateT<sdstate, GHMM_GMLTransition, GHMM_SWDiscreteModel>(my_model, my_index, attrs)
{
  emission          = NULL;
  tag               = "node"; // state
  hasData["ngeom"]  = 0;
  m_countme         = 0;
}

/************
GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, sstate* my_state) :
GHMM_StateT<sdstate, GHMM_GMLTransition>(my_model, my_index, my_state)
{
  tag               = "node"; // state
  hasData["ngeom"]  = false;
}

GHMM_GMLState::GHMM_GMLState(GHMM_AbstractModel* my_model, int my_index, state* my_state) :
GHMM_StateT<sdstate, GHMM_GMLTransition>(my_model, my_index, my_state)
{
  tag               = "node"; // state
  hasData["ngeom"]  = false;
}****************/


GHMM_GMLState::GHMM_GMLState(GHMM_SWDiscreteModel* my_model, int my_index, sdstate* my_state) :
GHMM_StateT<sdstate, GHMM_GMLTransition, GHMM_SWDiscreteModel>(my_model, my_index, my_state)
{
  emission          = NULL;
  tag               = "node"; // state
  hasData["ngeom"]  = 0;
  m_countme         = 0;
}


GHMM_GMLState::~GHMM_GMLState() {
  if (emission != NULL) {
    delete emission;
  }
}

const char* GHMM_GMLState::toString() const {
  return "GHMM_GMLState";
}


void GHMM_GMLState::fillState(sdstate* s) {
  GHMM_SWDiscreteModel* m = (GHMM_SWDiscreteModel*) parent_model;
  sdmodel* c_model        = m->c_model;

  /* store current c representation of state. */
  c_sdstate = s;

  s->b   = (double*) malloc(sizeof(double) * c_model->M);

  vector<GHMM_GMLTransition*> out_edges;
  vector<GHMM_GMLTransition*> in_edges;
  unsigned int i;
  
  for (i = 0; i < m->transitions.size(); ++i) {
    if (m->transitions[i]->source == id)
      out_edges.push_back(m->transitions[i]);
    if (m->transitions[i]->target == id)
      in_edges.push_back(m->transitions[i]);
  }

  if (out_edges.size() > 0) {
    assert ( out_edges.size() <= c_model->N );
    s->out_id = (int*) malloc(sizeof(int) * out_edges.size());
    s->out_a  = matrix_d_alloc(c_model->cos, out_edges.size());
  }

  if (in_edges.size() > 0) {
    assert ( in_edges.size() <= c_model->N );
    s->in_id = (int*) malloc(sizeof(int) * in_edges.size());
    s->in_a  = matrix_d_alloc(c_model->cos, in_edges.size());
  }

  /* now fill with useful data. */
  s->pi = initial;

  for (i = 0; i < out_edges.size(); ++i)
    s->out_id[i] = m->getStateIndex(out_edges[i]->target);

  for (i = 0; i < in_edges.size(); ++i)
    s->in_id[i] = m->getStateIndex(in_edges[i]->source);

  int cos;
  for (cos = 0; cos < c_model->cos; ++cos) {
    for (i = 0; i < out_edges.size(); ++i)
      s->out_a[cos][i] = out_edges[i]->probs[cos];

    for (i = 0; i < in_edges.size(); ++i)
      s->in_a[cos][i] = in_edges[i]->probs[cos];
  }

  s->out_states = out_edges.size();
  s->in_states  = in_edges.size();

  if (c_model->M != (int) emission->weights.size()) {
    fprintf(stderr,"M == %d, but just %d weights found in GHMM_State.cpp\n",c_model->M,(int) emission->weights.size());
    exit(1);
  }

  int silent_flag=0;
  for (i = 0; (int) i < c_model->M; ++i)
    {
      s->b[i]=emission->weights[i];
      silent_flag= silent_flag || (s->b[i] == 0.0);
    }
  if (silent_flag) c_model->model_type=kSilentStates;
  s->countme = m_countme;
  s->fix = 0;

  s->label = (char *)malloc(strlen(label.c_str()) + 1);
  strcpy(s->label, label.c_str());
}


GHMM_GMLTransition* GHMM_GMLState::createTransition(int edge_index, sdmodel *mo) {
  if (c_sdstate)
    {	 
      vector<double> vtprobs;
      for(int i=0; i < mo->cos; i++) {
	vtprobs.push_back( c_sdstate->out_a[i][edge_index] );
      }
      return new GHMM_GMLTransition(this,parent_model->getState(c_sdstate->out_id[edge_index]), vtprobs);
    }
  
  return NULL;
}


void GHMM_GMLState::setOutputProbability(int index, double prob) {
  if (!c_sdstate) {
    fprintf(stderr,"GHMM_State::setOutputProbability(int,double): object does not contain state.\n");
    exit(1);
  } else
    {
      GHMM_SWDiscreteModel* m = (GHMM_SWDiscreteModel*) parent_model;
      if (index >= m->c_model->M) {
	fprintf(stderr,"GHMM_State::setOutputProbability(int,double): symbol No. %d requested, but maximum number is %d.\n",index,m->c_model->M - 1);
	exit(1);
      }

      c_sdstate->b[index] = prob;
    }

  /******************
  if (!c_state) {
    fprintf(stderr,"GHMM_State::setOutputProbability(int,double): object does not contain state.\n");
    exit(1);
  } else
    {
      GHMM_SWDiscreteModel* m = (GHMM_SWDiscreteModel*) parent_model;
      if (index >= m->c_model->M) {
	fprintf(stderr,"GHMM_State::setOutputProbability(int,double): symbol No. %d requested, but maximum number is %d.\n",index,m->c_model->M - 1);
	exit(1);
      }

      c_state->b[index] = prob;
    }
  ***************/
}

void GHMM_GMLState::setOutputProbability(const string& id, double prob) {
  setOutputProbability(((GHMM_SWDiscreteModel*)(parent_model))->alphabet->getIndex(id),prob);
}


XMLIO_Element* GHMM_GMLState::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  
  if (tag == "data") /* For GraphML */
    {
      if (attrs["key"] == "initial") {
	reading = GHMM_STATE_INITIAL;
	hasData["initial"] = 1;
	return this;
      }

      if (attrs["key"] == "label") {
	reading = GHMM_STATE_LABEL;
	hasData["label"] = 1;
	return this;
      }

      if (attrs["key"] == "countme") {
	reading = GHMM_STATE_COUNTME;
	hasData["countme"] = 1;
	return this;
      }

      if (attrs["key"] == "emissions") {
		hasData["emissions"] = 1;
		attributes["key"] = "emissions";
		SAFE_DELETE(emission);
		emission = new GHMM_GMLEmission(this);
		return emission;
      }

      if (attrs["key"] == "ngeom") {
		attributes["key"] = "ngeom";
		return this;
      }
    }
  else
    if (tag == "pos")
      {
	hasData["ngeom"] = 1;
	vPosition[0] = (float)atof(attrs["x"].c_str());
	vPosition[1] = (float)atof(attrs["y"].c_str());
	return this;
      }
    else
    { 
      fprintf(stderr,"GMLState :: tag '%s' not recognized in state element\n",tag.c_str());
      exit(1);
    }

  return NULL;
}


void GHMM_GMLState::XMLIO_endTag(const string& tag) {
  reading = GHMM_STATE_NONE;
}

float GHMM_GMLState::get2DPosition(int index) {
  return vPosition[index];
}


void GHMM_GMLState::XMLIO_getCharacters(const string& characters) {
  switch (reading) {

  case GHMM_STATE_INITIAL:
    initial = atof(characters.c_str());
    break;

  case GHMM_STATE_LABEL:
    label = characters.c_str();
    break;

  case GHMM_STATE_COUNTME:
    m_countme = atoi(characters.c_str());
    break;

  case GHMM_STATE_NONE:
    break;
  }
}

XMLIO_Attributes& GHMM_GMLState::XMLIO_getAttributes() {
  attributes["id"] = id;
  return attributes;
}


const int GHMM_GMLState::XMLIO_writeContent(XMLIO_Document& writer) {
  int total_bytes = 0;
  int result = 0;
  
  double initial = getInitial();

  writer.changeIndent(2);

  //  if (initial > 0) {
  writer.writeEndlIndent();
  
  if (result < 0)
    return result;
  total_bytes += result;
  
  result = writer.writef("<data key=\"initial\">%f</data>",initial);
  if (result < 0)
    return result;
  total_bytes += result;

  result = writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;

  GHMM_GMLDataNode* my_label = new GHMM_GMLDataNode(this, "label");
  result = writer.writeElement(my_label);
  SAFE_DELETE(my_label);
  result += writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;

  //  }

  if ( hasData["ngeom"] )
    {
      GHMM_GMLDataNode* my_pos = new GHMM_GMLDataNode(this, "ngeom");
      result = writer.writeElement(my_pos);
      SAFE_DELETE(my_pos);
    }
  else
    {
      // Warning or put in some random positions
      // HMMEd expected a file to contain coordinate, or else we should change it
    }

  /* write emission */
  result += writer.writeEndl();
  if (result < 0)
    return result;
  total_bytes += result;
  
  GHMM_GMLDataNode* my_emission = new GHMM_GMLDataNode(this, "emissions");
  result = writer.writeElement(my_emission);
  SAFE_DELETE(my_emission);
  
  if (result < 0)
    return result;
  total_bytes += result;
  /* end write emission */

  result = writer.writeEndl();

  if (result < 0)
    return result;
  total_bytes += result;

  return total_bytes;
}








