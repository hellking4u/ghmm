/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include <xmlio/XMLIO_Document.h>
#include "ghmm++/GHMM_Emission.h"
#include "ghmm++/GHMM_State.h"
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_DiscreteModel.h"
#include "ghmm++/GHMM_Alphabet.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Emission::GHMM_Emission(GHMM_State* my_state) {
  int i;

  state             = my_state;
  tag               = "emission";
  xmlio_indent_type = XMLIO_INDENT_BOTH;

  switch (getModelType()) {

  case GHMM_CONTINOUS:
    /* only read data if state has been initialized */
    if (my_state->c_sstate) {
      mue.push_back(my_state->c_sstate->mue[0]);
      variance.push_back(my_state->c_sstate->u[0]);
      weights.push_back(my_state->c_sstate->c[0]);
      density  = ((GHMM_ContinuousModel*)(my_state->getModel()))->c_model->density;
    }
    break;

  case GHMM_DISCRETE:
    /* only read data if state has been initialized */
    if (my_state->c_state)
      for (i = 0; i < ((GHMM_DiscreteModel*)(my_state->getModel()))->c_model->M; ++i)
	weights.push_back(my_state->c_state->b[i]);
    break;
  }
}


GHMM_Emission::~GHMM_Emission() {
}


const char* GHMM_Emission::toString() const {
  return "GHMM_Emission";
}


XMLIO_Element* GHMM_Emission::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  bool found = false;

  switch (getModelType()) {

  case GHMM_CONTINOUS:
    if (tag == "gauss") {
      mue.push_back(atof(attrs["mue"].c_str()));
      variance.push_back(atof(attrs["variance"].c_str()));
      density  = normal;
      found    = true;
    }
    if (tag == "gauss-positive") {
      mue.push_back(atof(attrs["mue"].c_str()));
      variance.push_back(atof(attrs["variance"].c_str()));
      density  = normal_pos;
      found    = true;
    }
    if (tag == "gauss-approximated") {
      mue.push_back(atof(attrs["mue"].c_str()));
      variance.push_back(atof(attrs["variance"].c_str()));
      density  = normal_approx;
      found    = true;
    }
    break;
  }

  if (! found) {
    fprintf(stderr,"<emission> element of state '%s' has unrecognized tag '%s'.\n",
	    state->id.c_str(),tag.c_str());
    exit(1);
  }

  if (attrs["mue"] == "") {
    fprintf(stderr,"<%s> element of state '%s' lacks mue attribute.\n",
	    tag.c_str(),state->id.c_str());
    exit(1);
  }

  if (attrs["variance"] == "") {
    fprintf(stderr,"<%s> element of state '%s' lacks variance attribute.\n",
	    tag.c_str(),state->id.c_str());
    exit(1);
  }

  return NULL;
}


void GHMM_Emission::XMLIO_getCharacters(const string& characters) {
  unsigned int pos;
  for (pos = 0; pos < characters.size(); ++pos) {
    while (pos < characters.size() && isspace(characters[pos]))
      ++pos;

    if (pos < characters.size())
      weights.push_back(atof(characters.substr(pos).c_str()));

    while (pos < characters.size() && !isspace(characters[pos]))
      ++pos;
  }
}


void GHMM_Emission::XMLIO_finishedReading() {
  /* continuous model */
  if (state->c_sstate)
    if (weights.size() != mue.size()) {
      fprintf(stderr,"Different number of weights and density functions in state '%s'.\n",state->id.c_str());
      exit(1);
    }
}


const int GHMM_Emission::XMLIO_writeContent(XMLIO_Document& writer) {
  int result = 0;
  int i;

  writer.changeIndent(2);

  /* continuous model */
  if (state->c_sstate) {
    result = writer.writef("1 <");
    switch (density) {
    case normal:
      result += writer.write("gauss");
      break;
    case normal_pos:
      result += writer.write("gauss-positive");
      break;
    case normal_approx:
      result += writer.write("gauss-approximated");
      break;
    }
    result += writer.writef(" mue=\"%f\" variance=\"%f\">",mue[0],variance[0]);
  }

  /* discrete model */
  if (state->c_state) {
    GHMM_Alphabet* alphabet = state->getModel()->getAlphabet();
    GHMM_DiscreteModel* model = (GHMM_DiscreteModel*) state->getModel();

    result += writer.writeEndl();
    for (i = 0; i < model->c_model->M; ++i) {
      result += writer.writef("%s%.2f",writer.indent,state->c_state->b[i]);
      if (alphabet)
	result += writer.writef(" <!-- %s -->",alphabet->getSymbol(i).c_str());
      result += writer.writeEndl();
    }
  }

  return result;
}


GHMM_ModelType GHMM_Emission::getModelType() const {
  return state->getModelType();
}
