/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ppghmm++/GHMM_Document.h"
#include "ppghmm++/GHMM_DiscreteModel.h"
#include "ppghmm++/GHMM_ContinuousModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Document::GHMM_Document() {
  discrete_model   = NULL;
  continuous_model = NULL;
  sequences        = NULL;
  reading_ghmm     = false;
}


GHMM_Document::~GHMM_Document() {
  //  SAFE_DELETE(discrete_model);
  //  SAFE_DELETE(continuous_model);
}


const char* GHMM_Document::toString() const {
  return "GHMM_Document";
}


GHMM_DiscreteModel* GHMM_Document::getDiscreteModel() const {
  return discrete_model;
}


GHMM_ContinuousModel* GHMM_Document::getContinuousModel() const {
  return continuous_model;
}


GHMM_Sequences* GHMM_Document::getSequences() const {
  return sequences;
}


XMLIO_Element* GHMM_Document::XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs) {
  if (tag == "ghmm") {
    reading_ghmm = true;

    return this;
  }

  if (reading_ghmm) {
    if (tag == "hmm") {
      if (attrs["type"] == "continuous") {
	continuous_model = new GHMM_ContinuousModel(attrs);
	return continuous_model;
      }
      
      fprintf(stderr,"hmm type '%s' not recognized\n",attrs["type"].c_str());
      exit(1);
    }

    fprintf(stderr,"tag '%s' not recognized in ghmm element\n",tag.c_str());
    exit(1);
  }

  return this;
}


void GHMM_Document::XMLIO_endTag(const string& tag) {
  if (tag == "ghmm")
    reading_ghmm = false;
}
