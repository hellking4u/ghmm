/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ppghmm++/GHMM_Document.h"
#include "ppghmm++/GHMM_DiscreteModel.h"
#include "ppghmm++/GHMM_ContinuousModel.h"
#include "ppghmm++/GHMM_Sequences.h"


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
	continuous_model = new GHMM_ContinuousModel();
	return continuous_model;
      }
     
      /* error message if no valid hmm type is specified. */
      fprintf(stderr,"Error reading file: %s\n",xml_filename.c_str());
      if (attrs["type"] == "")
	fprintf(stderr,"No HMM type specified, ");
      else
	fprintf(stderr,"HMM type '%s' not recognized, ",attrs["type"].c_str());

      fprintf(stderr,"choose out of: continuous\n");

      exit(1);
    }
    
    if (tag == "sequences") {
      fprintf(stderr,"<sequences>\n");
      
      sequences = new GHMM_Sequences(attrs);
      
      return sequences;
    }
    
    fprintf(stderr,"Tag '%s' not recognized in ghmm element.\n",tag.c_str());
    exit(1);
  }
  
  return this;
}


void GHMM_Document::XMLIO_endTag(const string& tag) {
  if (tag == "ghmm")
    reading_ghmm = false;
}
