/*
 * created: 14 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ghmm++/GHMM_Document.h"
#include "ghmm++/GHMM_DiscreteModel.h"
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_Sequences.h"

#include <iostream>

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
  cerr << "Delete " << toString() << endl;
  //SAFE_DELETE(discrete_model);
  //SAFE_DELETE(continuous_model);
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


XMLIO_Element* GHMM_Document::XMLIO_startTag(const string& my_tag, XMLIO_Attributes &attrs) {
  if (my_tag == "ghmm") {
    reading_ghmm = true;

    return this;
  }

  if (reading_ghmm) {
    if (my_tag == "hmm") {
      if (attrs["type"] == "continuous") {
	continuous_model = new GHMM_ContinuousModel();
	return continuous_model;
      }

      if (attrs["type"] == "discrete") {
	discrete_model = new GHMM_DiscreteModel();
	return discrete_model;
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
    
    if (my_tag == "sequences") {
      sequences = new GHMM_Sequences(attrs);
      
      return sequences;
    }
    
    fprintf(stderr,"Tag '%s' not recognized in ghmm element.\n",my_tag.c_str());
    exit(1);
  }
  
  return this;
}


void GHMM_Document::XMLIO_endTag(const string& my_tag) {
  if (my_tag == "ghmm")
    reading_ghmm = false;
}


int GHMM_Document::XMLIO_writeTrailer() {
  return writef("</ghmm>\n");
}


int GHMM_Document::XMLIO_writeProlog() {
  int this_result;
  int return_result = 0;

  this_result = XMLIO_Document::XMLIO_writeProlog();

  /* Returns error code if an error occured. */
  if (this_result < 0)
    return this_result;

  return_result += this_result;

  this_result = writef("<ghmm version=\"1.0\">\n");

  /* Returns error code if an error occured. */
  if (this_result < 0)
    return this_result;
  return_result += this_result;

  changeIndent(2);

  return return_result;
}
