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

#include "ghmm++/GHMM_GMLDoc.h"
#include "ghmm++/GHMM_GMLDiscreteModel.h"
#include "ghmm++/GHMM_GMLContinuousModel.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_GMLAlphabet.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GraphMLDoc::GHMM_GraphMLDoc() {

  discrete_model    = NULL;
  continuous_model  = NULL;

  sequences        = NULL;
  reading_ghmm     = false;
  tmp_alphabets    = NULL;
}

GHMM_GraphMLDoc::~GHMM_GraphMLDoc() {
  //  SAFE_DELETE(discrete_model);
  //  SAFE_DELETE(continuous_model);
}


const char* GHMM_GraphMLDoc::toString() const {
  return "GHMM_GraphMLDoc";
}


GHMM_GMLDiscreteModel* GHMM_GraphMLDoc::getDiscreteModel() const {
  return discrete_model;
}


GHMM_GMLContinuousModel* GHMM_GraphMLDoc::getContinuousModel() const {
  return continuous_model;
}


GHMM_Sequences* GHMM_GraphMLDoc::getSequences() const {
  return sequences;
}

XMLIO_Element* GHMM_GraphMLDoc::XMLIO_startTag(const string& my_tag, XMLIO_Attributes &attrs) {

  if (my_tag == "graphml") {
    reading_ghmm = true;
printf("GHMM_GraphMLDoc::XMLIO_startTag\n");
  cout << "\t\t" << my_tag << endl;

    return this;
  }

  if (reading_ghmm) {

    if (my_tag == "desc") {
      return this;
    }

    if (my_tag == "hmm:class")
     { return this; }

    if (my_tag == "paint")
      {return this;
      }
    if (my_tag == "point" || my_tag == "line")
      {  return this;}
    
    if (my_tag == "map" || my_tag == "symbol")
      { return this; }

    if (my_tag == "hmm:alphabet") {
printf("GHMM_GraphMLDoc::XMLIO_startTag\n");
  cout << "\t\t" << my_tag << endl;


      tmp_alphabets = new GHMM_GMLAlphabet();
      return tmp_alphabets;
    }

    if (my_tag == "key") {
printf("GHMM_GraphMLDoc::XMLIO_startTag\n");
  cout << "\t\t" << my_tag << endl;

       if (attrs["id"] == "emissions")
	 {
	   if ( attrs["gd:type"] == "HigherDiscreteProbDist" )
	     {    
	       model_type = GHMM_DISCRETE;
	       return this;
	     }
	   
	   else if ( attrs["gd:type"] == "ContinuousProbDist" )
	     {
	       model_type = GHMM_CONTINUOUS;	
	       return this; 
	     } else
	       {
		 /* error message if no valid hmm type is specified. */
		 fprintf(stderr, "Need to know the type of the HMM\n");	     
		 exit(-1);
	       }
	 } else 
	   return this;
    }
    

    if (my_tag == "graph")
      {
	printf("GHMM_GraphMLDoc::XMLIO_startTag:"); cout << my_tag << endl;

	if ( model_type != NONE )
	  {
	    if ( model_type == GHMM_DISCRETE )
	      {
		printf("Discrete model found\n");

		GHMM_Alphabet *alphas = tmp_alphabets;
		assert( alphas != NULL );
		for(int i=0; i < alphas->size(); i++)
		  {
		    cout << "Symbol " << i << ":";
		    cout << alphas->getSymbol(i) << endl;
		  }

		discrete_model = new GHMM_GMLDiscreteModel(tmp_alphabets);
		return discrete_model;
	      }
	    
	    if ( model_type == GHMM_CONTINUOUS )
	      {
		printf("Continous model found\n");		
		continuous_model = new GHMM_GMLContinuousModel();
		return continuous_model;
	      }
	  }
	else
	  {
	    /* error message if no valid hmm type is specified. */
	    fprintf(stderr, "Need to know the type of the HMM\n");
	    exit(-1);
	  }
      } else
	{
	  fprintf(stderr,"Tag '%s' not recognized in graphml element.\n",my_tag.c_str());
	  exit(1);
	}
  }
  return NULL;
}


void GHMM_GraphMLDoc::XMLIO_endTag(const string& my_tag) {
  if (my_tag == "graphml")
    reading_ghmm = false;
}


int GHMM_GraphMLDoc::XMLIO_writeTrailer() {
  return writef("</graphml>\n");
}


int GHMM_GraphMLDoc::XMLIO_writeProlog() {
  int this_result;
  int return_result = 0;

  this_result = XMLIO_Document::XMLIO_writeProlog();

  /* Returns error code if an error occured. */
  if (this_result < 0)
    return this_result;

  return_result += this_result;

  this_result = writef("<graphml version=\"1.0\">\n");
     
  /* Returns error code if an error occured. */
  if (this_result < 0)
    return this_result;
  return_result += this_result;

  changeIndent(2);

  XMLIO_Element modeltype;
  char tmp[15];
  modeltype.tag = "key";
  modeltype.attributes["for"] ="node";

  switch( model_type )
    {
    case GHMM_DISCRETE:
      modeltype.attributes["gd:type"] = "HigherDiscreteProbDist";
      break;
    case GHMM_CONTINUOUS:
      modeltype.attributes["gd:type"] = "ContinuousProbDist";
      break;
    }
  modeltype.attributes["id"]      = "emissions";

  writeEndl();

  this_result = writeElement(&modeltype);
  writeEndl();
  /* Returns error code if an error occured. */
  if (this_result < 0)
    return this_result;
  return_result += this_result;
      
  return return_result;
}
