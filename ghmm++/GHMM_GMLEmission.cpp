/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#include <xmlio/XMLIO_Document.h>

#include "ghmm++/GHMM_Emission.h"
#include "ghmm++/GHMM_Alphabet.h"
#include "ghmm++/GHMM_ContinuousModel.h"
#include "ghmm++/GHMM_DiscreteModelT.hh"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_GMLEmission::~GHMM_GMLEmission() {
}


const int GHMM_GMLEmission::XMLIO_writeContent(XMLIO_Document& writer) {
  int result = 0;
  int i;

  writer.changeIndent(2);

  /* continuous model */
  /*******
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
    default:
      break;
    }
    result += writer.writef(" mue=\"%f\" variance=\"%f\">",mue[0],variance[0]);
    } else 
  **********/
  {
	 /* discrete model */
    GHMM_Alphabet* alphabet      = state->getModel()->getAlphabet();
    GHMM_DiscreteModelT*model   = (GHMM_DiscreteModelT*)(state->getModel());

    result += writer.writeEndl();
    for (i = 0; i < model->c_model->M; ++i) {
      result += writer.writef("%s%.2f",writer.indent,state->getEmissionFrom(i)); // c_state->b[i]
      if (alphabet)
	result += writer.writef(" <!-- %s -->",alphabet->getSymbol(i).c_str());
      result += writer.writeEndl();
    }
  }
  
  return result;
}

