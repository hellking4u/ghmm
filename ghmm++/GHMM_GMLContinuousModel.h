/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */


#ifndef _GHMM_GMLCONTINUOUSMODEL_H
#define _GHMM_GMLCONTINUOUSMODEL_H 1

#include "ghmm++/GHMM_ContinuousModel.h"

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif


class GHMM_GMLContinuousModel;

class GHMM_GMLContinuousModel : public GHMM_ContinuousModel {

 public:
  
  /** Just for reading from xml file. c_data is left uninitialized. */
  GHMM_GMLContinuousModel();
  /** 
      Constructor.
      @param N       Number of states
      @param M       Number of output densities per state
      @param cos     smodel includes continuous model with one transition matrix 
                     (cos  is set to 1) and an extension for models with several matrices
                     (cos is set to a positive integer value > 1).
      @param density Flag for density function. 0: normal density, 1: truncated normal 
                     density, 2: approximated normal density.
      @param prior   prior for a priori probability of the model. -1 means no prior specified (all
                     models have equal probability a priori. 
  */
  GHMM_GMLContinuousModel(int N, int M, int cos, density_t density, double prior=-1);

};


#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* */
