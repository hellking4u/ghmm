/*
 * created: 29 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_MODELFACTORY_H
#define _GHMM_MODELFACTORY_H 1

#include <vector>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_ModelFactory;
class GHMM_DiscreteModel;

/** */
class GHMM_ModelFactory {

 public:

  /**
     Reads in ASCII data to initialize a vector of models.
     @return vector of pointers to the models
     @param filename:   the ASCII input file
  */
  static vector<GHMM_DiscreteModel*> read_models(char *filename);

};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_MODELFACTORY_H */
