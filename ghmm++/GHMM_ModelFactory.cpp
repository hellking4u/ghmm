/*
 * created: 29 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ghmm++/GHMM_ModelFactory.h"
#include "ghmm++/GHMM_DiscreteModel.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


vector<GHMM_DiscreteModel*> GHMM_ModelFactory::read_models(char *filename) {
  model** models;
  int mo_number;
  vector<GHMM_DiscreteModel*> v;

  models = model_read(filename,&mo_number);

  for (int i = 0; i < mo_number; ++i)
    v.push_back(new GHMM_DiscreteModel(models[i]));

  return v;
}


/**
   Produces simple left-right models given sequences. 
   The function "model_generate_from_sequence" is called for each 
   model that should be made. The sequences are read in from the
   ASCII file and thrown away again when leaving the function.
   @return vector of models
   @param s:          scanner
   @param new_models: number of models to produce */
//model **GHMM_ModelFactory::model_from_sequence_ascii(scanner_t *s, long *mo_number) {
//}


/** 
    Produces simple left-right models given sequences. The sequences
    are not read in from file, but exists already as a structur.
    @return vector of models
    @param s:          scanner
    @param new_models: number of models to produce */
//model **GHMM_ModelFactory::model_from_sequence(sequence_t *sq, long *mo_number) {
//}


/**
   Tests if number of states and number of outputs in the models match.
   @return 0 for succes; -1 for error
   @param mo:           vector of models
   @param model_number: numbr of models */
//int GHMM_ModelFactory::model_check_compatibility(model **mo, int model_number) {
//}


/**
   Reads one or several arrays of sequences. 
   @param filename    input filename
*/
//void GHMM_ModelFactory::sequence_read(char *filename) {
//  int number;

//  clean();

//  if (sequence_type == GHMM_INT)
//    c_t_sequence = sequence_read(filename,&number);
//  if (sequence_type == GHMM_DOUBLE)
//    c_d_sequence = sequence_d_read(filename,&number);
//}

