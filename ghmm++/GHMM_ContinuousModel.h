/*
 * created: 05 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_CONTINUOUSMODEL_H
#define _GHMM_CONTINUOUSMODEL_H 1

#include <vector>
#include <map>
#include <ghmm/smodel.h>
#include <ghmm++/GHMM_AbstractModel.h>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_ContinuousModel;
class GHMM_Sequences;
class GHMM_State;
class GHMM_Transition;

/** */
class GHMM_ContinuousModel: public GHMM_AbstractModel {

 public:

  /** Just for reading from xml file. c_data is left uninitialized. */
  GHMM_ContinuousModel();
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
  GHMM_ContinuousModel(int N, int M, int cos, density_t density, double prior=-1);
  /** Destructor. */
  virtual ~GHMM_ContinuousModel();

  /** Returns name of class. */
  virtual const char* toString() const;

  /**
     Tests if all standardization requirements of model are fulfilled. 
     (That is, if the sum of the probabilities is 1).
     @return 0 for succes; -1 for error. 
  */
  virtual int check() const;
  /* Returns state with given index. */
  sstate* getCState(int index) const;
  /** 
      Produces sequences to a given model. All memory that is needed for the 
      sequences is allocated inside the function. It is possible to define
      the length of the sequences global (global_len > 0) or it can be set 
      inside the function, when a final state in the model is reach (a state
      with no output). If the model has no final state, the sequences will
      have length MAX_SEQ_LEN.
      @return             pointer to an array of sequences
      @param seed:        initial parameter for the random value generator
                          (an integer). If seed == 0, then the random value
			  generator is not initialized.
      @param global_len:  length of sequences (=0: automatically via final states)
      @param seq_number:  number of sequences
      @param label:       label tag
      @param Tmax:        maximal sequence length, set to MAX_SEQ_LEN if -1 
  */
  GHMM_Sequences* generate_sequences(int seed, int global_len, long seq_number, long label, int Tmax=-1) const;
  /**
     Writes the model in matrix format.
     @param file: output file
  */
  virtual void print(FILE *file) const;
  /**
     Baum-Welch Algorithm for SHMMs.
     Training of model parameter with multiple double sequences (incl. scaling).
     New parameters set directly in hmm (no storage of previous values!).
     @return           0/-1 success/error
     @param seq        sequences used for training.
     @param logp       calculated log likelihood.
     @param eps        leave reestimation loop if diff. between successive logp values 
                       is smaller than eps.
     @param max_iter   max. no of iterations. 
  */
  int reestimate_baum_welch(GHMM_Sequences* seq, double* logp, double eps, int max_iter);
  /** Clean model. */
  virtual void clean();
  /** Copies c model into this object. */
  void copyFromModel(smodel* smo);
  /** */
  virtual int getNumberOfTransitionMatrices() const;
  /** */
  void read(const string& filename);

  /** C Model. */
  smodel* c_model;


 private:

  /** */
  virtual void XMLIO_finishedReading();
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  /** Build c++ data from c_model. */
  void buildCppData();
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_CONTINUOUSMODEL_H */
