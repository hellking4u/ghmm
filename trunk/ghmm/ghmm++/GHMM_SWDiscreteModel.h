/*
  created: 2. April 2003
  authors: Wasinee Rungsarityotin
  file   : $Source$
  $Id$

  __copyright__
*/

#ifndef _GHMM_SWDISCRETEMODEL_H
#define _GHMM_SWDISCRETEMODEL_H 1

#include <ghmm/model.h>
#include <ghmm/sdmodel.h>

#include <ghmm++/GHMM_StateT.hh>
#include <ghmm++/GHMM_GMLTransition.h>
#include <ghmm++/GHMM_Alphabet.h>

#include <ghmm++/GHMM_AbstractModelT.hh> // Template

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Sequences;
class GHMM_IntVector;
class GHMM_DoubleVector;
class GHMM_DoubleMatrix;
class GHMM_Alphabet;
class GHMM_GMLState;
class GHMM_GMLTransition;

/** Discrete HMM model (wrapper around model in C data structure). */
class GHMM_SWDiscreteModel : public GHMM_AbstractModelT<GHMM_GMLState, GHMM_GMLTransition> {

public:

  GHMM_SWDiscreteModel(GHMM_Alphabet* my_alphabet); 


  GHMM_SWDiscreteModel(GHMM_Alphabet* my_alphabet, int no_klass);


  /** Constructor. */
  GHMM_SWDiscreteModel( sdmodel* my_model);


    /** Destructor. */

  ~GHMM_SWDiscreteModel();

  /** Returns name of class. */
  const char* toString() const;

 /**
     Writes the model in matrix format.
     @param file: output file
  */
  void print(FILE *file);

   /** 
      Writes transition matrix of a model.
      @param file: output file
      @param tab:  format: leading tabs
      @param separator: format: seperator for columns
      @param ending:    format: end of a row 
  */
  void A_print(FILE *file, char *tab, char *separator, char *ending) const;
   
  /**
     Writes output matrix of a model.
     @param file: output file
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  void B_print(FILE *file, char *tab, char *separator, char *ending) const;
    
  void Pi_print(FILE *file, char *tab, char *separator, char *ending) const;
  
  /**
     Tests if all standardization requirements of model are fulfilled. 
     (That is, if the sum of the probabilities is 1).
     @return 0 on succes; -1 on error. 
  */
  int check() const { return 0; }

  GHMM_SWDiscreteModel* copy();

  GHMM_Sequences* generate_sequences(int seed, int global_len, long seq_number) const
  {  return NULL; }

  GHMM_Sequences* generate_sequences(int seed, int global_len, long seq_number, int Tmax) const;

  /** Returns alphabet of model or NULL, if no such alphabet exists. */
  GHMM_Alphabet* getAlphabet() const;

  /** Returns model type. */
  GHMM_ModelType getModelType() const;

  /** */
  int getNumberOfTransitionMatrices() const;

   /** */
  void XMLIO_finishedReading();
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  void createTransitions();

  const int      XMLIO_writeContent(XMLIO_Document& writer);

  /* Returns state with given index. */
  sdstate* getCState(int index) const;

  /** Forward-Algorithm.
      Calculates alpha[t][i], scaling factors scale[t] and log( P(O|lambda) ) for
      a given double sequence and a given model.
      @param seq:    sequences
      @param index:  index of sequence to take
      @param alpha:  alpha[t][i]
      @param scale:  scale factors (return value, if scale is NULL this vector will not be returned)
      @param log\_p:  log likelihood log( P(O|lambda) )
      @return NULL on error, alpha matrix on success
  */
  GHMM_DoubleMatrix* foba_forward(GHMM_Sequences* seq, int index, GHMM_DoubleVector* scale, 
				  double *log_p) const;

  /**
     Viterbi algorithm. Calculates the Viterbi path (the optimal path trough
     the model) and the Viterbi probability to a given model and a given 
     sequence.
     @return Viterbi path
     @param sequences  sequences structure
     @param index      index of sequence to take
     @param log_p      probability of the sequence in the Viterbi path 
                       (return value).
  */
  GHMM_IntVector* viterbi(GHMM_Sequences* sequences, int index, double* log_p = NULL) const;


  /** C Model. */
  sdmodel* c_model;

  /** Alphabet of model. */
  GHMM_Alphabet* alphabet;

  int no_classes;

protected:

  /** */
  bool own_alphabet;

  void init();

  void init(GHMM_Alphabet *alphabet);

  void init(int number_of_states, int alphabet_size, double prior=-1);

  void setNodeTag (const string& tag);
  void setTransitionTag (const string& tag);

  void buildCppData();

  void cleanCPP(); 
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_DISCRETEMODEL_H */
