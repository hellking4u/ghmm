/*
  created: 4. April 2003 by Wasinee Rungsarityotin
  authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
  file   : $Source$
  $Id$

  __copyright__
*/

#ifndef _GHMM_DISCRETEMODELT_H
#define _GHMM_DISCRETEMODELT_H 1

#include <vector>
#include <ghmm/sdmodel.h>
#include <ghmm++/GHMM_Alphabet.h>

#include <ghmm++/GHMM_AbstractModelT.hh> // Template

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
class GHMM_DiscreteModelT : public GHMM_AbstractModelT<GHMM_GMLState, GHMM_GMLTransition> {
 public:
  
  /** Constructor. */
  GHMM_DiscreteModelT ()
    {
      init();
    }

  /** Constructor. 
      @param number_of_states Number of the states.
      @param alphabet_size Size of the alphabet. 
      @param prior Prior for the a priori probability for the model 
                   (-1 for none). */

  GHMM_DiscreteModelT(int number_of_states, int alphabet_size, double prior=-1)
  /** Constructor. Construct from c model. Object now is owner of this model. 
      @param my_model model as C data structure. */
    {
    }

  /** Constructor. */
  GHMM_DiscreteModelT( sdmodel* my_model)
    {
      init();
      c_model = my_model;
      buildCppData();
    }


  /** Constructor. */
  GHMM_DiscreteModelT(GHMM_Alphabet* my_alphabet)
    {
      init(0,my_alphabet->size());
      alphabet = my_alphabet;
    }

  /** Destructor. */

  virtual ~GHMM_DiscreteModelT() {;}

  /** Returns name of class. */
  virtual const char* toString() const {;}

  /** Adds state with given id to model. */
  void addState(const string& my_id);

  /** 
      Writes transition matrix of a model.
      @param file: output file
      @param tab:  format: leading tabs
      @param separator: format: seperator for columns
      @param ending:    format: end of a row 
  */
  virtual void A_print(FILE *file, char *tab, char *separator, char *ending)
    { ;  }

  /**
     Writes output matrix of a model.
     @param file: output file
     @param tab:  format: leading tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  virtual void B_print(FILE *file, char *tab, char *separator, char *ending)
    { ;  }

  /**
     Tests if all standardization requirements of model are fulfilled. 
     (That is, if the sum of the probabilities is 1).
     @return 0 on succes; -1 on error. 
  */
  virtual int check() { return 0; }

  /**
     Copies a given model. Allocates the necessary memory.
     @return copy of the model
  */
  virtual GHMM_DiscreteModelT* copy() = 0;


  /** 
      Backward-Algorithm. 
      Calculates beta[t][i] given an integer sequence and a model. Scale factors 
      given as parameter (come from sfoba\_forward).
      @param seq:     sequences
      @param index:   index of sequence to take
      @param beta     beta[t][i]
      @param scale    scale factors
      @return 0 for success, -1 for error
  */
  //int foba_backward(GHMM_Sequences* seq, int index, double **beta, const double *scale) const;

  
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
     Calculation of  log( P(O|lambda) ). 
     Done by calling sfoba\_forward. Use this function if only the
     log likelihood and not alpha[t][i] is needed.
     @param seq:      sequences
     @param index:    index of sequence to take
     @param log\_p    log likelihood log( P(O|lambda) )
     @return 0 for success, -1 for error
  */
  // int foba_logp(GHMM_Sequences* seq, int index, double *log_p) const;

  /**
     Writes fix vector of a matrix.
     @param file: output file
     @param tab:  format: leading Tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  // void fix_print(FILE *file, char *tab, char *separator, char *ending) const;
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
  */
  virtual GHMM_Sequences* generate_sequences(int seed, int global_len, long seq_number)
    {
      return NULL; 
    }

  /** Returns alphabet of model or NULL, if no such alphabet exists. */
  virtual GHMM_Alphabet* getAlphabet();

  /** Returns model type. */
  GHMM_ModelType getModelType() { return GHMM_DISCRETE; }

  /** */
  virtual int getNumberOfTransitionMatrices();

  /* Returns state with given index. */
  sdstate* getCState(int index) const
    {
      if (index >= c_model->N) {
	fprintf(stderr,"GHMM_SWDiscreteModel::getCState(int):\n");
	fprintf(stderr,"State no. %d does not exist. Model has %d states.\n",index,c_model->N);
	exit(1);
      }
      
      return &c_model->s[index]; 
    }

  /**
     Writes initial allocation vector of a matrix.
     @param file: output file
     @param tab:  format: leading Tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  virtual void Pi_print(FILE *file, char *tab, char *separator, char *ending);
  
  /**
     Writes transposed initial allocation vector of a matrix.
     @param file: output file
     @param tab:  format: leading Tabs
     @param separator: format: seperator for columns
     @param ending:    format: end of a row  
  */
  //void Pi_print_transp(FILE *file, char *tab, char *ending);
  /**
     Writes the model in matrix format.
     @param file: output file
  */
  //virtual void print(FILE *file) const;
  /** Computes probabilistic distance of this model to a second model. 
      This model is used to generate random output. The second model
      is compared to these sequences.
      @return probabilistic distance
      @param m  model to compare with
      @param maxT  maximum output length (for HMMs with absorbing states multiple
                   sequences with a toal langth of at least maxT will be 
                   generated)
      @param symmetric  flag, whether to symmetrize distance (not implemented yet)
      @param verbose  flag, whether to monitor distance in 40 steps. 
                      Prints to stdout (yuk!)
  */
  //double prob_distance(GHMM_DiscreteModel* m, int maxT, int symmetric, int verbose);
  /** 
      Writes the parameters of this model sorted by states. 
      Is not very concise.   
      @param file: output file
  */
  //void states_print(FILE *file);

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


  /**
     Calculates the logarithmic probability to a given path through the 
     states (does not have to be the Viterbi path), given sequence and
     a model.
     @param seq:       sequence
     @param index:     index of sequence to take
     @param state_seq: path through the states
     @return log P
  */
  // double viterbi_logp(GHMM_Sequences* seq, int index, int* state_seq);
  /**
     Baum-Welch Algorithm for HMMs.
     Training of model parameter with multiple integer sequences (incl. scaling).
     New parameters set directly in hmm (no storage of previous values!).
     @return           0/-1 success/error
     @param seq        sequences used for training.
  */
  //int reestimate_baum_welch(GHMM_Sequences* seq);

  /** C Model. */
  sdmodel* c_model;

  /** Alphabet of model. */
  GHMM_Alphabet* alphabet;

 protected:
  virtual void setNodeTag      (const string& tag);
  virtual void setTransitionTag(const string& tag);


 protected:

  /** */
  bool own_alphabet;

  /** Build c++ data from c_model. */
  virtual void buildCppData() {;}

  /** */
  void cleanCPP();

  /** Init function. */
  virtual void init()
    { 
      //attributes["type"] = "discrete";
      alphabet           = NULL;
      c_model            = NULL;
      own_alphabet       = false;      
    }

  /** Init function.
      @param number_of_states Number of the states.
      @param alphabet_size Size of the alphabet. 
      @param prior Prior for the a priori probability for the model 
                   (-1 for none). */
  virtual void init(int number_of_states, int alphabet_size, double prior=-1)
    {
      int i;
      int j;

      init();
      c_model = (sdmodel*) calloc(1,sizeof(sdmodel));
      if (!c_model) {
		fprintf(stderr,"GHMM_DiscreteModel::GHMM_DiscreteModel() could not allocate c_model\n");
		exit(1);
      }

      c_model->N       = number_of_states;
      c_model->M       = alphabet_size;
      c_model->prior   = prior;
      c_model->s       = (sdstate*) malloc(sizeof(sdstate) * max(c_model->N,1));
      /* initialize all states. */
      
      for (i = 0; i < number_of_states; ++i) {
	  c_model->s[i].pi         = 0;
	  c_model->s[i].b          = (double*) malloc(sizeof(double) * alphabet_size);
	/* output probabilities are initialized with 0. */
	for (j = 0; j < alphabet_size; ++j)
	  c_model->s[i].b[j] = 0;
	c_model->s[i].out_id     = NULL;
	c_model->s[i].in_id      = NULL;
	c_model->s[i].out_a      = NULL;
	c_model->s[i].in_a       = NULL;
	c_model->s[i].out_states = 0;
	c_model->s[i].in_states  = 0;
	c_model->s[i].fix        = 0;
      }

      buildCppData();
    }

  /** */
  virtual void init(GHMM_Alphabet *my_alphabet)
    {
      attributes["type"] = "discrete";
      alphabet           = my_alphabet;
      c_model            = NULL;
      own_alphabet       = true; 
      
      c_model = (sdmodel*) calloc(1,sizeof(sdmodel));
      if (!c_model) {
	fprintf(stderr,"GHMM_DiscreteModel::GHMM_DiscreteModel() could not allocate c_model\n");
	exit(1);
      }

      c_model->N        = 0;
      c_model->M        = alphabet->size();
      buildCppData();
    }

  /** */
  virtual void XMLIO_finishedReading();
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& tag, XMLIO_Attributes &attrs);

  virtual void createTransitions();

  virtual const int      XMLIO_writeContent(XMLIO_Document& writer); 
};

#ifdef HAVE_NAMESPACES
}
#endif


#endif /* _GHMM_DISCRETEMODEL_H */
