%{
#include <ghmm/model.h>
#include <ghmm/foba.h>
#include <ghmm/viterbi.h>
#include <ghmm/reestimate.h>
#include <ghmm/gradescent.h>
#include <ghmm/kbest.h>
%}

/*==========================================================================
  ===== discrete background distribution =================================== */
typedef struct {
  /** Number of distributions */
  int n;
  /** Number of symbols in alphabet */
  int m;
  /** Order of the respective distribution */
  int *order;
  /** The probabilities */
  double **b;
  /** string ids of the background distributions */
  char * * name;
} ghmm_dbackground;

ghmm_dbackground* ghmm_dbackground_alloc (int n, int m, int* orders, double** B);

%extend ghmm_dbackground {
        ghmm_dbackground(int n, int m, int* orders, double** B)
            { return ghmm_dbackground_alloc(n, m, orders, B); }
        ~ghmm_dbackground() { ghmm_dbackground_free(self); }

        int getOrder(size_t index) { return self->order[index]; }
}


/*==========================================================================
  ===== discrete emission states =========================================== */
typedef struct {
  /** Initial probability */
  double pi;
  /** Output probability */
  double *b;

  /** IDs of the following states */
  int *out_id;
  /** IDs of the previous states */
  int *in_id;

  /** transition probabilities to successor states. */
  double *out_a;
  /** transition probabilities from predecessor states. */
  double *in_a;

  /** Number of successor states */
  int out_states;
  /** Number of precursor states */
  int in_states;
  /** if fix == 1 --> b stays fix during the training */
  int fix;
  /** contains a description of the state (null terminated utf-8)*/
  unsigned char * desc;
  /** x coordinate position for graph representation plotting **/
  int xPosition;
  /** y coordinate position for graph representation plotting **/
  int yPosition;
} ghmm_dstate;

%extend ghmm_dstate {
        void clean();

        int getInState(size_t index) { return self->in_id[index]; }

        int getOutState(size_t index) { return self->out_id[index]; }
}

STRUCT_ARRAY(ghmm_dstate, dstate)
REFERENCE_ARRAY(ghmm_dstate, dstate_ptr)


/*==========================================================================
  ===== discrete emission models =========================================== */
typedef struct {
  /** Number of states */
  int N;
  /** Number of outputs */
  int M;
  /** Vector of the states */
  ghmm_dstate *s;
  /** The a priori probability for the model.
      A value of -1 indicates that no prior is defined. 
      Note: this is not to be confused with priors on emission
      distributions*/
  double prior;

  /* contains a arbitrary name for the model (null terminated utf-8) */
  /*unsigned */char * name;

  /** Contains bit flags for varios model extensions such as
      kSilentStates, kTiedEmissions (see ghmm.h for a complete list)
  */
  int model_type;

  /** Flag variables for each state indicating whether it is emitting
      or not. 
      Note: silent != NULL iff (model_type & kSilentStates) == 1  */
  int *silent;
  /*AS*/
  /** Int variable for the maximum level of higher order emissions */
  int maxorder;
  /** saves the history of emissions as int, 
      the nth-last emission is (emission_history * |alphabet|^n+1) % |alphabet|
      see ...*/
  int emission_history;
  
  /** Flag variables for each state indicating whether the states emissions
      are tied to another state. Groups of tied states are represented
      by their tie group leader (the lowest numbered member of the group).
      
      tied_to[s] == kUntied  : s is not a tied state
      
      tied_to[s] == s        : s is a tie group leader
      
      tied_to[t] == s        : t is tied to state s (t>s)
      
      Note: tied_to != NULL iff (model_type & kTiedEmissions) != 0  */
  int *tied_to;

  /** Note: State store order information of the emissions.
      Classic HMMS have emission order 0, that is the emission probability
      is conditioned only on the state emitting the symbol.

      For higher order emissions, the emission are conditioned on the state s
      as well as the previous emission_order[s] observed symbols.

      The emissions are stored in the state's usual double* b. The order is
      set order.

      Note: order != NULL iff (model_type & kHigherOrderEmissions) != 0  */
  int * order;

  /** ghmm_dbackground is a pointer to a
      ghmm_dbackground structure, which holds (essentially) an
      array of background distributions (which are just vectors of floating
      point numbers like state.b).
      
      For each state the array background_id indicates which of the background
      distributions to use in parameter estimation. A value of kNoBackgroundDistribution
      indicates that none should be used.
      
      Note: background_id != NULL iff (model_type & kHasBackgroundDistributions) != 0  */
  int *background_id;
  ghmm_dbackground *bp;

  /** (WR) added these variables for topological ordering of silent states 
      Condition: topo_order != NULL iff (model_type & kSilentStates) != 0  */
  int *topo_order;
  int topo_order_length;
  
  /** pow_lookup is a array of precomputed powers
      It contains in the i-th entry M (alphabet size) to the power of i
      The last entry is maxorder+1 */
  int *pow_lookup;

  /** Store for each state a class label. Limits the possibly state sequence

      Note: label != NULL iff (model_type & kLabeledStates) != 0  */
  int * label;
  ghmm_alphabet * labelAlphabet;

  ghmm_alphabet * alphabet;
} ghmm_dmodel;

extern int ghmm_dmodel_free(ghmm_dmodel **mo);

%newobject ghmm_dmodel::generate_sequences;
%newobject ghmm_dmodel_generate_sequences;

%extend ghmm_dmodel {
        ghmm_dmodel() { return calloc(1, sizeof(ghmm_dmodel)); }
        ~ghmm_dmodel() { ghmm_dmodel_free(&self); }

        int forward_init(double *alpha_1, int symb, double *scale);

        int forward(const int *O, int len, double **alpha, double *scale,
                    double *log_p);

        int backward(const int *O, int len, double **beta, const double *scale);

        int backward_termination(const int *O, int length, double **beta,
                                 double *scale, double *log_p);

        int logp(const int *O, int len, double *log_p);

        int forward_lean(const int *O, int len, double *log_p);

        //int check_compatibility(ghmm_dmodel **mo, int model_number);

        ghmm_dseq* generate_sequences(int seed, int global_len, long seq_number,
                                      int Tmax);

        double likelihood(ghmm_dseq *sq);

        double get_transition(int i, int j);
        void   set_transition(int i, int j, double prob);

        double prob_distance(ghmm_dmodel *m, int maxT, int symmetric, int verbose);

        int normalize();

        int add_noise(double level, int seed);

        int duration_apply(int cur, int times);

        int background_apply(double *background_weight);

        int get_uniform_background(ghmm_dseq *sq);

        void order_topological();

        void update_tie_groups();

        int baum_welch(ghmm_dseq *sq);

        int baum_welch_nstep(ghmm_dseq *sq, int max_step, double likelihood_delta);

        int* viterbi(int *o, int len, double *log_p);

        double viterbi_logp(int *o, int len, int *state_seq);

        ghmm_dstate* getState(size_t index) { return self->s + index; }

        int  getSilent(size_t index) { return self->silent[index]; }
        void setSilent(size_t index, int value) { self->silent[index] = value; }
}

STRUCT_ARRAY(ghmm_dmodel, dmodel)
REFERENCE_ARRAY(ghmm_dmodel, dmodel_ptr)

%newobject ghmm_dmodel_label_generate_sequences;
%newobject ghmm_dmodel::label_generate_sequences;

/* ====== labeled =========================================================== */
%extend ghmm_dmodel{
        int label_forward(const int *O, const int *label, int len, double **alpha, double *scale, double *log_p);
        int label_logp(const int *O, const int *label, int len, double *log_p);
        int label_backward(const int *O, const int *label, int len, double **beta, double *scale, double *log_p);
        //int label_forward_lean(const int *O, const int *label, int len, double *log_p);
        int label_gradient_descent(ghmm_dseq *sq, double eta, int no_steps);
        int* label_kbest(int *o_seq, int seq_len, int k, double *log_p);
        int label_baum_welch(ghmm_dseq *sq);
        int label_baum_welch_nstep(ghmm_dseq *sq, int max_step, double likelihood_delta);
        ghmm_dseq* label_generate_sequences(int seed, int global_len, long seq_number, int Tmax);

        int getLabel(size_t index) { return self->label[index]; }
}

extern double ghmm_dmodel_label_discrim_perf(ghmm_dmodel** mo, ghmm_dseq** sqs, int noC);
extern int ghmm_dmodel_label_discriminative(ghmm_dmodel** mo, ghmm_dseq** sqs, int noC, int max_steps, int gradient);


/* obsolete stuff */
#ifdef GHMM_OBSOLETE
extern ghmm_dmodel** ghmm_dmodel_read(char* filename, int* mo_number);
#endif
