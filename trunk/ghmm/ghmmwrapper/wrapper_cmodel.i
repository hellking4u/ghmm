%{
#include "ghmm/smodel.h"
#include "ghmm/sfoba.h"
#include "ghmm/sviterbi.h"
#include "ghmm/sreestimate.h"
%}
/*==========================================================================
  ===== continous emission density types =================================== */
typedef enum {
  normal,
  normal_right,
  normal_approx,
  normal_left,
  uniform,
  density_number
} ghmm_density_t;

%inline %{
        ghmm_density_t* density_array_alloc(size_t length) { return malloc(length*sizeof(ghmm_density_t)); }
        ghmm_density_t  density_array_getitem(ghmm_density_t* self, size_t index) { return self[index]; }
        void            density_array_setitem(ghmm_density_t* self, size_t index, ghmm_density_t value) { self[index] = value; }
%}


/*==========================================================================
  ===== continous emission states ========================================== */
typedef struct ghmm_cstate {
  /** Number of output densities per state */
  int M;
  /** initial prob. */
  double pi;
  /** IDs of successor states */
  int *out_id;
  /** IDs of predecessor states */
  int *in_id;
  /** transition probs to successor states. It is a
      matrix in case of mult. transition matrices (COS > 1)*/
  double **out_a;
  /** transition probs from predecessor states. It is a
      matrix in case of mult. transition matrices (COS > 1) */
  double **in_a;
  /** number of  successor states */
  int out_states;
  /** number of  predecessor states */
  int in_states;
  /** weight vector for output function components */
  double *c;
  /** mean vector for output functions (normal density and truncated normal
      density) or max value for uniform distribution */
  double *mue;
  /** variance vector for output functions or min value for uniform distribution */
  double *u;
  /** value where the normal is truncated (only for truncated distributions) */
  double *a;
  /** flag for fixation of parameter. If fix = 1 do not change parameters of
      output functions, if fix = 0 do normal training. Default is 0. */
  int fix;
  /** Flag for density function for each component of the mixture
      0: normal density, 1: truncated normal (right side) 
      density, 2: approximated normal density, 3: truncated normal (left side)
      4: uniform distribution */
  ghmm_density_t *density;
  /**  array of flags for fixing mixture components in the reestimation
        mixture_fix[i] = 1 means mu and sigma of component i are fixed.  **/
  int *mixture_fix;
  /** contains a description of the state (null terminated utf-8)*/
  char* desc;
  /** x coordinate position for graph representation plotting **/
  int xPosition;
  /** y coordinate position for graph representation plotting **/
  int yPosition;
} ghmm_cstate;

%extend ghmm_cstate {
        int alloc(int M, int in_states, int out_states, int cos);

        void  setDensity(size_t i, ghmm_density_t type) { self->density[i] = type; }
        void   setWeight(size_t i, double value)        { self->c[i]       = value; }
        void     setMean(size_t i, double value)        { self->mue[i]     = value; }
        void   setStdDev(size_t i, double value)        { self->u[i]       = value; }
        void setTruncate(size_t i, double value)        { self->a[i]       = value; }

        ghmm_density_t getDensity(size_t i) { return self->density[i]; }
        double          getWeight(size_t i) { return self->c[i]; }
        double            getMean(size_t i) { return self->mue[i]; }
        double          getStdDev(size_t i) { return self->u[i]; }
        double        getTruncate(size_t i) { return self->a[i]; }


        int getInState(size_t index) { return self->in_id[index]; }
        int getOutState(size_t index) { return self->out_id[index]; }

        double getInProb(size_t index)            { return self->in_a[0][index]; }
        double getInProb(size_t index, size_t c)  { return self->in_a[c][index]; }
        double getOutProb(size_t index)           { return self->out_a[0][index]; }
        double getOutProb(size_t index, size_t c) { return self->out_a[c][index]; }
}

STRUCT_ARRAY(ghmm_cstate, cstate)
REFERENCE_ARRAY(ghmm_cstate, cstate_ptr)


/*==========================================================================
  ===== continous emission class change context ============================ */
struct ghmm_cmodel;
typedef struct ghmm_cmodel_class_change_context {

  /* Names of class change module/function (for python callback) */
  char *python_module;
  char *python_function;

  /* index of current sequence */
  int k;

  /** pointer to class function */
  int (*get_class) (struct ghmm_cmodel *, double *, int, int);

  /* space for any data necessary for class switch, USER is RESPONSIBLE */
  void *user_data;
  } ghmm_cmodel_class_change_context;


/*==========================================================================
  ===== continous emission models ========================================== */
  typedef struct {
  /** pointer of continuous model*/
    ghmm_cmodel *smo;
  /** sequence\_d\_t pointer */
    ghmm_cseq *sqd;
  /** calculated log likelihood */
    double *logp;
  /** leave reestimation loop if diff. between successive logp values 
      is smaller than eps */
    double eps;
  /** max. no of iterations */
    int max_iter;
  } ghmm_cmodel_baum_welch_context;

%extend ghmm_cmodel_baum_welch_context{
        ghmm_cmodel_baum_welch_context(ghmm_cmodel *smo, ghmm_cseq *sqd)
        {
                ghmm_cmodel_baum_welch_context *bwc = malloc(sizeof(ghmm_cmodel_baum_welch_context));
                bwc->smo = smo;
                bwc->sqd = sqd;
                bwc->logp = malloc(sizeof(*bwc->logp));
                return bwc;
        }
        ~ghmm_cmodel_baum_welch_context() {
                free(self->logp);
                free(self);
        }
}


/*==========================================================================
  ===== continous emission models ========================================== */
%newobject ghmm_cmodel::generate_sequences;

typedef struct ghmm_cmodel {
  /** Number of states */
  int N;
  /** Maximun number of components in the states */
  int M;
  /** ghmm_cmodel includes continuous model with one transition matrix 
      (cos  is set to 1) and an extension for models with several matrices
      (cos is set to a positive integer value > 1).*/
  int cos;
  /** prior for a priori prob. of the model. -1 means no prior specified (all
      models have equal prob. a priori. */
  double prior;

  /* contains a arbitrary name for the model (null terminated utf-8) */
  char * name;

  /** Contains bit flags for varios model extensions such as
      kSilentStates (see ghmm.h for a complete list) */
  int model_type;

  /** All states of the model. Transition probs are part of the states. */
  ghmm_cstate *s;

  /* pointer to a ghmm_cmodel_class_change_context struct necessary for multiple transition
     classes */
  ghmm_cmodel_class_change_context *class_change;

} ghmm_cmodel;

extern int ghmm_cmodel_free(ghmm_cmodel **smo);

%extend ghmm_cmodel {
        ghmm_cmodel() { return calloc(1, sizeof(ghmm_cmodel)); }
%apply SWIGTYPE* DISOWN {ghmm_cmodel *smo};
        ghmm_cmodel(ghmm_cmodel *smo) { return smo; }
%clear ghmm_cmodel *smo;
        ghmm_cmodel(int no_states, int cos) {
                ghmm_cmodel *mo = calloc(1, sizeof(ghmm_cmodel));
                mo->model_type = kContinuousHMM;
                mo->N = no_states;
                mo->M = 1;
                mo->cos = cos;
                mo->prior = -1;
                if (cos > 1)
                    ghmm_cmodel_class_change_alloc(mo);
                mo->s = calloc(mo->N, sizeof(ghmm_cstate));
                return mo;
        }
        ~ghmm_cmodel() { ghmm_cmodel_free(&self); }

        int write_xml(char* filename) { return ghmm_cmodel_xml_write(&self, filename, 1); }

        int forward(double *O, int T, double ***b, double **alpha, double *scale, double *log_p);

        int backward(double *O, int T, double ***b, double **beta, const double *scale);

        int logp(double *O, int T, double *log_p);

        int logp_joint(const double *O, int len, const int *S, int slen, double *log_p);

        int class_change_alloc(void);

        //int check_compatibility(ghmm_cmodel **smo, int smodel_number);

        double get_random_var(int state, int m);

        ghmm_cseq* generate_sequences(int seed, int global_len, long seq_number, int Tmax);

        int likelihood(ghmm_cseq *sqd, double *log_p);

        int individual_likelihoods(ghmm_cseq *sqd, double *log_ps);

        double calc_cmbm(int state, int m, double omega);

        double calc_b(int state, double omega);

        double prob_distance(ghmm_cmodel *cm, int maxT, int symmetric, int verbose);

        double calc_cmBm(int state, int m, double omega);

        double calc_B(int state, double omega);

        //int count_free_parameter(ghmm_cmodel **smo, int smo_number);

        void get_interval_B(int state, double *a, double *b);

        int normalize();

        double get_transition(int i, int j, int c);
        int    check_transition(int i, int j, int c);
        void   set_transition(int i, int j, int c, double prob);

        int* viterbi(double *O, int T, double *log_p);

        ghmm_cstate* getState(size_t index) { return self->s + index; }

        void setModelTypeFlag(unsigned int flag) { return; }
}

extern int ghmm_cmodel_xml_write(ghmm_cmodel** smo, const char* file, int smo_number);
extern int ghmm_cmodel_baum_welch(ghmm_cmodel_baum_welch_context* cs);

STRUCT_ARRAY(ghmm_cmodel, cmodel)
REFERENCE_ARRAY(ghmm_cmodel, cmodel_ptr)


/* obsolete stuff */
#ifdef GHMM_OBSOLETE

// ignore the input value for ghmm_cmodel array to python list conversion
%typemap(in, numinputs=0) int *smo_number (int temp) {
    $1 = &temp;
}
// convert array of ghmm_cmodels to python list
%typemap(argout) (int *smo_number) {
    int i;
    PyObject *obj;
    Py_XDECREF($result);   /* Blow away any previous result */
    if (result) {
        $result = PyList_New(*$1);
        for (i=0; i<*$1; i++) {
            obj = SWIG_NewPointerObj(result[i], SWIGTYPE_p_ghmm_cmodel, SWIG_POINTER_NEW);
            PyList_SetItem($result, i, obj);
        }
        free(result);
    }
    else {
        PyErr_SetString(PyExc_ValueError,"got a null pointer");
        return NULL;
    }
}
extern ghmm_cmodel **ghmm_cmodel_read(const char *filename, int *smo_number);

#endif /* GHMM_OBSOLETE */
