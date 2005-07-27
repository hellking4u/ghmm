#ifndef PVITERBI_PROPAGATE_H
#define PVITERBI_PROPAGATE_H
#ifdef __cplusplus
extern "C" {
#endif

#include "pmodel.h"
#include "psequence.h"
#include "linkedlist.h"

/*------------        Here comes the Propagate stuff          ------------- */

#define PROP_EPS 10e-6

/** Matrix cell **/
typedef struct cell {
  /** x coordinate of the cell **/
  int x;
  /** y coordinate of the cell **/
  int y;
  /** state coordinate of the cell **/
  int state;
  /** state that was inhabited before **/
  int previous_state;
  /** probability up to the previous cell **/
  double log_p;
  /** transition from the provious state to this state **/
  double log_a;
  int ref_count;
} cell;

typedef struct plocal_propagate_store_t {
  /** precomputed log probabilities for transitions into the states 
      for each transition class of the source state **/
  double *** log_in_a;
  /** precomputed log probabilities for each state for the emissions **/
  double **log_b;
  /** lookback matrix for the last offset steps **/
  double *** phi;
  /** log probabilities for the current u, v and every state **/
  double *phi_new;
  /** traceback matrix of cell pointers **/
  cell **** end_of_first;
  /** for convinience store a pointer to the model **/
  pmodel * mo;
  /** for the debug mode store information of matrix sizes **/
  /** length of sequence X determines size of psi **/
  int len_x;
  /** length of sequence Y determines size of phi and psi **/
  int len_y;
  /** for the recursion reuse memory here is the start index **/
  int start_x;
  int start_y;
  /** non functional stuff **/
  int    *topo_order;
  int    topo_order_length;
} plocal_propagate_store_t;

cell * init_cell(int x, int y, int state, int previous_state, double log_p, double log_a);

int * pviterbi_propagate(pmodel *mo, mysequence * X, mysequence * Y, double *log_p, int *path_length, double max_size);

int * pviterbi_propagate_recursion(pmodel *mo, mysequence * X, mysequence * Y, double *log_p, int *path_length, cell *start, cell *stop, double max_size, plocal_propagate_store_t * pv);

cell * pviterbi_propagate_step(pmodel *mo, mysequence * X, mysequence * Y, cell * start, cell * stop, double * log_p, plocal_propagate_store_t * pv);

plocal_propagate_store_t * pviterbi_propagate_alloc(pmodel *mo, int len_y);
  
int pviterbi_propagate_free(plocal_propagate_store_t **v, int n, int max_offset_x, int max_offset_y, int len_y);

#ifdef __cplusplus
}
#endif


#endif
