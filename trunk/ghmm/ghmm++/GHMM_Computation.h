/*
  created: 
  authors: 
  file   : $Source$
  $Id$

  __copyright__
*/


#include <ghmm/sdmodel.h>
#include <ghmm++/begin_code.h>

#include <vector>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class Viterbi {
private:  
  enum DFSFLAG { DONE, NOTVISITED, VISITED };

  double ***log_in_a;
  double **log_b;
  double *phi;
  double *phi_new;
  int    **psi;
  vector<int>  stateind;

  DFSFLAG *colors;
  int     *doneTime;
  vector<int>  topo_order;

  
public:
  int    *state_seq;

  Viterbi() 
  {
    log_in_a = NULL;
    log_b    = NULL;
    psi      = NULL;
  }

  ~Viterbi();

  Viterbi( sdmodel *mo, int len);

  void Viterbi_precompute( sdmodel *mo, int *o, int len);

  /** Return the log score of the sequence */
  double Viterbi_runme          ( sdmodel *mo, int *o, int len);
  double Viterbi_runme_silent   ( sdmodel *mo, int *o, int len);

  void   print_path      (int len, char *ts);

  void  DFSVisit( sdmodel *mo, int j, int &counter);
  void  DFS( sdmodel *mo );

};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

/* _GHMM_VITERBI_H */
