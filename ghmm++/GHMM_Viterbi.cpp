#include <float.h>
#include <math.h>
#include <assert.h>

#include <ghmm/matrix.h>
#include <ghmm/sdmodel.h>
#include <ghmm++/GHMM_Computation.h>
#include <iostream>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

static double **mat_d_alloc(int n, int m)
{
  int i,j;
  double** matout = new double*[n];
  for (i=0;i<n;i++)
    {
      matout[i] = new double[m];
    }
  return matout;
}

static void mat_d_free(double ***mat, int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      delete [] (*mat)[i];
    }
  delete *mat;
}

static int **mat_i_alloc(int n, int m)
{
  int i,j;
  int** matout = new int*[n];
  for (i=0;i<n;i++)
    {
      matout[i] = new int[m];
    }
  return matout;
}

static void mat_i_free(int ***mat, int n)
{
  int i;
  for (i=0;i<n;i++)
    {
      delete [] (*mat)[i];
    }
  delete *mat;
}

Viterbi::Viterbi( sdmodel *mo, int len)  
/**  
     double ***log_in_a;
     double **log_b;
     double *phi;
     double *phi_new;
     int    **psi;
*/
{
  int i,j,t,k;
  /* Allocate the log_in_a's -> individal lenghts */

  log_in_a = new double**[mo->N];
  for(i = 0; i < mo->N; i++) 
    {
      log_in_a[i] = mat_d_alloc(mo->cos, mo->s[i].in_states);
      assert(log_in_a[i]);
    }
      
  log_b    = mat_d_alloc(mo->N, len);
  psi      = mat_i_alloc(len, mo->N);

  assert(log_in_a);
  assert(log_b);
  assert(psi);

  phi     = new double[mo->N];
  phi_new = new double[mo->N];  
  for(i = 0; i < mo->N; i++) phi[i]     = 0.0;
  for(i = 0; i < mo->N; i++) phi_new[i] = 0.0;
}


Viterbi::~Viterbi()
{
  delete [] log_in_a;
  delete [] log_b;
  delete psi;
  delete [] phi;
  delete [] phi_new;
}

void Viterbi::DFSVisit(sdmodel *mo, int j, int &counter)
{
  colors[j] = VISITED; 
  for(int i = 0; i < mo->s[j].out_states; i++)	
    {
      if (colors[mo->s[j].out_id[i]] == NOTVISITED) /* looping back taken care of */
	DFSVisit(mo,mo->s[j].out_id[i],counter);
    }
  colors[j]   = DONE;
  doneTime[j] = ++counter;
  if ( mo->silent[j] ) topo_order.push_back(j);
}

void Viterbi::DFS(sdmodel *mo)
{
  int i;
  colors   = new DFSFLAG[mo->N];
  doneTime = new int[mo->N];

  for(i=0; i < mo->N; i++)
    {
      colors[i] = NOTVISITED;
    }
  int counter = 0;
  for(i=0; i < mo->N; i++)
  {
    if (colors[i] == NOTVISITED) DFSVisit(mo, i, counter);
  }

  vector<int>::iterator it;
  for( it = topo_order.end() - 1; it != topo_order.begin(); it-- )
    cout << *it << " , ";
  cout << endl;
}

void Viterbi::Viterbi_precompute( sdmodel *mo, int *o, int len )
{
  int i, j, k, t, osc;
  double log_p = +1;

  /* Precomputing the log(a_ij) */

  //for (j = 0; j < mo->N; j++)
  // for (i = 0; i < mo->s[j].in_states; i++)
  //  if ( mo->s[j].in_a[i] == 0.0 )   /* DBL_EPSILON ? */
  //log_in_a[j][i] = +1; /* Not used any further in the calculations */
  //  else
  //log_in_a[j][i] = log( mo->s[j].in_a[i] );

  /* Collecting emitting states */
  for(j = 0; j < mo->N; j++) 
    {
      if (!mo->silent[j]) stateind.push_back(j);
    }

  for(j = 0; j < mo->N; j++) 
    {
      for(k = 0; k < mo->cos; k++)
	for(i = 0; i < mo->s[j].in_states; i++)	
	  if ( mo->s[j].in_a[k][i] == 0.0 )   /* DBL_EPSILON ? */
	    log_in_a[j][k][i] = +1; /* Not used any further in the calculations */
	  else 
	    log_in_a[j][k][i] = log( mo->s[j].in_a[k][i] );	
    }


  /* Precomputing the log(bj(ot)) */
  for (j = 0; j < mo->N; j++)
    for (t = 0; t < len; t++)
      if ( mo->s[j].b[o[t]] == 0.0 )   /* DBL_EPSILON ? */ 
	log_b[j][t] = +1; 
      else
	log_b[j][t] = log( mo->s[j].b[o[t]] );
}/* viterbi_precompute */


/** Return the log score of the sequence */
double Viterbi::Viterbi_runme   ( sdmodel *mo, int *o, int len)
{
  int t, j, i;
  int last_osc = -1;
  double value, max_value;  
  double log_p, sum, osum = 0.0;
  double dummy = 0.0;

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  // Done in the constructor Viterbi

  /* Precomputing the log(a_ij) and log(bj(ot)) */
  Viterbi_precompute(mo, o, len);

  // The Viterbi path
  state_seq = new int[len];

  /* Initialization, that is t = 0 */
  for (j = 0; j < mo->N; j++) {
    if ( mo->s[j].pi == 0.0 || log_b[j][0] == +1 ) /* instead of 0, DBL_EPS.? */
      phi[j] = +1;
    else
      phi[j] = log(mo->s[j].pi) + log_b[j][0];

    printf("phi[%d],%f\n", j, phi[j]);
  }
  /* psi[0][i] = 0, also unneccessary here */

  /* Recursion step */
  for (t = 1; t < len; t++) {

    int osc = mo->get_class(&dummy,t,&osum);

    for (j = 0; j < mo->N; j++) {
      /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      max_value = -DBL_MAX;
      psi[t][j] = -1;
      for (i = 0; i < mo->s[j].in_states; i++) {

	if ( last_osc != -1 && 
	     last_osc != osc ) { /* just switch class */
	  
	  if (   phi[ mo->s[j].in_id[i] ] != +1 ) {
	    if ( phi[ mo->s[j].in_id[i] ] > max_value ) {
	      max_value = phi[ mo->s[j].in_id[i] ];
	      psi[t][j] = mo->s[j].in_id[i];			       
	    }
	  }

	} else {
	  if ( phi[ mo->s[j].in_id[i] ] != +1 &&
	       log_in_a[j][osc][i]    != +1) {
	    value = phi[ mo->s[j].in_id[i] ] + log_in_a[j][osc][i];
	    if (value > max_value) {
	      max_value = value;
	      psi[t][j] = mo->s[j].in_id[i];
	    }
	  }
	  else
	    fprintf(stderr, " %d --> %d = %f, \n", i,j,log_in_a[j][osc][i]);
	}
      }

      last_osc = osc; /* save last transition class */

      /* No maximum found (that is, state never reached)
         or the output O[t] = 0.0: */
      if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	  log_b[j][t] == +1 ) {
	phi_new[j] = +1;
	fprintf(stderr, "b[%d][%d]=0.0, \n", j, t);
      }
      else
	phi_new[j] = max_value + log_b[j][t];
    } /* complete time step */

    /* First now replace the old phi with the new phi */
    for (j = 0; j < mo->N; j++) 
      {      
	phi[j] = phi_new[j];
	printf("\npsi[%d],%d, phi, %f\n", t, psi[t][j], phi[j]);
      }
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[len-1] = -1;
  for (j = 0; j < mo->N; j++)
    if ( phi[j] != +1 && phi[j] > max_value) { 
      max_value = phi[j];
      state_seq[len-1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* Sequence can't be generated from the model! */
    log_p = +1;
    /* Backtracing doesn't work, because state_seq[*] allocated with -1 */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = -1;    
  }
  else {
    log_p = max_value;
    /* Backtracing */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = psi[t+1][state_seq[t+1]];
  }

  return log_p;
} /* viterbi */



/** Return the log score of the sequence */
double Viterbi::Viterbi_runme_silent   ( sdmodel *mo, int *o, int len)
{
  int t, k, i, j;
  int last_osc = -1;
  double value, max_value;  
  double log_p, sum, osum = 0.0;
  double dummy = 0.0;

  /* Allocate the matrices log_in_a, log_b,Vektor phi, phi_new, Matrix psi */
  // Done in the constructor Viterbi

  /* Precomputing the log(a_ij) and log(bj(ot)) */
  Viterbi_precompute(mo, o, len);

  // The Viterbi path
  vector<int> seq_local;
  state_seq = new int[len];

  /* Initialization for all states (emitting and non-silent), that is t = 0 */
  for (j = 0; j < mo->N; j++) {
    if ( mo->s[j].pi == 0.0 || log_b[j][0] == +1 ) /* instead of 0, DBL_EPS.? */
      phi[j] = +1;
    else
      phi[j] = log(mo->s[j].pi) + log_b[j][0];

    printf("phi[%d],%f\n", j, phi[j]);
  }
  /* psi[0][i] = 0, also unneccessary here */

  
  /* Recursion step */
  for (t = 1; t < len; t++) {

    int osc = mo->get_class(&dummy,t,&osum);
    int St;
    for (k = 0; k < stateind.size(); k++) {

      /* Determine the maximum */
      /* max_phi = phi[i] + log_in_a[j][i] ... */
      St        = stateind[k];
      max_value = -DBL_MAX;
      psi[t][St] = -1;
      for (i = 0; i < mo->s[St].in_states; i++) {

	if ( last_osc != -1 && 
	     last_osc != osc ) { /* just switch class */
	  
	  if (   phi[ mo->s[St].in_id[i] ] != +1 ) {
	    if ( phi[ mo->s[St].in_id[i] ] > max_value ) {
	      max_value = phi[ mo->s[St].in_id[i] ];
	      psi[t][St] = mo->s[St].in_id[i];			       
	    }
	  }

	} else {
	  if ( phi[ mo->s[St].in_id[i] ] != +1 &&
	         log_in_a[St][osc][i]    != +1) {
	    value = phi[ mo->s[St].in_id[i] ] + log_in_a[St][osc][i];
	    if (value > max_value) {
	      max_value = value;
	      psi[St][j] = mo->s[St].in_id[i];
	    }
	  }
	  else
	    fprintf(stderr, " %d --> %d = %f, \n", i,St,log_in_a[St][osc][i]);
	}
      }

      last_osc = osc; /* save the last transition class */

      /* No maximum found (that is, state never reached)
         or the output O[t] = 0.0: */
      if (max_value == -DBL_MAX || /* and then also: (v->psi[t][j] == -1) */
	  log_b[St][t] == +1 ) {
	phi_new[St] = +1;
	fprintf(stderr, "b[%d][%d]=0.0, \n", j, t);
      }
      else
	phi_new[St] = max_value + log_b[St][t];
    } /* complete time step for emitting states */

    // 
    // Silent states, visited in topological order
    //
    
    /* first compute only in-coming from emitting states, then non-silent */
    for ( k = 0; k < mo->N; k++) 
      {
	if ( mo->silent[k] ) /* Silent states */ 
	  {
	    /* Determine the maximum */
	    /* max_phi = phi[i] + log_in_a[j][i] ... */
	    max_value = -DBL_MAX;
	    psi[t][k] = -1;
	    for (i = 0; i < mo->s[k].in_states; i++) {	      
	      if ( phi[ mo->s[k].in_id[i] ] != +1 &&
		     log_in_a[k][osc][i]    != +1) {
		value = phi[ mo->s[k].in_id[i] ] + log_in_a[k][osc][i];
		if (value > max_value) {
		  max_value = value;
		  psi[k][j] = mo->s[k].in_id[i];
		}
	      }
	    }

	    /* No maximum found (that is, state never reached)
	       or the output O[t] = 0.0: */
	    if (max_value    == -DBL_MAX)
	      {
		phi_new[k] = +1;
		fprintf(stderr, "b[%d][%d]=0.0, \n", j, t);
	      }
	    else
	      phi_new[k] = max_value;
	  }
      } /* complete time step for silent states */
    
    /* First now replace the old phi with the new phi */
    for (k = 0; k < mo->N; k++) 
      {      
	phi[k] = phi_new[k];
	printf("\npsi[%d],%d, phi, %f\n", t, psi[t][k], phi[k]);
      }
  }

  /* Termination */
  max_value = -DBL_MAX;
  state_seq[len-1] = -1;
  for (j = 0; j < mo->N; j++)
    if ( phi[j] != +1 && phi[j] > max_value) { 
      max_value = phi[j];
      state_seq[len-1] = j;
    }
  if (max_value == -DBL_MAX) {
    /* Sequence can't be generated from the model! */
    log_p = +1;
    /* Backtracing doesn't work, because state_seq[*] allocated with -1 */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = -1;    
  }
  else {
    log_p = max_value;
    /* Backtracing */
    for (t = len - 2; t >= 0; t--)
      state_seq[t] = psi[t+1][state_seq[t+1]];
  }

  return log_p;
} /* viterbi */

void Viterbi::print_path(int T, char *ts)
{
  printf("%s{ ", ts);
  for(int i=0; i < T; i++)
    {
      printf(" %d, ", state_seq[i] );
    }
  printf("%s} ", ts);
}
