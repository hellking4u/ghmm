/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/cluster.c
  created      : TIME: 10:55:33     DATE: Tue 02. June 1998
  $Id$

__copyright__

*******************************************************************************/

#include <stdio.h>
#include "mprintf.h"
#include "mes.h"
#include "cluster.h"
#include "rng.h"
#include "model.h"
#include "sequence.h"
#include "reestimate.h"
#include "foba.h"
#include "const.h"
#include "matrix.h"

/*============================================================================*/
int cluster_hmm(char *seq_file, char *mo_file, char *out_filename)  {
# define CUR_PROC "cluster_hmm"
  int res = -1, i, iter = 0, sq_number;
  sequence_t *sq = NULL, **sq_vec = NULL;
  long j, changes = 1; 
  long *oldlabel;
  double log_p;
  FILE *outfile = NULL;
  cluster_t cl;
  cl.mo = NULL;
  cl.mo_seq = NULL;  
  if(!(outfile = mes_fopen(out_filename, "wt"))) {mes_proc(); goto STOP;}

/*----------------------------------------------------------------------------*/

  /* Allocate memory and read in data: Sequences and initial model */
  sq_vec = sequence_read(seq_file, &sq_number);
  if (!sq_vec[0]) {mes_proc(); goto STOP;}
  if (sq_number > 1) 
    mes_prot("Warning: seq. file contains multiple seq. arrays. \
                      Only first array is used for clustering\n");
  sq = sq_vec[0];
  fprintf(outfile, "Cluster Sequences\n");
  sequence_print(outfile, sq);
  cl.mo = model_read(mo_file, &cl.mo_number);
  if (!cl.mo) {mes_proc(); goto STOP;} 
  if(!m_calloc(oldlabel, sq->seq_number)) {mes_proc(); goto STOP;}
  for (i = 0; i < sq->seq_number; i++)
    oldlabel[i] = (-1);
  if(!m_calloc(cl.mo_seq, cl.mo_number)) {mes_proc(); goto STOP;}
  for (i = 0; i < cl.mo_number; i++)
    cl.mo_seq[i] = NULL;
  if (model_check_compatibility(cl.mo, cl.mo_number)) { 
    mes_proc(); goto STOP; 
  } 

/*----------------------------------------------------------------------------*/
  fprintf(outfile, "\nInitial Models:\n");
  for (i = 0; i < cl.mo_number; i++)
    model_print(outfile, cl.mo[i]);
/*----------------------------------------------------------------------------*/

  while(changes > 0) { 
    iter++;

    /* Associate the sequences */
    fprintf(outfile, "\nSequence, Best Model, logP of generating Seq.:\n");
    for (j = 0; j < sq->seq_number; j++) {
      sq->seq_label[j] = 
	sequence_best_model(cl.mo, cl.mo_number, sq->seq[j], sq->seq_len[j],
			    &log_p);    
      fprintf(outfile, "seq %ld, mo %ld, log p %.4f\n", j, 
	      sq->seq_label[j], log_p);
      if (sq->seq_label[j] == -1 || sq->seq_label[j] >= cl.mo_number) { 
	/* No model fits! */ 
	char *str =
	  mprintf(NULL, 0, "Seq. %ld: sequence_best_model gives %d\n",
		  j, sq->seq_label[j]); 
	mes_prot(str); m_free(str); goto STOP;
      }
    }
    if (cluster_avoid_empty_model(sq->seq_label, sq->seq_number, cl.mo_number))
      { mes_proc(); goto STOP; }
    changes = cluster_update_label(oldlabel, sq->seq_label, sq->seq_number);
    fprintf(outfile, "%ld changes\n", changes);
    fprintf(stdout, "\n*** %ld changes in iteration %d ***\n\n", changes, iter);

    /* Reestimate models with the associated sequences */
    if (changes > 0) {
      if (cluster_update(&cl, sq)) {
	mes_proc(); goto STOP; 
      }
      fprintf(outfile, "\nGes. WS VOR %d.Reestimate:\n", iter);
      cluster_print_likelihood(outfile, &cl);
      for (i = 0; i < cl.mo_number; i++) {
	if (reestimate_baum_welch(cl.mo[i], cl.mo_seq[i])) {
	  char *str =  
	    mprintf(NULL,0,"%d.reestimate false, mo[%d]\n", iter, i); 
	  mes_prot(str);	m_free(str);
	  /* model_print(stdout, cl.mo[i]); */
	  goto STOP; 
	}
      }
      fprintf(outfile, "\nGes. WS NACH %d.Reestimate:\n", iter);
      cluster_print_likelihood(outfile, &cl);
    } /* if changes */
  } /* while */

/*----------------------------------------------------------------------------*/

  if (!cluster_ausgabe(&cl, sq, outfile, out_filename)) {
    mes_proc(); goto STOP;
  }
      
  res = 0;
STOP:
  /* ...noch div. free! */
  if (outfile) fclose(outfile);
  return(res);
# undef CUR_PROC
}/* cluster_hmm */

/*============================================================================*/
int cluster_update(cluster_t *cl, sequence_t *sq) {
#define CUR_PROC "cluster_update"
  int i, res = -1;
  long *seq_counter;
  sequence_t *seq_t;
  if(!m_calloc(seq_counter, cl->mo_number)) {mes_proc(); goto STOP;}
  /* Fix the number of associated sequences */
  for (i = 0; i < sq->seq_number; i++) 
    seq_counter[sq->seq_label[i]]++;
  /* allocate memory block for block */
  for (i = 0; i < cl->mo_number; i++) {
    if (cl->mo_seq[i]) {
      /* Important: No sequence_free here, otherwise are the original 
	 sequences gone! */
      sequence_clean(cl->mo_seq[i]);
      m_free(cl->mo_seq[i]);
    }
    cl->mo_seq[i] = sequence_calloc(seq_counter[i]);
    cl->mo_seq[i]->seq_number = 0; /* Counted later */
  }
  /* Set entries */
  for (i = 0; i < sq->seq_number; i++) {
    seq_t = cl->mo_seq[sq->seq_label[i]];
    seq_t->seq_len[seq_t->seq_number] = sq->seq_len[i];
    seq_t->seq[seq_t->seq_number] = sq->seq[i]; /* Pointer!!! */
    seq_t->seq_label[seq_t->seq_number] = sq->seq_label[i];
    seq_t->seq_number++;
  }
  res = 0;
 STOP:
  m_free(seq_counter);
  return(res);
# undef CUR_PROC
} /* cluster_update */

/*============================================================================*/
void cluster_print_likelihood(FILE *outfile, cluster_t *cl) {
  double ges_prob = 0.0, mo_prob;
  int i;
  for (i = 0; i < cl->mo_number; i++) {
    mo_prob = model_likelihood(cl->mo[i], cl->mo_seq[i]); 
    ges_prob += mo_prob;
    fprintf(outfile, "mo %d (#Seq. %ld): %.4f\n", i, 
	    cl->mo_seq[i]->seq_number, mo_prob);
  }
  fprintf(outfile, "Summe: %.4f\n\n", ges_prob);
} /* cluster_print_likelihood */

/*============================================================================*/
/* Prevents that empty models are sent out (no associated seqences) by 
   associating a random sequence. Since it's possible to produce an empty model
   in this way, we swap the sequences until a nonempty model is produced. (This 
   could be a never-ending process and therefore it's only done 100 times). */
int cluster_avoid_empty_model(long *seq_label, long seq_number, 
			       int mo_number) {
#define CUR_PROC "cluster_avoid_empty_model"
  int i;
  long j = 0;
  long *counter;
  char error = 1, change = 0;
  int iter = 0;
  /* Initialization */
  if(!m_calloc(counter, mo_number)) {
    mes_proc(); return (-1);
  }
  for (i = 0; i < mo_number; i++)
    counter[i] = 0;
  for (i = 0; i < seq_number; i++) 
    counter[seq_label[i]]++;

  while (error && iter < 100) {
    iter++;
    error = 0;
    /* Do we have empty models ? */
    for (i = 0; i < mo_number; i++) {
      if (counter[i] == 0) {
	change = 1;
	/* Choose a random sequence for an empty model */
	j = m_int(gsl_rng_uniform(RNG) * (seq_number - 1));
	/* The orginal model loses one sequence */
	printf("!!\"avoid_empty_model\":Verschiebe Seq. %ld: %ld --> %d !!\n", 
	       j, seq_label[j], i);
	counter[seq_label[j]] --;
	counter[i] = 1;
	seq_label[j] = i;	
      }  
    }
    /* Is now everything OK ? */
    if (change) {
      for (i = 0; i < mo_number; i++) {
	if (counter[i] <= 0) {
	  error = 1;
	  break;
	}
      }
    }
  } /* while */
  if (error) return (-1);
  else return 0;
#undef CUR_PROC
} /* cluster_avoid_empty_model */

/*============================================================================*/
int cluster_ausgabe(cluster_t *cl, sequence_t *sq, FILE *outfile,
		     char *out_filename) {
#define CUR_PROC "cluster_ausgabe"
  int res = -1;
  int i;
  sequence_d_t *sqd = NULL;

  fprintf(outfile, "\nFinal Models:\n");
  for (i = 0; i < cl->mo_number; i++)
    model_print(outfile, cl->mo[i]);

  res = 0;
STOP:
  sequence_d_free(&sqd);
  return(res);
#undef CUR_PROC
} /* cluster_ausgabe */


/*============================================================================*/
long cluster_update_label(long *oldlabel, long *seq_label, 
			  long seq_number) {
  long i, changes = 0; 
  for (i = 0; i < seq_number; i++) 
    if (oldlabel[i] != seq_label[i]) {
      changes++;
      oldlabel[i] = seq_label[i];
    }
  return changes;
} /* cluster_update_label */
