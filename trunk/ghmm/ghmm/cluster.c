/*******************************************************************************
  author       : Bernd Wichern
  filename     : /homes/hmm/bknab/src/cluster.c
  created      : TIME: 10:55:33     DATE: Tue 02. June 1998
  last-modified: TIME: 10:53:16     DATE: Wed 10. March 1999
*******************************************************************************/
/* $Id: */ 
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


#if 0
/*============================================================================*/
int main(int argc, char* argv[]) {
#define CUR_PROC "main"
  int exitcode = -1;
  if (argc == 4 || argc == 5) {
    if (argc == 5) {
      int j;
      gsl_rng_set(RNG,atoi(argv[4]));
      for (j = 0; j < 100; j++)
	printf("%d \n", m_int(gsl_rng_uniform(RNG) * 10));
      exit(1);
    }
    else
      gsl_rng_set(RNG,0);
    exitcode = cluster_hmm(argv[1], argv[2], argv[3]);
  }
  else {
    mes_prot
      ("Insufficient arguments. \
        Usage: cluster [sequence file][model file][outfile] <seed>\n"); 
  }
  return exitcode;
# undef CUR_PROC
} /* main */

#endif /* 0 */

/*============================================================================*/
int cluster_hmm(char *seq_file, char *mo_file, char *out_filename)  {
# define CUR_PROC "cluster_hmm"
  int res = -1, i, iter = 0;
  sequence_t *sq = NULL;
  long j, changes = 1; /* muss > 0 ! */
  long *oldlabel;
  double log_p;
  FILE *outfile = NULL;
  cluster_t cl;
  cl.mo = NULL;
  cl.mo_seq = NULL;  
  if(!(outfile = mes_fopen(out_filename, "wt"))) {mes_proc(); goto STOP;}

/*----------------------------------------------------------------------------*/

  /* Speicher alloc und Daten einlesen: Sequenzen und Initialmodelle */
  sq = sequence_read(seq_file);
  if (!sq) {mes_proc(); goto STOP;}
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

    /* Sequenzen zuordnen */
    fprintf(outfile, "\nSequence, Best Model, logP of generating Seq.:\n");
    for (j = 0; j < sq->seq_number; j++) {
      sq->seq_label[j] = 
	sequence_best_model(cl.mo, cl.mo_number, sq->seq[j], sq->seq_len[j],
			    &log_p);    
      fprintf(outfile, "seq %ld, mo %ld, log p %.4f\n", j, 
	      sq->seq_label[j], log_p);
      if (sq->seq_label[j] == -1 || sq->seq_label[j] >= cl.mo_number) { 
	/* kein Modell passt, was tun ? */
	char *str =
	  mprintf(NULL, 0, "Seq. %ld: sequence_best_model liefert %d\n",
		  j, sq->seq_label[j]); 
	mes_prot(str); m_free(str); goto STOP;
      }
    }
    if (cluster_avoid_empty_model(sq->seq_label, sq->seq_number, cl.mo_number))
      { mes_proc(); goto STOP; }
    changes = cluster_update_label(oldlabel, sq->seq_label, sq->seq_number);
    fprintf(outfile, "%ld changes\n", changes);
    fprintf(stdout, "\n*** %ld changes in iteration %d ***\n\n", changes, iter);

    /* Modelle mit den zugeordneten Sequenzen Reestimieren */
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
/* NEU: Speicher fuer Sequenzen fuer jedes Modell nur einmal allocieren und
   nicht wie vorher fuer jede Sequenz mit realloc arbeiten.*/
int cluster_update(cluster_t *cl, sequence_t *sq) {
#define CUR_PROC "cluster_update"
  int i, res = -1;
  long *seq_counter;
  sequence_t *seq_t;
  if(!m_calloc(seq_counter, cl->mo_number)) {mes_proc(); goto STOP;}
  /* Anzahl zugeordneter Sequenzen feststellen */
  for (i = 0; i < sq->seq_number; i++) 
    seq_counter[sq->seq_label[i]]++;
  /* Speicher blockweise allocieren */
  for (i = 0; i < cl->mo_number; i++) {
    if (cl->mo_seq[i]) {
      /* wichtig: hier kein sequence_free, sonst sind auch die Original
	 Sequenzen futsch */
      sequence_clean(cl->mo_seq[i]);
      m_free(cl->mo_seq[i]);
    }
    cl->mo_seq[i] = sequence_calloc(seq_counter[i]);
    cl->mo_seq[i]->seq_number = 0; /* wird unten hochgezaehlt */
  }
  /* Eintraege setzen */
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
/* Verhindert, dass Modelle leer ausgehen (keine Sequenz zugeordnet), 
   indem ihnen eine zufaellige Sequenz zugeordnet wird. Da hierdurch 
   evt. erneut leere Modelle erzeugt werden, werden Sequenzen getauscht,
   bis keine leeren Modelle vorliegen. (Gefahr einer Endlosschleife, daher
   Abbruch nach 100 Interationen) */
int cluster_avoid_empty_model(long *seq_label, long seq_number, 
			       int mo_number) {
#define CUR_PROC "cluster_avoid_empty_model"
  int i;
  long j = 0;
  long *counter;
  char error = 1, change = 0;
  int iter = 0;
  /* Initialisierungen */
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
    /* gibt es leere Modelle ? */
    for (i = 0; i < mo_number; i++) {
      if (counter[i] == 0) {
	change = 1;
	/* zuf. Sequenz fuer leeres Modell auswaehlen */
	j = m_int(gsl_rng_uniform(RNG) * (seq_number - 1));
	/* urspruengliches Modell verliert eine Sequenz */
	printf("!!\"avoid_empty_model\":Verschiebe Seq. %ld: %ld --> %d !!\n", 
	       j, seq_label[j], i);
	counter[seq_label[j]] --;
	counter[i] = 1;
	seq_label[j] = i;	
      }  
    }
    /* jetzt alles ok ? */
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
