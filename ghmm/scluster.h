/*-----------------------------------------------------------------------------
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/scluster.h
  created      : TIME: 15:53:53     DATE: Tue 16. November 1999
  last-modified: TIME: 11:25:33     DATE: Tue 09. January 2001
------------------------------------------------------------------------------*/


#ifndef SCLUSTER_H
#define SCLUSTER_H
#include "sequence.h"
#include "sreestimate.h"

/**
   @name scluster
 */

//@{

/* switsch for parallel mode: (1) = sequential, (0) = parallel */ 
#define POUT 1
/* number of parallel threads */
#define THREADS 4

#define CLASSIFY 0 /* Switch for Classificator: 0 == MD, 1 == MAW */

/**
   Clusterstruktur: Alle Modelle und Sequencen.. */
struct scluster_t{
  /** 
  Vektor von SHMM-Modell Pointern */
  smodel **smo;
  /** 
    Vektor von sequence_t Pointern: zur Speicherung von zu Modellen
    gehoerenden Sequenz Daten */
  sequence_d_t **smo_seq;
  /** 
    Anzahl der eingelesenen Modelle */
  int smo_number;
  /** 
    Anzahl der Sequenzen pro Modell */
  long *seq_counter;
  /** 
    log(P) der Modelle */
  double *smo_Z_MD;
  /** a posteriore WS der Modelle zur Berechnung der Zielfunktion
      im Fall des MAW-Klassifikators. Wird berechnet mit smap_bayes */
  double *smo_Z_MAW;
};
///
typedef struct scluster_t scluster_t;


///
int scluster_out(scluster_t *cl, sequence_d_t *sqd, FILE *outfile,
		 char *argv[]);

///
int scluster_avoid_empty_smodel(sequence_d_t *sqd, scluster_t *cl);

///
int scluster_hmm(char *argv[]);

///
int scluster_update(scluster_t *cl, sequence_d_t *sqd);

///
long scluster_update_label(long *oldlabel, long *seq_label, long seq_number, 
			   long *smo_changed);

///
void scluster_print_likelihood(FILE *outfile, scluster_t *cl);

///
int scluster_calc_print_indices(FILE *outfile, scluster_t *cl);

///
int scluster_best_model(scluster_t *cl, long seq_id, double **all_log_p,
			double *log_p);

///
void scluster_prob(smosqd_t *cs);

/** int scluster_labels_from_kmeans(sequence_d_t *sqd, int smo_number); */

int scluster_random_labels(sequence_d_t *sqd, int smo_number);

/** calculates the aposteriori prob. $\log(p(\lambda_best | O[seq\_id]))$, 
    where $\lambda_best$ is the model with highest apost. prob. 
*/
int  scluster_log_aposteriori(scluster_t *cl, sequence_d_t *sqd, int seq_id, 
			      double *log_apo);

///
void scluster_print_header(FILE *file, char* argv[]);

//@}

#endif /* SCLUSTER_H */
