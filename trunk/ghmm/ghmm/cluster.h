/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/cluster.h
  created      : TIME: 11:17:55     DATE: Tue 02. June 1998
  $Id$ 

__copyright__

*********************************************************************************/




#ifndef CLUSTER_H
#define CLUSTER_H
#include "sequence.h"
#include "model.h"

/**
 @name cluster
 */

//@{

/**
 */
struct cluster_t{
  /** 
    Vektor von HMM-Modell Pointern */
  model **mo;
  /** 
    Vektor von sequence_t Pointern: zur Speicherung von zu Modellen
    gehoerenden Sequenz Daten */
  sequence_t **mo_seq;
  /** 
    Anzahl der eingelesenen Modelle */
  int mo_number;
};
///
typedef struct cluster_t cluster_t;

///
int cluster_ausgabe(cluster_t *cl, sequence_t *sq, FILE *outfile,
		    char *out_filename);

///
int cluster_avoid_empty_model(long *best_model, long seq_number, int mo_number);

///
int cluster_hmm(char *seq_file, char *mo_file, char *out_file);

///
int cluster_update(cluster_t *cl, sequence_t *sq);

///
long cluster_update_label(long *oldlabel, long *seq_label, long seq_number);

///
void cluster_print_likelihood(FILE *outfile, cluster_t *cl);

//@} cluster documentation

#endif /* CLUSTER_H */
