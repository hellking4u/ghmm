/*******************************************************************************
  author       : Bernd Wichern
  filename     : /homes/hmm/wichern/hmm/src/smixturehmm.h
  created      : TIME: 18:07:43     DATE: Thu 20. July 2000
  last-modified: TIME: 11:17:22     DATE: Fri 22. December 2000
*******************************************************************************/

#ifndef SMIXTUREHMM_H
#define SMIXTUREHMM_H
#include "sequence.h"
#include "smodel.h"

/**
   @name mixture-shmm methods
 */

//@{

///
int smixturehmm_cluster(FILE *outfile, double **cp, sequence_d_t *sqd, 
			smodel **smo, int smo_number);

///
double smixturehmm_like(smodel **smo, int  smo_number, sequence_d_t *sqd_test,
			long *errors);

///
int smixturehmm_init(double **cp, sequence_d_t *sqd, smodel **smo,
			      int smo_number, int mode);

///
int smixturehmm_calc_priors(double **cp, sequence_d_t *sqd, smodel **smo,
			    int smo_number);

///
int smixturehmm_calc_cp(double **cp, sequence_d_t *sqd, smodel **smo, 
			int smo_number);

///
void smixture_calc_logp(double **logp, int **error, sequence_d_t *sqd, 
			smodel **smo,  int smo_number);

///
void smixturehmm_print_header(FILE *file, char *argv[], int flag);

///
double *smixturehmm_tilg_w(double **cp, sequence_d_t *sqd, 
			   smodel **smo, int smo_number);

///
double *smixturehmm_avg_like(double **cp, sequence_d_t *sqd, 
			     smodel **smo, int smo_number);

//@} smixturehmm section

#endif /* SMIXTUREHMM_H */
