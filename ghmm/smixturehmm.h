/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/smixturehmm.h
  created      : TIME: 18:07:43     DATE: Thu 20. July 2000
  $Id$

__copyright__

*******************************************************************************/

#ifndef SMIXTUREHMM_H
#define SMIXTUREHMM_H

#ifdef __cplusplus
extern "C" {
#endif

#include <ghmm/sequence.h>
#include <ghmm/smodel.h>

/**
   @name mixture-shmm methods
 */

/*@{ smixturehmm section
 */

/**
 */
int smixturehmm_cluster(FILE *outfile, double **cp, sequence_d_t *sqd, 
			smodel **smo, int smo_number);

/**
 */
double smixturehmm_like(smodel **smo, int  smo_number, sequence_d_t *sqd_test,
			long *errors);

/**
 */
int smixturehmm_init(double **cp, sequence_d_t *sqd, smodel **smo,
			      int smo_number, int mode);

/**
 */
int smixturehmm_calc_priors(double **cp, sequence_d_t *sqd, smodel **smo,
			    int smo_number);

/**
 */
int smixturehmm_calc_cp(double **cp, sequence_d_t *sqd, smodel **smo, 
			int smo_number, double *total_train_w);

/**
*/
void smixture_calc_logp(double **logp, int **error, sequence_d_t *sqd, 
			smodel **smo,  int smo_number);

/**
 */
void smixturehmm_print_header(FILE *file, char *argv[], int flag);


/**
 */
double *smixturehmm_avg_like(double **cp, sequence_d_t *sqd, 
			     smodel **smo, int smo_number);
#ifdef __cplusplus
}
#endif

/*@} smixturehmm section */

#endif /* SMIXTUREHMM_H */
