/* 
   File:        probdist.c

   synopsis:    probdist model1.hmm model2.hmm

   options:     -t <int>   Sequence length 
                -s         symmetric 

   description: computes a probabilistic distance between the two 
                models

                cf. Juang & Rabiner "A Probabilistic Distance Measure
                for Hidden Markov Models"

   History:

   3/18/99 Initial version AS
*/

#include <math.h>

#include "stdmacro.h"
#include "const.h"
#include "foba.h"
#include "matrix.h"
#include "model.h"
#include "sequence.h"
#include "viterbi.h"
#include "rng.h"
#include "cfoba.h"
#include "cmodel.h"
#include "smodel.h"


int main(int argc, char* argv[]) {
# define CUR_PROC "main"

  int mo_number, cmo_number, smo_number;
  int discrete = 0, smodelflag = 0;
  int T, i, j;
  model **mo = NULL;
  cmodel **cmo = NULL;
  smodel **smo = NULL;
  double d;
  /*  double dangle, ddiff, s1, s2; */

  if (argc != 5) {
    printf("Usage: main <Model File> <discrete-flag> <T> <smodel flag>\n");
    exit(1);
  }

  discrete = atoi(argv[2]);
  T = atoi(argv[3]);
  smodelflag = atoi(argv[4]);

  gsl_rng_init();
  if (smodelflag) {
    smo = smodel_read(argv[1], &smo_number);
    if (!smo)  {mes_proc(); return -1;}  
    if (smo_number < 2) {
      printf("Need at least two HMMs to compare (read %d)\n", smo_number);
      return -1;
    }
    for (i = 0; i < smo_number - 1; i++)
      for (j = i + 1; j < smo_number; j++) {
	printf("#----- mo[%d], mo[%d] \n", i , j);
	/* syntax prob_dist: (smo1, smo2, total seqlen., symmetric, verbose) */
	d = smodel_prob_distance(smo[i], smo[j], T, 1, 0);
	printf("probdist = %f\n",d);

      }

  }  
  else if (discrete) {

    mo = model_read(argv[1], &mo_number);
    if (!mo) {mes_proc(); return -1;}        

    if (mo_number < 2) {
      printf("Need at least two HMMs to compare\n");
      return -1;
    }

    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[1], mo[0] \n");
    d = model_prob_distance(mo[1],mo[0], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- mo[1], mo[0] \n");
    d = model_prob_distance(mo[1],mo[0], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 1, 1);
    printf("d=%f\n",d);
    
    printf("#----- mo[0], mo[1] \n");
    d = model_prob_distance(mo[0],mo[1], T, 1, 0);
    printf("d=%f\n",d);
    
  }    
  else {
    
    cmo = cmodel_read(argv[1], &cmo_number);
    if (!cmo) {mes_proc(); return -1;}        
    
    if (cmo_number < 2) {
      printf("Need at least two CHMMs to compare\n");
      return -1;
    }
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[1], cmo[0] \n");
    d = cmodel_prob_distance(cmo[1], cmo[0], T, 0, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- cmo[1], cmo[0] \n");
    d = cmodel_prob_distance(cmo[1], cmo[0], T, 0, 0);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 1, 1);
    printf("d=%f\n",d);
    
    printf("#----- cmo[0], cmo[1] \n");
    d = cmodel_prob_distance(cmo[0], cmo[1], T, 1, 0);
    printf("d=%f\n",d);
    
#if 0
    /* coemission likelihood */
    printf("#----- cmo[0]/cmo[1] \n");
    if (cmodel_coemission_likelihood(cmo[0], cmo[1], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    printf("#----- cmo[0]/cmo[0] \n");
    if (cmodel_coemission_likelihood(cmo[0], cmo[0], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    printf("#----- cmo[1]/cmo[1] \n");
    if (cmodel_coemission_likelihood(cmo[1], cmo[1], &d) == -1) d = -1;
    printf("Coemission Likelihood = %e\n",d); 
    
    printf("#----- D_angle, D_diff, S1, S2\n");
    cmodel_measures(cmo[0], cmo[1], &dangle, &ddiff, &s1, &s2);
    printf("D_angle = %e\n", dangle);
    printf("D_diff  = %e\n", ddiff); 
    printf("S1      = %e\n", s1);
    printf("S2      = %e\n", s2); 
#endif

  }
  
  return 0;
}
