/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/tools/cluster.c
  created      : DATE: March 2001 by Achim Gaedke from hmm/src/cluster.c
  $Id$
*****************************************************************************/

#include <stdio.h>
#include <ghmm/mes.h>
#include <ghmm/cluster.h>

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


