/*******************************************************************************
  author       : Bernhard Knab
  filename     : ghmm/tools/scluster.c
  created      : 2001-04-20 by Achim Gaedke from hmm/src/scluster.c
  $Id$

__copyright__

*******************************************************************************/

#ifdef WIN32
#  include "win_config.h"
#endif

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdio.h>
#include <ghmm/mes.h>
#include <ghmm/rng.h>
#include <ghmm/scluster.h>

int main(int argc, char* argv[]) {
#define CUR_PROC "scluster_main"
  int exitcode = -1;

  ghmm_rng_init();

  if (argc == 5 || argc == 6) {
    if (argc == 6)
      GHMM_RNG_SET(RNG,atoi(argv[5]));
    else {
      /* random init */
      ghmm_rng_timeseed(RNG);
    }

    printf("Clustering Sequences with start partition\n");
    switch(atoi(argv[4])) {
    case 0: printf("SP_BEST (best model)\n"); break;
    case 1: printf("NO_SP (no start partition)\n"); break;
    case 2: printf("SP_KM (partition from k-means)\n"); break;
    case 3: printf("SP_ZUF (random start partition)\n"); break;
    default: printf("argv[4] %d not valid. must be in [0, 3]\n", atoi(argv[4]));
      return exitcode;
    }
    exitcode = scluster_hmm(argv);
  }
  else {
    mes_prot
      ("Insufficient arguments. Usage: scluster [sequence file][model file][outfile][labels]<seed>\n"); 
  }
  /*------------------------------------------------------------------------*/
  mes(MES_WIN, "\n(%2.2T): Program finished with exitcode %d.\n", exitcode );
  mes_exit();
  return(exitcode);
# undef CUR_PROC
} /* main */
