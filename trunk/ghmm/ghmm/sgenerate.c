/*******************************************************************************
  author       : Bernhard Knab
  filename     : /homes/hmm/wichern/hmm/src/sgenerate.c
  created      : TIME: 09:21:51     DATE: Tue 16. November 1999
  last-modified: TIME: 09:52:15     DATE: Wed 14. March 2001
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "mprintf.h"
#include "mes.h"
#include "const.h"
#include "sequence.h"
#include "smodel.h"
#include "sgenerate.h"
#include "sfoba.h"
#include "matrix.h"
#include "rng.h"


/*============================================================================*/
sequence_d_t *sgenerate_sequences(smodel* smo, int seed, int global_len,
				  long seq_number, long label, int Tmax) {
# define CUR_PROC "sgenerate_sequences"

  /* Endzustand dadurch charakterisiert, dass es keine Ausgangswahrsch. gibt */

  sequence_d_t *sq = NULL;
  int state, n, i, j, m, reject_os, reject_tmax, badseq, class, tilgphase = 0;
  double p, sum, osum = 0.0;
  int len = global_len, up = 0, stillbadseq = 0, reject_os_tmp = 0;

  sq = sequence_d_calloc(seq_number);
  if (!sq) { mes_proc(); goto STOP; }
  if (len <= 0)
    /* Keine spezielle Vorgabe der Sequenzlaenge, Modell sollte 
       Endzustand haben, Verwendung der Konstanten MAX_SEQ_LEN: */
    len = (int)MAX_SEQ_LEN;
  if (Tmax <= 0)
    /* Keine spezielle Vorgabe einer MAX. Sequenzlaenge */
    Tmax = (int)MAX_SEQ_LEN;
  
  /* gsl wird auch von randvar_std_normal "benutzt" 
     seed == -1: Initialisierung fand schon von aussen statt */
  if (seed >= 0) {
    gsl_rng_init();
    if (seed > 0)
      gsl_rng_set(RNG,seed);
    else /* Zufallsinitialisierung! */
      gsl_rng_timeseed(RNG);
  }

  n = 0;
  reject_os = reject_tmax = 0;

  while (n < seq_number) {
    /* Test : fuer jede Seq einen neuen Seed */
    /*   gsl_rng_timeseed(RNG); */
    stillbadseq = badseq = 0;
    if(!m_calloc(sq->seq[n], len)) {mes_proc(); goto STOP;}

    /* Startzustand i wuerfeln */
    p = gsl_rng_uniform(RNG);
    sum = 0.0;
    for (i = 0; i < smo->N; i++) {
      sum += smo->s[i].pi;
      if (sum >= p)
	break;
    }
    if (i == smo->N) { /* kann durch Rundungsfehler in der Eingabe passieren */
      i--;
      while (i > 0 && smo->s[i].pi == 0.0) i--;
    }

    /* Startausgabe erzeugen
       -> m "auswuerfeln", dann nach entspr. pdf omega "auswuerfeln" */
    p = gsl_rng_uniform(RNG);
    sum = 0.0;   
    for (m = 0; m < smo->M; m++) {
      sum += smo->s[i].c[m];
      if (sum >= p)
	break;
    }
    if (m == smo->M) m--;
    /* Zufallszahl nach entsprechender Dichtefunktion */
    sq->seq[n][0] = smodel_get_random_var(smo, i, m);
    state = 1;

    /* Anfangsklasse nach erstem Symbol */
    class = sequence_d_class(sq->seq[n], 0, &osum, &tilgphase);
    while (state < len) {
      /* neuen Zustand i wuerfeln: */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (j = 0; j < smo->s[i].out_states; j++) {
	sum += smo->s[i].out_a[class][j];   
	if (sum >= p)
	  break;
      }
      if (j == smo->s[i].out_states) {/* kann durch Rundungsfehler passieren */
	j--;
	while (j > 0 && smo->s[i].out_a[class][j] == 0.0) j--;
      }
      if (sum == 0.0) {
	if (smo->s[i].out_states > 0) {
	  /* Sequenz verwerfen, falls alle smo->s[i].out_a[class][.] == 0, 
	     d.h. Klasse class ist in Orginaldaten nicht belegt:
	     aus while-Schleife raus, n soll nicht hochgezaehlt werden */
	  /* printf("Zustand %d, class %d, len %d out_states %d \n", i, class,
	     state, smo->s[i].out_states); */
	  badseq = 1;
	  /* break; */
	  
	  /* Versuch: bei "leerer" class die Nachbarklassen probieren; 
	     erst sweep down bis null; falls immer noch ohne Erfolg sweep 
	     up bis COS - 1. Falls immer noch kein Erfolg --> Sequenz
	     verwerfen.
	  */
	  if (class > 0 && up == 0) {
	    class--;
	    continue;
	  }
	  else if (class < COS - 1) {
	    class++;
	    up = 1;
	    continue;
	  }
	  else {
	    stillbadseq = 1;
	    break;
	  }
	}
	else
	  /* Finalzustand erreicht, aus while-Schleife raus: */
	  break;
      }
      i = smo->s[i].out_id[j];

      /* Printing commented out 04.02.01 by Disa */
      /* fprintf(stderr, "%d\n", i); */
      /*      fprintf(stderr, "%d\n", i); */

      /* Ausgabe von Zustand i erzeugen */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (m = 0; m < smo->M; m++) {
	sum += smo->s[i].c[m];
	if (sum >= p)
	  break;
      }
      if (m == smo->M) {
	m--;
	while (m > 0 && smo->s[i].c[m] == 0.0) m--;
      }
      /* Zufallszahl nach entsprechender Dichtefunktion */
      sq->seq[n][state] = smodel_get_random_var(smo, i, m);

      /* class fuer naechsten Schritt bestimmen */
      class = sequence_d_class(sq->seq[n], state, &osum, &tilgphase);
      up = 0;
      state++;

    }  /* while (state < len) */
    if (badseq) {
      reject_os_tmp++;
    }
        
    if (stillbadseq) {
      reject_os++;
      m_free(sq->seq[n]);
      //      printf("cl %d, s %d, %d\n", class, i, n);
    }
    else if (state > Tmax) {
      reject_tmax++;
      m_free(sq->seq[n]);
    }
    else {
      if (state < len)
	if(m_realloc(sq->seq[n], state)) {mes_proc(); goto STOP;}
      sq->seq_len[n] = state;
      sq->seq_label[n] = label;
      //vector_d_print(stdout, sq->seq[n], sq->seq_len[n]," "," ","");
      n++;
    }
    /*    printf("reject_os %d, reject_tmax %d\n", reject_os, reject_tmax); */
    if (reject_os > 10000) {
      mes_prot("Reached max. no. of rejections\n");
      break;
    }
    if (!(n%1000)) printf("%d Seqs. generated\n", n);
  } /* n-loop */

  if (reject_os > 0) printf("%d sequences rejected (os)!\n", reject_os);
  if (reject_os_tmp > 0) printf("%d sequences changed class\n", 
				reject_os_tmp - reject_os);	
  if (reject_tmax > 0) printf("%d sequences rejected (Tmax, %d)!\n", 
			      reject_tmax, Tmax);

  return(sq);
STOP:
  sequence_d_free(&sq);
  return(NULL);
# undef CUR_PROC
} /* sgenerate_sequences */

/*============================================================================*/

/* Extend given sequences on the basis of the model.
   Following modes are possible:
   mode = 0: Initial state: Viterbi, Extension: Viterbi-Path
        = 1: Initial state: Viterbi, Extension: all paths possible
	= 2: Initial state: probability(i) ~= alpha_t(i), 
	     Extension: Viterbi-Path
	= 3: Initial state: probability(i) ~= alpha_t(i),
	     Extension: all paths possible   
	     FRAGE: macht Extension Viterbi ueberhaupt Sinn???
*/



sequence_d_t *sgenerate_extensions(smodel *smo, sequence_d_t *sqd_short, 
				   int seed, int global_len,
				   sgeneration_mode_t mode) {
#define CUR_PROC "sgenerate_extensions"  
  sequence_d_t *sq = NULL;  
  int i, j, t, n, m, len = global_len, short_len, max_short_len = 0, 
    tilgphase = 0, up = 0;
  /* int *v_path = NULL; */
  double log_p, *initial_distribution, **alpha, *scale, p, sum, osum = 0.0;
  /* aicj */
  int class = -1;

  /* TEMP */
  if (mode == all_viterbi || mode == viterbi_viterbi || mode == viterbi_all) {
    mes_prot("Error: mode not implemented yet\n");
    goto STOP;
  }

  if (len <= 0)
    /* no global length; model should have a final state */
    len = (int)MAX_SEQ_LEN;
  max_short_len = sequence_d_max_len(sqd_short);

  /*---------------alloc-------------------------------------------------*/
  sq = sequence_d_calloc(sqd_short->seq_number);
  if (!sq) { mes_proc(); goto STOP; }
  if (!m_calloc(initial_distribution, smo->N)) { mes_proc(); goto STOP; }
  /* is needed in cfoba_forward() */
  alpha = matrix_d_alloc(max_short_len, smo->N);
  if (!alpha) {mes_proc(); goto STOP;}
  if (!m_calloc(scale, max_short_len)) {mes_proc(); goto STOP;}
  gsl_rng_init();
  gsl_rng_set(RNG,seed); 

  /*---------------main loop over all seqs-------------------------------*/
  for (n = 0; n < sqd_short->seq_number; n++) {
    if(!m_calloc(sq->seq[n], len)) {mes_proc(); goto STOP;}
    short_len = sqd_short->seq_len[n];
    if (len < short_len) {
      mes_prot("Error: given sequence is too long\n");
      goto STOP;
    }
    sequence_d_copy(sq->seq[n], sqd_short->seq[n], short_len);
    sq->seq_label[n] = sqd_short->seq_label[n];

    /* Initial distribution */
    /* 1. Viterbi-state */
#if 0
    /* wieder aktivieren, wenn sviterbi realisiert */
    if (mode == viterbi_all || mode == viterbi_viterbi) {
      v_path = cviterbi(smo, sqd_short->seq[n], short_len, &log_p);
      if (v_path[short_len - 1] < 0 || v_path[short_len - 1] >= smo->N) {
	mes_prot("Warning:Error: from viterbi()\n");
	sq->seq_len[n] = short_len; m_realloc(sq->seq[n], short_len);
	continue;
      }
      m_memset(initial_distribution, 0, smo->N);
      initial_distribution[v_path[short_len - 1]] = 1.0; /* all other 0 */
      m_free(v_path);
    }
#endif

    /* 2. Initial Distribution ???
       Pi(i) = alpha_t(i)/P(O|lambda) */
    if (mode == all_all || mode == all_viterbi) {
      if (short_len > 0) {
	if (sfoba_forward(smo, sqd_short->seq[n], short_len, NULL /* ?? */,
			  alpha, scale, &log_p)) { mes_proc(); goto STOP; }	
	sum = 0.0;
	for (i = 0; i < smo->N; i++) {
	  /* alpha ist skaliert! */
	  initial_distribution[i] = alpha[short_len - 1][i]; 
	  sum += initial_distribution[i];
	}
	/* nicht ok.? auf eins skalieren? */
	for (i = 0; i < smo->N; i++)
	  initial_distribution[i] /= sum;      
      }
      else {
	for (i = 0; i < smo->N; i++)
	  initial_distribution[i] = smo->s[i].pi;
      }
    }
    /* if short_len > 0:
       Initial state == final state from sqd_short; no output here
       else
       choose inittial state according to pi and do output
    */
    p = gsl_rng_uniform(RNG);
    sum = 0.0;
    for (i = 0; i < smo->N; i++) {
      sum += initial_distribution[i];
      if (sum >= p)
	break;
    }
    /* error due to incorrect normalization ?? */
    if (i == smo->N) {
      i--;
      while (i > 0 && initial_distribution[i] == 0.0) i--;
    }
    t = 0;
    if (short_len == 0) {
      /* Output in state i */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (m = 0; m < smo->M; m++) {
	sum += smo->s[i].c[m];
	if (sum >= p)
	  break;
      }
      /* error due to incorrect normalization ?? */
      if (m == smo->M) {
	m--;
	while (m > 0 && smo->s[i].c[m] == 0.0) m--;
      }
      sq->seq[n][t] = smodel_get_random_var(smo, i, m);
      class = sequence_d_class(sq->seq[n], t, &osum, &tilgphase); 
      t++;
    }
    /* generate completion for sequence */   
    else {
      for (t = 0; t < short_len; t++)
	class = sequence_d_class(sq->seq[n], t, &osum, &tilgphase); 
      t = short_len;
    }
    while (t < len) {      
      if (smo->s[i].out_states == 0) 
	/* reached final state, exit while loop */
	break;
      sum = 0.0;   
      for (j = 0; j < smo->s[i].out_states; j++) {
	sum += smo->s[i].out_a[class][j];   
	if (sum >= p)
	  break;
      } 
      /* error due to incorrect normalization ?? */
      if (j == smo->s[i].out_states) {
	j--;
	while (j > 0 && smo->s[i].out_a[class][j] == 0.0) j--;
      }
      if (sum == 0.0) {
	/* Versuch: bei "leerer" class die Nachbarklassen probieren; 
	   erst sweep down bis null; falls immer noch ohne Erfolg sweep 
	   up bis COS - 1. Falls immer noch kein Erfolg --> Sequenz
	   verwerfen.
	*/
	if (class > 0 && up == 0) {
	  class--;
	  continue;
	}
	else if (class < COS - 1) {
	  class++;
	  up = 1;
	  continue;
	}
	else {
	  break;
	}
      }
      /* new state */
      i = smo->s[i].out_id[j];

      /* Output in state i */
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (m = 0; m < smo->M; m++) {
	sum += smo->s[i].c[m];
	if (sum >= p)
	  break;
      }
      if (m == smo->M) {
	m--;
	while (m > 0 && smo->s[i].c[m] == 0.0) m--;
      }
      /* random variable from density function */
      sq->seq[n][t] = smodel_get_random_var(smo, i, m);
      class = sequence_d_class(sq->seq[n], t, &osum, &tilgphase); 
      up = 0;
      t++;
    }  /* while (t < len) */    
    if (t < len)
      if(m_realloc(sq->seq[n], t)) {mes_proc(); goto STOP;}
    sq->seq_len[n] = t;    

  } /* for n .. < seq_number */

  matrix_d_free(&alpha, max_short_len);
  m_free(scale);
  return sq;
STOP:
  matrix_d_free(&alpha, max_short_len);
  sequence_d_free(&sq);
  return(NULL);
# undef CUR_PROC
} /* sgenerate_extensions */
 
/*============================================================================*/

double *sgenerate_single_ext(smodel *smo, double *O, const int len, 
			     int *new_len, double **alpha,
			     sgeneration_mode_t mode) {
# define CUR_PROC "sgenerate_single_ext"
  int i, j, m, t, phase, class, up = 0;
  double *new_O = NULL, *scale = NULL, *initial_distribution = NULL;
  double log_p, sum, p, osum;
  int k;
  /* TEMP */
  if (mode == all_viterbi || mode == viterbi_viterbi || mode == viterbi_all) {
    mes_prot("Error: mode not implemented yet\n");
    goto STOP;
  }
  if (len <= 0) {
    mes_prot("Error: sequence with zero or negativ length\n");
    goto STOP;
  }
  if(!m_calloc(new_O, (int)MAX_SEQ_LEN)) {mes_prot("calloc new_O\n"); goto STOP;}
  if(!m_calloc(scale, len)) {mes_prot("calloc scale\n"); goto STOP;}
  if (!m_calloc(initial_distribution, smo->N)) { 
    mes_prot("initial_distribution\n"); goto STOP; } 
  sequence_d_copy(new_O, O, len);
  *new_len = len;
  /* no further extension possible */
  if (O[len - 1] == symbol_kuend || O[len - 1] == symbol_dverz ||
      O[len - 1] == symbol_tilgende) {
    if(m_realloc(new_O, *new_len)) {mes_proc(); goto STOP;}
    return new_O;
  }
    /* Initial Distribution ???
       Pi(i) = alpha_t(i)/P(O|lambda) */
  if (mode == all_all || mode == all_viterbi) {
    if (sfoba_forward(smo, O, len, NULL /* ?? */, alpha, scale, &log_p)) {
      mes_prot("error from sfoba_forward, unable to extend\n"); 
      if(m_realloc(new_O, *new_len)) {mes_proc(); goto STOP;}
      return new_O;
    }
    sum = 0.0;


    for (i = 0; i < smo->N; i++) {
      /* alpha ist skaliert! */
      initial_distribution[i] = alpha[len - 1][i]; 
      sum += initial_distribution[i];
    }
    /* nicht ok.? auf eins skalieren? */
    for (i = 0; i < smo->N; i++) {
      initial_distribution[i] /= sum;                
    }
  }

  p = gsl_rng_uniform(RNG);
  sum = 0.0;
  for (i = 0; i < smo->N; i++) {
    sum += initial_distribution[i];
    if (sum >= p)
      break;
  }
  /* error due to incorrect normalization ?? */
  if (i == smo->N) {
    i--;
    while (i > 0 && initial_distribution[i] == 0.0) i--;
  }  
  /* TEST */
  /* bereits zu Beginn in einem Endzustand? Wie ist das moeglich? */
  if (smo->s[i].out_states == 0) {
    printf("Beginn: Endzustand, State %d\n", i);
    for (k = 0; k < len; k++)
      printf("%.2f ", O[k]);
    printf("\n");
  }
  /* End Test */

  for (t = 0; t < len; t++)
    class = sequence_d_class(O, t, &osum, &phase); 
  t = len;
  while (t < (int)MAX_SEQ_LEN) {  
    if (smo->s[i].out_states == 0) 
	/* reached final state, exit while loop */
	break;
    p = gsl_rng_uniform(RNG);
    sum = 0.0;   
    for (j = 0; j < smo->s[i].out_states; j++) {
      sum += smo->s[i].out_a[class][j];   
      if (sum >= p)
	break;
    } 
    /* error due to incorrect normalization ?? */
    if (j == smo->s[i].out_states) {
      j--;
      while (j > 0 && smo->s[i].out_a[class][j] == 0.0) j--;
    }
    if (sum == 0.0) {
      /* Versuch: bei "leerer" class die Nachbarklassen probieren; 
	 erst sweep down bis null; falls immer noch ohne Erfolg sweep 
	 up bis COS - 1. Falls immer noch kein Erfolg --> Sequenz
	 verwerfen.
      */
      if (class > 0 && up == 0) {
	class--;
	continue;
      }
      else if (class < COS - 1) {
	class++;
	up = 1;
	continue;
      }
      else {
	char *str = mprintf(NULL,0,"unable to extend seq (all osc empty)\n"); 
	mes_prot(str); m_free(str);
	goto STOP;
      }      
    } /* sum == 0 */

    /* new state */
    i = smo->s[i].out_id[j];

    if (smo->M == 1)
      m = 0;
    else {            
      p = gsl_rng_uniform(RNG);
      sum = 0.0;   
      for (m = 0; m < smo->M; m++) {
	sum += smo->s[i].c[m];
	if (sum >= p)
	  break;
      }
      if (m == smo->M) {
	m--;
	while (m > 0 && smo->s[i].c[m] == 0.0) m--;
      }
    }
    /* Output in state i, komp. m */
    /* random variable from density function */
    new_O[t] = smodel_get_random_var(smo, i, m);
    class = sequence_d_class(new_O, t, &osum, &phase); 
    t++;
    up = 0;
  }  /* while (t < MAX_SEQ_LEN) */ 
  if (t < (int) MAX_SEQ_LEN) 
    if(m_realloc(new_O, t)) {mes_proc(); goto STOP;}
  
  *new_len = t;    
  

  m_free(scale);
  m_free(initial_distribution);

  return new_O;
STOP:
  m_free(new_O);
  m_free(scale);
  m_free(initial_distribution); 
  return NULL;
# undef CUR_PROC
} /* sgenerate_single_ext */


/* generate a single next value bases on a trained model and on a seq und
   to length "len". Use the most prob. state given the seq as an initial state
   and determin the next state und the symbol with the RNG.

   wird diese Fkt. fuer EIN O und aufeinanderfolg. "len" aufgerufen, so kann
   noch deutlich optimiert werden
*/

double sgenerate_next_value(smodel *smo, double *O, const int len) {
# define CUR_PROC "sgenerate_next_value"
  double **alpha = NULL;
  double res = -1.0, sum, p;
  double *scale = NULL, log_p, max_val = -1000000;
  int i, j, m, init_state = -1;

  if (COS > 1) {
    mes_prot("sgenerate_next_value only for COS == 1\n");
    goto STOP;
  }

  alpha = matrix_d_alloc(len, smo->N);
  if (!alpha) {mes_proc(); goto STOP;}
  if(!m_calloc(scale, len)) {mes_prot("calloc scale\n"); goto STOP;}
  if (sfoba_forward(smo, O, len, NULL /* ?? */, alpha, scale, &log_p)) {
      mes_prot("error from sfoba_forward\n"); 
      goto STOP;
  }

  /* find inititial state */
  sum = 0.0;
  for (i = 0; i < smo->N; i++) 
    sum += alpha[len - 1][i];
  if ( sum < 0.9 || sum > 1.1) {
    printf("Error sum = %.4f (!= 1)\n", sum);
    goto STOP;
  }
  /* max state */  
  for (i = 0; i < smo->N; i++) {
    if (alpha[len - 1][i] > max_val) {
      init_state = i;
      max_val = alpha[len - 1][i];     
    }
  }

  /* random state */
  /*
    p = gsl_rng_uniform(RNG);
    sum = 0.0;
    for (i = 0; i < smo->N; i++) {
    sum += alpha[len - 1][i];
    if (sum >= p)
    break;    
    }   
    if (i == smo->N) {
    i--;
    while (i > 0 && alpha[len - 1][i] == 0.0) i--;
    }  
    init_state = i;
  */


  if (init_state == -1 || smo->s[init_state].out_states == 0) goto STOP;

  p = gsl_rng_uniform(RNG);
  sum = 0.0;   
  for (j = 0; j < smo->s[init_state].out_states; j++) {
    sum += smo->s[init_state].out_a[0][j];   
    if (sum >= p)
      break;
  } 
  /* error due to incorrect normalization ?? */
  if (j == smo->s[init_state].out_states) {
    j--;
    while (j > 0 && smo->s[init_state].out_a[0][j] == 0.0) j--;
  }

  /* new state */
  i = smo->s[init_state].out_id[j];
  
  if (smo->M == 1)
    m = 0;
  else {            
    p = gsl_rng_uniform(RNG);
    sum = 0.0;   
    for (m = 0; m < smo->M; m++) {
      sum += smo->s[i].c[m];
      if (sum >= p)
	break;
    }
    if (m == smo->M) {
      m--;
      while (m > 0 && smo->s[i].c[m] == 0.0) m--;
    }
  }
  /* Output in state i, komp. m */
  /* random variable from density function */
  res = smodel_get_random_var(smo, i, m);

 STOP:
 matrix_d_free(&alpha, len);
 m_free(scale);
  return res;
# undef CUR_PROC
}

