/*******************************************************************************
  author       : Bernhard Knab
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/scluster.c
  created      : TIME: 15:47:54     DATE: Tue 16. November 1999
*******************************************************************************/

#include <stdio.h>
#include <pthread.h>
#include <float.h>
#include <math.h>
#include "mprintf.h"
#include "mes.h"
#include "scluster.h"
#include "smodel.h"
#include "sreestimate.h"
#include "sfoba.h"
#include "rng.h"
#include "sequence.h"
#include "const.h"
#include "matrix.h"
#include "vector.h"


#if 0
/*============================================================================*/
int main(int argc, char* argv[]) {
#define CUR_PROC "scluster_main"
  int exitcode = -1;

  gsl_rng_init();

  if (argc == 5 || argc == 6) {
    if (argc == 6)
      gsl_rng_set(RNG,atoi(argv[5]));
    else {
      /* random init */
      gsl_rng_timeseed(RNG); //previously: gsl_rng_set(RNG,0);
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

#endif /* 0 */

/*============================================================================*/

int scluster_hmm(char* argv[]) {
# define CUR_PROC "scluster_hmm"
  //  int res = -1, i, k, iter = 0, sqd_number;
  char *seq_file = argv[1],  *smo_file = argv[2], *out_filename = argv[3];
  int labels = atoi(argv[4]);
  int res = -1, i, k, iter = 0, sqd_number, idummy;
  sequence_d_t *sqd = NULL;
  sequence_d_t **sqd_vec = NULL; /* only temp. pointer */
  long j, changes = 1; 
  long *oldlabel, *smo_changed;
  double log_p, log_apo;
  double **all_log_p = NULL; /* for storing all log_p */
  FILE *outfile = NULL;
  char *str;
  char filename[1024];
  scluster_t cl;
  double eps_bw, q;
  int max_iter_bw;
  /* sreestimate_baum_welch needs this structure (introduced for parallel mode) */
  smosqd_t *cs; 
  double *tilgw, *model_weight;
  int *tid, *return_value;
  pthread_attr_t Attribute;
  cl.smo = NULL;
  cl.smo_seq = NULL;
  cl.seq_counter = NULL;
  cl.smo_Z_MD = NULL;
  cl.smo_Z_MAW = NULL;

  /* -----------------initialization -------------------------*/
  fprintf(stdout, "Clustering seqs. \"%s\"\nwith ", seq_file);
  if (CLASSIFY == 0)    fprintf(stdout, "MD-Classification\n");
  else                              fprintf(stdout, "MAW-Classification\n");
  fflush(stdout);
  /*-----------open output file----------------------------------*/
  sprintf(filename, "%s%s", out_filename, ".cl");
  if(!(outfile = mes_fopen(filename, "wt"))) {mes_proc(); goto STOP;}

/*--------------Speicher alloc und Daten einlesen----------------------------*/

  /* 1 Sequenzfeld und Initialmodelle */
  scluster_print_header(outfile, argv);
  /*--- memory alloc and read data ----------------------------*/
  sqd_vec = sequence_d_read(seq_file, &sqd_number);
  if (!sqd_vec) {mes_proc(); goto STOP;}
  sqd = sqd_vec[0];
  cl.smo = smodel_read(smo_file, &cl.smo_number);
  if (!cl.smo) {mes_proc(); goto STOP;}
  /* if labels == 0 || labels == 1: keine Startlabels benoetigt. */
  /* if `labels` =  3: random start labels */
  if (labels == 3) {
    fprintf(outfile, "(random start partition)n");
    fprintf(stdout, "(Random start partition...)\n");
    scluster_random_labels(sqd, cl.smo_number);
  }  

  if(!m_calloc(oldlabel, sqd->seq_number)) {mes_proc(); goto STOP;}
  for (i = 0; i < sqd->seq_number; i++)
    oldlabel[i] = (-1);
  if(!m_calloc(cl.smo_seq, cl.smo_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(cl.seq_counter, cl.smo_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(cl.smo_Z_MD, cl.smo_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(cl.smo_Z_MAW, cl.smo_number)) {mes_proc(); goto STOP;}
  all_log_p = matrix_d_alloc(cl.smo_number, (int) sqd->seq_number);
  if (!all_log_p) {mes_proc(); goto STOP;}
  if (smodel_check_compatibility(cl.smo, cl.smo_number)) { 
    mes_proc(); goto STOP; 
  } 
  if(!m_calloc(smo_changed, cl.smo_number)) {mes_proc(); goto STOP;}
  for (i = 0; i < cl.smo_number; i++) {
    cl.smo_seq[i] = NULL;
    smo_changed[i] = 1;
  }
  /* uninformative model prior */
  if (CLASSIFY == 1) 
    for (i = 0; i < cl.smo_number; i++) 
      cl.smo[i]->prior = 1/(double) cl.smo_number;
  /*--------for parallel mode  --------------------*/
  /* id for threads */
  if(!m_calloc(tid, cl.smo_number)) {mes_proc(); goto STOP;} 
  /* data structure for  threads */
  if(!m_calloc(cs, cl.smo_number)) {mes_proc(); goto STOP;} 
  for (i = 0; i < cl.smo_number; i++)
    cs[i].smo = cl.smo[i]; 
  /* returnvalues for each thread */
  if(!m_calloc(return_value, cl.smo_number)) {mes_proc(); goto STOP;} 
  pthread_attr_init(&Attribute);
  pthread_attr_setscope(&Attribute, PTHREAD_SCOPE_SYSTEM);	

  /* eps_bw: stopping criterion in baum-welch
     max_iter_bw: max. number of baum-welch iterations */
  eps_bw = 0.0;    // not used
  q = 0.1;         // ???
  max_iter_bw = 20; //MAX_ITER_BW;
  
/*------------------------main loop-------------------------------------------*/
  /* do it until no sequence changes model; 
     attention: error function values are not used directly as a stopping criterion */
  while(changes > 0 || eps_bw > EPS_ITER_BW || max_iter_bw < MAX_ITER_BW) { 
    iter++;
    /* reset error functions and counters */
    for (i = 0; i < cl.smo_number; i++) {
      cl.seq_counter[i] = 0;
      cl.smo_Z_MD[i] = 0.0;
      cl.smo_Z_MAW[i] = 0.0;
    }

    /* set pointer of cs to sqd-field and matrix of all_logp values */
    for (i = 0; i < cl.smo_number; i++) {      
      cs[i].sqd = sqd; /* all seqs. ! */
      cs[i].logp = all_log_p[i];
    }
    /* ------------calculate logp for all seqs. and all models --------------*/
#if POUT
    /* sequential version */
    for (i = 0; i < cl.smo_number; i++) {
      if (!smo_changed[i]) continue;
      scluster_prob(&cs[i]);
    }
#else 
    /* parallel version */
    i = 0;
    while (i < cl.smo_number) {
      for (j = i; j < m_min(cl.smo_number,i+THREADS); j++) {
	if (!smo_changed[j]) continue;
	if (perror = pthread_create(&tid[j],&Attribute,scluster_prob,&cs[j])){
	  str = mprintf(NULL,0,
		"pthread_create returns false (step %d, smo[%d])\n",iter,j); 

	  mes_prot(str); m_free(str); goto STOP; 
	}
      }
      for (j = i; j < m_min(cl.smo_number,i+THREADS); j++) {
	if (smo_changed[j]) pthread_join(tid[j], NULL);
      }
      i = j;
    }
#endif


    if (iter == 1 && labels > 1) {
      /* use sequence labels as start labels */
      fprintf(outfile, "(sequences with start labels!)\n");
      fprintf(stdout, "(sequences with start labels!)\n");
    } 
    
    /*----------- classification,  sequence labeling and error function -----------*/
    if (iter == 1 && labels ==1) {
      for (i = 0; i < cl.smo_number; i++) {
	cl.seq_counter[i] = sqd->seq_number;
      }
    }
    else {
    for (j = 0; j < sqd->seq_number; j++) {
      if (iter > 1 || labels == 0) 
	/* classification: set seq_label to ID of best_model  */
	sqd->seq_label[j] = scluster_best_model(&cl, j, all_log_p, &log_p);
      if (sqd->seq_label[j] == -1 || sqd->seq_label[j] >= cl.smo_number) { 
	/* no model fits! What to do?  hack: use arbitrary model ! */
	str =  mprintf(NULL, 0, "Warning: seq. %ld, ID %.0f: scluster_best_model returns %d\n",
		       j, sqd->seq_id[j], sqd->seq_label[j]); 
	mes_prot(str); m_free(str);
	sqd->seq_label[j] = j % cl.smo_number;
	/* goto STOP; */
      }
      cl.seq_counter[sqd->seq_label[j]]++; 
      /* add to error function value with seq. weights  */
      /* 1. Z_MD */
      cl.smo_Z_MD[sqd->seq_label[j]] += sqd->seq_w[j] * 
	all_log_p[sqd->seq_label[j]][j];  
      /* 2. Z_MAW */
      if (CLASSIFY == 1) { 
	idummy = scluster_log_aposteriori(&cl, sqd, j, &log_apo);
	if (idummy == -1) {
	  str = mprintf(NULL, 0 , "Warn: no model fits to Seq %10.0f, use  PENALTY_LOGP\n",
			sqd->seq_id[j]);	     
	  mes_prot(str); m_free(str);
	  cl.smo_Z_MAW[sqd->seq_label[j]] += sqd->seq_w[j] * PENALTY_LOGP;
	  continue;
	}
	if (idummy != sqd->seq_label[j]) {
	  printf("Warn: best_model = %ld; best_bayes(logapo) = %d\n", 
		  sqd->seq_label[j], idummy);	
	}
	cl.smo_Z_MAW[sqd->seq_label[j]] += sqd->seq_w[j] * log_apo;
      }
    }   /* for (j = 0; j < sqd->seq_number; j++) */
    }  /* else */

    if (!(iter == 1 && labels == 1)) {
      if (scluster_avoid_empty_smodel(sqd, &cl)== -1) {mes_proc(); goto STOP;}
      for (i = 0; i < cl.smo_number; i++) smo_changed[i] = 0;
      changes = scluster_update_label(oldlabel, sqd->seq_label, sqd->seq_number,
				      smo_changed);
    }
    fprintf(outfile, "%ld changes in iteration %d \n", changes, iter);
    fprintf(stdout, "\n*** %ld changes in iteration %d ***\n", changes, iter);

    /* NEU: Falls keine Wechsel mehr: max_iter erhoehen 
       (Alternative: eps veraendern */
    if (changes == 0 && max_iter_bw < MAX_ITER_BW) {
      max_iter_bw = MAX_ITER_BW;
      changes = 1; // damit wieder reestimiert und zugeordnet wird 
      for (i = 0; i < cl.smo_number; i++) smo_changed[i] = 1; // " 
      fprintf(outfile, "max_iter_bw := %d\n", max_iter_bw);
      fprintf(stdout, "*** max_iter_bw := %d ***\n", max_iter_bw);
    }
    
    /* -------Reestimate all models with assigned sequences---------------- */
    if (changes > 0) {
      if (!(iter == 1 && labels == 1))
	if (scluster_update(&cl, sqd)) { mes_proc(); goto STOP; }
      
      /* TEST: k-Means-Start-Clusterung rausschreiben */
      if (iter == 1 && labels == 2) {
	for (i = 0; i < cl.smo_number; i++) {
	  if (cl.smo_seq[i] != NULL)
	    sequence_d_print(stdout, cl.smo_seq[i], 0);  
	}
      }

      /*if (labels == -1) {
	fprintf(outfile, "(only determining the best model for each seq.!)\n");
	fprintf(stdout, "(only determining the best model for each seq.!)\n");
	break;
      }
	*/      
      
      fprintf(outfile, "\ntotal Prob. before %d.reestimate:\n", iter);
      scluster_print_likelihood(outfile, &cl);
      
      for (i = 0; i < cl.smo_number; i++) {
	cs[i].smo = cl.smo[i]; 
	if (!(iter == 1 && labels == 1))
	  cs[i].sqd = cl.smo_seq[i];
	cs[i].logp = &cl.smo_Z_MD[i]; 
	cs[i].eps = eps_bw;
	cs[i].max_iter = max_iter_bw;
      }


#if POUT
      /* sequential version */
      for (i = 0; i < cl.smo_number; i++) {	
	printf("SMO %d\n", i);
	if (!smo_changed[i]) continue;
	if (cs[i].sqd == NULL)
	  mes(MES_WIN, "cluster %d empty, no reestimate!\n", i);
	else if (sreestimate_baum_welch(&cs[i]) == -1) {
	  str = mprintf(NULL,0,"%d.reestimate false, smo[%d]\n", iter, i); 
	  mes_prot(str); m_free(str);
	  /* model_print(stdout, cl.mo[i]); */
	  goto STOP; 
	}
      }
#else
      /* parallel version */     
      i = 0;
      while (i < cl.smo_number) {
	for (j = i; j < m_min(cl.smo_number, i+THREADS); j++) {
	  if (!smo_changed[j]) continue;
	  if (cs[j].sqd == NULL)
	    mes(MES_WIN, "cluster %d empty, no reestimate!\n", j);
	  else if (perror = pthread_create(&tid[j], &Attribute, 
					   sreestimate_baum_welch, &cs[j])) {
	    str = mprintf(NULL,0, 
			  "pthread_create returns false (step %d, smo[%d])\n",
			  iter, j); 
	    mes_prot(str); m_free(str);
	    goto STOP; 
	  }
	}
	for (j = i; j < m_min(cl.smo_number, i+THREADS); j++) {
	  if (!smo_changed[j]) continue;
	  if (cs[j].sqd != NULL) {
	    pthread_join(tid[j], &return_value[j]);
	    if (return_value[j] == -1) {
	      str = mprintf(NULL,0,"%d.reestimate false, smo[%d]\n", iter, j); 
	      mes_prot(str); m_free(str);
	      goto STOP; 
	    }
	  }
	}
	i = j;
      }
#endif

      /* update model priors */
      if (CLASSIFY == 1) {
	if (sqd->total_w == 0) {
	  mes_prot("weights = 0!\n"); goto STOP;
	}
	for (i = 0; i < cl.smo_number; i++) 
	  cl.smo[i]->prior = cl.smo_seq[i]->total_w / sqd->total_w;      
      }

      fprintf(outfile, "\ntotal Prob. after %d.reestimate:\n", iter);
      scluster_print_likelihood(outfile, &cl);
    } /* if changes ... end reestimate */
  } /* while */
  
/*------------------- OUTPUT   -----------------------------------------------*/
  
  scluster_calc_print_indices(outfile, &cl);

  if (scluster_out(&cl, sqd, outfile, argv) == -1) 
    { mes_proc(); goto STOP; }


/*--------------------------------------------------------------------------*/
  res = 0;
STOP:
  pthread_attr_destroy(&Attribute);
  /* ...noch div. free! */
  if (outfile) fclose(outfile);
  return(res);
# undef CUR_PROC
}/* scluster_hmm */

/*============================================================================*/
/* Speicher fuer Sequenzen fuer jedes Modell nur einmal allocieren und
   nicht wie vorher fuer jede Sequenz mit realloc arbeiten. */
int scluster_update(scluster_t *cl, sequence_d_t *sqd) {
#define CUR_PROC "scluster_update"
  int i;
  sequence_d_t *seq_ptr;
  /* Speicher blockweise allocieren */
  for (i = 0; i < cl->smo_number; i++) {
    if (cl->smo_seq[i]) {
      /* wichtig: hier kein sequence_free, sonst sind auch die Original
	 Sequenzen futsch */
      sequence_d_clean(cl->smo_seq[i]);
      m_free(cl->smo_seq[i]);
    }
    if (cl->seq_counter[i] != 0) {
      cl->smo_seq[i] = sequence_d_calloc(cl->seq_counter[i]);
      cl->smo_seq[i]->seq_number = 0; /* wird unten hochgezaehlt */
    }
    else
      cl->smo_seq[i] = NULL;
  }
  /* Eintraege setzen */
  for (i = 0; i < sqd->seq_number; i++) {
    seq_ptr = cl->smo_seq[sqd->seq_label[i]];
    seq_ptr->seq_len[seq_ptr->seq_number] = sqd->seq_len[i];
    seq_ptr->seq[seq_ptr->seq_number] = sqd->seq[i]; /* Pointer!!! */
    seq_ptr->seq_label[seq_ptr->seq_number] = sqd->seq_label[i];
    seq_ptr->seq_w[seq_ptr->seq_number] = sqd->seq_w[i];
    seq_ptr->seq_number++;
    seq_ptr->total_w += sqd->seq_w[i];
  }
  return(0);
# undef CUR_PROC
} /* scluster_update */

/*============================================================================*/
void scluster_print_likelihood(FILE *outfile, scluster_t *cl) {
  double  total_Z_MD = 0.0, total_Z_MAW = 0.0;
  int i;
  for (i = 0; i < cl->smo_number; i++) {
    total_Z_MD += cl->smo_Z_MD[i];
    total_Z_MAW += cl->smo_Z_MAW[i];
    fprintf(outfile, "smo %d\t(#Seq. %4ld):\t", i, cl->seq_counter[i]);
    if (CLASSIFY == 0)
      fprintf(outfile, "ZMD:%8.2f", cl->smo_Z_MD[i]);
    else
      fprintf(outfile, "ZMAW:%8.2f", cl->smo_Z_MAW[i]);
    fprintf(outfile, "\n");      
  }
  fprintf(outfile, "sum\t");
  if  (CLASSIFY == 0)
    fprintf(outfile, "ZMD: %12.5f", total_Z_MD);
  else
    fprintf(outfile, "ZMAW: %12.5f", total_Z_MAW);
  fprintf(outfile, "\n\n");      
  if (CLASSIFY == 0)
    printf("total error function (ZMD): %15.4f\n", total_Z_MD);
  else 
    printf("total error function (ZMAW): %15.4f\n", total_Z_MAW);
} /* scluster_print_likelihood */

/*============================================================================*/
/* Verhindert, dass Modelle leer ausgehen, bzw nur eine
   Sequenz zugeordnet haben, 
   indem ihnen eine zufaellige Sequenz zugeordnet wird. Da hierdurch 
   evt. erneut leere Modelle erzeugt werden, werden Sequenzen getauscht,
   bis keine leeren Modelle vorliegen. (Gefahr einer Endlosschleife, daher
   Abbruch nach 100 Interationen) */
int scluster_avoid_empty_smodel(sequence_d_t *sqd, scluster_t *cl){
#define CUR_PROC "scluster_avoid_empty_smodel"
  int i, best_smo;
  long i_old, j = 0;
  char error = 1, change = 0;
  int iter = 0;
  double log_p_minus, log_p_plus;
  double *result = NULL;

  /* falls es nur ein Modell oder eine Sequenz gibt: nix zu tun */
  if (sqd->seq_number < 2 || cl->smo_number < 2)
    return(0);
  if (CLASSIFY == 1)
    if(!m_calloc(result, cl->smo_number)) {mes_proc(); goto STOP;}

  while (error && iter < 100) {
    iter++;
    error = 0;
    for (i = 0; i < cl->smo_number; i++) {
      if (cl->seq_counter[i] <= 1) {
	change = 1;
	/* zuf. Sequenz fuer leeres Modell auswaehlen */
	j = m_int(gsl_rng_uniform(RNG) * (sqd->seq_number - 1));
	/* falls Seq. von leerem Modell nicht erzeugt werden kann: abbrechen*/
	if (sfoba_logp(cl->smo[i], sqd->seq[j], sqd->seq_len[j],
		       &log_p_plus) == -1) 
	  continue;
	if (CLASSIFY == 1) {
	  best_smo = smap_bayes(cl->smo, result, cl->smo_number, sqd->seq[j],
				sqd->seq_len[j]);
	  if (best_smo == -1) continue;
	}
	i_old = sqd->seq_label[j];
	/* updates  */
	printf("!!\"avoid_empty_model\": move seq. %ld: %ld --> %d !!\n", 
	       j, i_old, i);
	cl->seq_counter[i_old]--;
	cl->seq_counter[i]++;
	sqd->seq_label[j] = i;
	/* smo_Z_MD update */
	if (sfoba_logp(cl->smo[i_old], sqd->seq[j], sqd->seq_len[j],
		       &log_p_minus) == -1) 
	  { mes_prot("sfoba_logp returns -1!\n"); goto STOP; };
	cl->smo_Z_MD[i_old] -= sqd->seq_w[j] * log_p_minus;
	cl->smo_Z_MD[i] += sqd->seq_w[j] * log_p_plus;
	/* smo_Z_MAW update */
	if (CLASSIFY == 1) {
	  cl->smo_Z_MAW[i_old] -= sqd->seq_w[j] * result[i_old];
	  cl->smo_Z_MAW[i] += sqd->seq_w[j] * result[i];
	}
      }
    }
    /* jetzt alles ok ? */
    if (change) {
      for (i = 0; i < cl->smo_number; i++) {
	if (cl->seq_counter[i] <= 1) {
	  error = 1;
	  break;
	}
      }
    }
  } /* while */
  if (!error) return (0);
STOP: 
  if (result) m_free(result);
  return(-1);
#undef CUR_PROC
} /* scluster_avoid_empty_smodel */

/*============================================================================*/
long scluster_update_label(long *oldlabel, long *seq_label, long seq_number,
			   long *smo_changed) {
  long i, changes = 0; 
  for (i = 0; i < seq_number; i++) 
    if (oldlabel[i] != seq_label[i]) {
      changes++;
      smo_changed[oldlabel[i]] = smo_changed[seq_label[i]] = 1;
      oldlabel[i] = seq_label[i];      
    }
  return changes;
} /* scluster_update_label */

/*============================================================================*/

/* scluster_best_model ermittelt aus einer vorher berechneten
   Wahrscheinlichkeitsmatrix, welches Modell am besten zur Sequenz mit der 
   ID seq_id passt. */
int scluster_best_model(scluster_t *cl, long seq_id, double **all_log_p, double *log_p) {
#define CUR_PROC "scluster_best_model"
  int i, smo_id;
  double save = -DBL_MAX;
  
  *log_p = -DBL_MAX;
  smo_id = -1;
  /* MD-Classify: argmax_i (log(p(O | \lambda_i ))) */
  if (CLASSIFY == 0) {
    for (i = 0; i < cl->smo_number; i++) {
      if (all_log_p[i][seq_id] == PENALTY_LOGP) continue; /* model and seq. don't fit */
      if (*log_p < all_log_p[i][seq_id]) {
	*log_p = all_log_p[i][seq_id];
	smo_id = i;
      }  
    }
  }
  /* MAW-Classify: argmax_i (log(p(O | \lambda_i )) + log(p( \lambda_i))) */
  else {
    for (i = 0; i < cl->smo_number; i++) {
      if (cl->smo[i]->prior == -1) {
	mes_prot("Error! Need model prior for MAW-classification\n");
      }
      if (cl->smo[i]->prior != 0) {
	if (all_log_p[i][seq_id] == PENALTY_LOGP) continue; /* model and seq. don't fit */
	if (save < all_log_p[i][seq_id] + log(cl->smo[i]->prior)) {
	  save = all_log_p[i][seq_id] + log(cl->smo[i]->prior);
	  *log_p = all_log_p[i][seq_id];
	  smo_id = i;
	} 
      }
    } 
  }


  return (smo_id);
#undef CUR_PROC
} /* scluster_best_model */

/*============================================================================*/

void scluster_prob(smosqd_t *cs) {
  int i;
  for (i = 0; i < cs->sqd->seq_number; i++)
    if (sfoba_logp(cs->smo, cs->sqd->seq[i], cs->sqd->seq_len[i], 
		   &(cs->logp[i])) == -1)
      cs->logp[i] = (double) PENALTY_LOGP;  /*  Strafkosten */
} /* scluster_prob */



/*============================================================================*/
int scluster_calc_print_indices(FILE *outfile, scluster_t *cl) {
#define CUR_PROC "scluster_calc_print_indices"
  int res = -1;
  double **x = NULL;
  double **w = NULL;
  int *len = NULL;
  double *wi = NULL;
  double logp, total_logp;
  int k, i, j_id, j_ct, j_ks, n1, n2, n3;
  //int Tk,j;
  //double *Pk_inv, *Pk, **Rjk, *Rj, Djk, DB, **D1, **D2, **S1, **S2,
  //  dangle, ddiff, s1, s2;
  double *idiff, ID_max, ID_mean, ID_smax, ID_sum, id_smax,
    *cprob, *ksprob, CT_mean, CT_smin, CT_min, ct_smin, KS_mean,
    sum1, sum2, sum3, wsum1, wsum2, wsum3;

  fprintf(outfile, "\nIndices of cluster validity:\n");
  fprintf(outfile, "# ID_mean, ID_sum, CT_mean und KS_mean weighted");
  fprintf(outfile, " over states\n");

  /* LogLikelihood */
  total_logp = 0.0;
  for (k = 0; k < cl->smo_number; k++) {
    logp = cl->smo_Z_MD[k];
    total_logp += logp;
  }
  fprintf(outfile, "Likelihood-Sum: %.4f\n", total_logp);

#if 0
  /*** Inneres Mass fuer jedes Cluster: Pk = exp(1/Tk*logP) */

  if(!m_calloc(Pk, cl->smo_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(Pk_inv, cl->smo_number)) {mes_proc(); goto STOP;}
  for (k = 0; k < cl->smo_number; k++) {
    Tk = 0; /* summierte Laengen aller Sequenzen eines Clusters */
    for (i = 0; i < cl->smo_seq[k]->seq_number; i++)
      Tk += cl->smo_seq[k]->seq_len[i];
    Pk[k] = exp(1/(double)Tk * cl->smo_Z_MD[k]);
    Pk_inv[k] = 1/Pk[k];
  }
  fprintf(outfile, "Vector pk[k]:\n");
  vector_d_print_prec(outfile, Pk, cl->smo_number, 8,6,"","","");
  fprintf(outfile, "Vector 1/pk[k]:\n");
  vector_d_print_prec(outfile, Pk_inv, cl->smo_number, 8,6,"","","");
  
  /*** Abstaende zwischen 2 Clustern: Dangle und Ddiff aus Lyngso et al. */ 

  D1 = matrix_d_alloc(cl->smo_number, cl->smo_number);
  if (!D1) {mes_proc(); goto STOP;}
  D2 = matrix_d_alloc(cl->smo_number, cl->smo_number);
  if (!D2) {mes_proc(); goto STOP;}
  S1 = matrix_d_alloc(cl->smo_number, cl->smo_number);
  if (!S1) {mes_proc(); goto STOP;}
  S2 = matrix_d_alloc(cl->smo_number, cl->smo_number);
  if (!S2) {mes_proc(); goto STOP;}

  for (j = 0; j < cl->smo_number; j++) {
    for(k = j+1; k < cl->smo_number; k++) {
      /* TODO: compactness-measure for smodels */
      smodel_measures(cl->smo[j], cl->smo[k], &dangle, &ddiff, &s1, &s2); 
      D1[j][k] = D1[k][j] = dangle;
      D2[j][k] = D2[k][j] = ddiff;
      S1[j][k] = S1[k][j] = s1;
      S2[j][k] = S2[k][j] = s2;
    }
  }
  fprintf(outfile, "Matrix D_angle[j][k]:\n");
  matrix_d_print_prec(outfile, D1, cl->smo_number, cl->smo_number,8,6,"","","");
  fprintf(outfile, "Matrix D_diff[j][k]:\n");
  matrix_d_print_prec(outfile, D2, cl->smo_number, cl->smo_number,8,6,"","","");
  fprintf(outfile, "Matrix S1[j][k]:\n");
  matrix_d_print_prec(outfile, S1, cl->smo_number, cl->smo_number,8,6,"","","");
  fprintf(outfile, "Matrix S2[j][k]:\n");
  matrix_d_print_prec(outfile, S2, cl->smo_number, cl->smo_number,8,6,"","","");


  /*** DB mit D_angle als Abstand zweier Cluster */
  
  Rjk = matrix_d_alloc(cl->smo_number, cl->smo_number);
  if (!Rjk) {mes_proc(); goto STOP;}
  if(!m_calloc(Rj, cl->smo_number)) {mes_proc(); goto STOP;}
  for (j = 0; j < cl->smo_number; j++) {
    Rj[j] = 0.0;
    for(k = 0; k <= j; k++) {
      if (Rjk[k][j] > Rj[j])
	Rj[j] = Rjk[k][j];
    }
    for(k = j+1; k < cl->smo_number; k++) {
      //Djk = smodel_prob_distance(cl->smo[j], cl->smo[k], Tdist, 1, 0);
      Djk = D1[j][k];
      if (Djk != -1)
	Rjk[j][k] = Rjk[k][j] = (Pk_inv[j] + Pk_inv[k]) / Djk; 
      if (Rjk[j][k] > Rj[j])
	Rj[j] = Rjk[j][k];
    }
  }
  DB = 0.0;
  for (j = 0; j < cl->smo_number; j++)
    DB += Rj[j];
  DB /= (double)cl->smo_number;

  fprintf(outfile, "Matrix R[j][k]:\n");
  matrix_d_print_prec(outfile, Rjk, cl->smo_number,cl->smo_number,8,2,"","","");
  fprintf(outfile, "Vector R[j]:\n");
  vector_d_print_prec(outfile, Rj, cl->smo_number,8,2,"","",""); 
  fprintf(outfile, "Davies-Bouldin-Verschnitt(D_angle): %.2f\n", DB);
#endif


  /*** Integral-Differenzen der Data-/analyt.Verteilungen pro state */
  ID_max = 0.0;  /* Max ueber alle models und states */
  ID_smax = 0.0; /* Max aller states EINES models, dann mitteln ueber models */ 
  ID_mean = 0.0; /* Mittel ueber alle models und states */
  ID_sum = 0.0; /* Summe ueber alle models und states */
  /*** Chi-Square-Test-Probabilities pro state */
  CT_min = DBL_MAX;  /* Min ueber alle models und states */
  CT_smin = 0.0; /* Min aller states EINES models, dann mitteln ueber models */ 
  CT_mean = 0.0; /* Mittel ueber alle models und states */
  /*** K-S-Test-Probabilities pro state */
  KS_mean = 0.0; /* Mittel ueber alle models und states */
  j_id = j_ct = j_ks = 0;
  for (k = 0; k < cl->smo_number; k++) {
    if (cl->smo_seq[k] != NULL) { /* moeglich wenn nur Zuordnung (labels -1) */
      if (smodel_statedatacdf(cl->smo[k],cl->smo_seq[k],&x,&w,&len,&wi) == -1)
	{mes_proc(); goto STOP;}
      if (smodel_idiff_index(cl->smo[k],cl->smo_seq[k],x,w,len,wi,&idiff)==-1)
	{mes_proc(); goto STOP;}
      if (smodel_chisquare_index(cl->smo[k], cl->smo_seq[k], x, w, len, wi,
				 &cprob) == -1)
	{mes_proc(); goto STOP;}
      if (smodel_kstest_index(cl->smo[k], cl->smo_seq[k], x, w, len, wi,
			      &ksprob) == -1)
	{mes_proc(); goto STOP;}
      id_smax = 0.0;
      ct_smin = DBL_MAX;      
      sum1 = sum2 = sum3 = wsum1 = wsum2 = wsum3 = 0.0;
      n1 = n2 = n3 = 0;
      for (i = 0; i < cl->smo[k]->N - 1; i++) { // end state omitted
	if (idiff[i] != -1.0) {
	  sum1 += idiff[i] * wi[i];
	  wsum1 += wi[i];
	  n1++;
	  /* Maxima jeweils ungewichtet nehmen */
	  if (idiff[i] > id_smax) id_smax =  idiff[i]; 
	  if (idiff[i] > ID_max) ID_max = idiff[i];
	}
	if (cprob[i] != -1.0) {
	  sum2 += cprob[i] * wi[i];
	  wsum2 += wi[i];
	  n2++;	  
	  /* Minima jeweils ungewichtet nehmen */
	  if (cprob[i] < ct_smin) ct_smin = cprob[i]; 
	  if (cprob[i] < CT_min) CT_min = cprob[i];
	}
	if (ksprob[i] != -1.0) {
	  sum3 += ksprob[i] * wi[i];
	  wsum3 += wi[i];
	  n3++;	  
	}
      }
      if (wsum1 > DBL_MIN) {
	ID_sum += sum1/wsum1 * n1; /* gewichtete Summe ueber die Zustaende */ 
	j_id += n1;
	ID_smax += id_smax; 
      }
      if (wsum2 > DBL_MIN) {
	CT_mean += sum2/wsum2 * n2; /* gewichtete Summe ueber die Zustaende */ 
	j_ct += n2;
	CT_smin += ct_smin; 
      }
      if (wsum3 > DBL_MIN) {
	KS_mean += sum3/wsum3 * n3; /* gewichtete Summe ueber die Zustaende */ 
	j_ks += n3;
      }
      /* free all fields and vectors */
      for (i = 0; i < cl->smo[k]->N; i++) { m_free(x[i]); m_free(w[i]); }
      m_free(x); 
      m_free(w); 
      m_free(len); 
      m_free(wi); 
      m_free(idiff);
      m_free(cprob);
      m_free(ksprob);
    }
  }
  ID_mean = ID_sum;
  ID_mean /= j_id;
  ID_smax /= cl->smo_number; 
  CT_mean /= j_ct;
  CT_smin /= cl->smo_number;
  KS_mean /= j_ks;
  fprintf(outfile, "ID_mean: %.2f\n", ID_mean);
  fprintf(outfile, "ID_smax: %.2f\n", ID_smax);
  fprintf(outfile, "ID_max:  %.2f\n", ID_max);
  fprintf(outfile, "ID_sum:  %.2f\n", ID_sum);
  fprintf(outfile, "CT_mean: %.6f\n", CT_mean);
  fprintf(outfile, "CT_smin: %.6f\n", CT_smin);
  fprintf(outfile, "CT_min:  %.6f\n", CT_min);
  fprintf(outfile, "KS_mean: %.6f\n", KS_mean);

  res = 0;
STOP:
  return res;
# undef CUR_PROC
} /* scluster_calc_print_indices */


/*============================================================================*/
int scluster_random_labels(sequence_d_t *sqd, int smo_number) {
#define CUR_PROC "scluster_random_labels"
  int j, label;
  for (j = 0; j < sqd->seq_number; j++) {
    label = m_int(gsl_rng_uniform(RNG) * (smo_number - 1));
    sqd->seq_label[j] = label;
  }
  return(0);
# undef CUR_PROC
} /* scluster_random_labels */

/*============================================================================*/

int  scluster_log_aposteriori(scluster_t *cl, sequence_d_t *sqd, int seq_id, 
			      double *log_apo) {
#define CUR_PROC "scluster_log_aposteriori"
  double *result = NULL;
  int best_smo = -1;
  if(!m_calloc(result, cl->smo_number)) {mes_proc(); goto STOP;}
  best_smo = smap_bayes(cl->smo, result, cl->smo_number, sqd->seq[seq_id],
			sqd->seq_len[seq_id]);
  if (best_smo == -1) {
    mes_prot("best_smo == -1 !\n");
    goto STOP;
  } 
  else 
    *log_apo = result[best_smo]; 
  /* TEST */
  /*  printf("result: ");
      for (i = 0; i< cl->smo_number; i++)
      printf("%.2f ", result[i]); printf("\n");
  */
STOP:
  m_free(result);
  return best_smo;
# undef CUR_PROC
} /* scluster_log_aposteriori */

/*============================================================================*/

void scluster_print_header(FILE *file, char* argv[]) {
  time_t zeit;
  time(&zeit);
  fprintf(file, "# Date: %s# scluster:\n", ctime((const time_t*)&zeit));
  fprintf(file, "# Sequence File: %s\n", argv[1]);
  fprintf(file, "# Model File: %s\n", argv[2]);
  fprintf(file, "# Start Partion Label: ");
  switch(atoi(argv[4])) {
  case 0: fprintf(file, "SP_BEST (best model)\n\n"); break;
  case 1: fprintf(file, "NO_SP (no start partition)\n\n"); break;
  case 2: fprintf(file, "SP_KM (partition from k-means)\n\n"); break;
  case 3: fprintf(file, "SP_ZUF (random start partition)\n\n"); break;
  default: fprintf(file, "???\n\n");
  }
}


#undef POUT
#undef THREADS

