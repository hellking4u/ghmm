/*******************************************************************************
  author       : Bernd Wichern and extended by Andrea Weisse and Utz J. Pape
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm/sequence.c
  created      : TIME: 11:29:02     DATE: Thu 12. February 1998
  $Id$

__copyright__

*****************************************************************************/

#include "mprintf.h"
#include "mes.h"
#include "sequence.h"
#include <math.h>
#include <float.h>
#include "mprintf.h"
#include "mes.h"
#include "matrix.h"
#include "vector.h"
#include "const.h"
#include "model.h"
#include "foba.h"
#include "sfoba.h"
#include "vector.h"
#include "rng.h"
#include "string.h"

/*============================================================================*/
sequence_t** sequence_read(const char *filename, int *sq_number) {
#define CUR_PROC "sequence_read"
  int i;
  sequence_t **sequence = NULL;
  scanner_t *s = NULL;
  *sq_number = 0;
  s = scanner_alloc(filename);  if(!s) {mes_proc(); goto STOP;}

  while(!s->err && !s->eof && s->c - '}') {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP; 
    /* sequence file */
    if (!strcmp(s->id, "SEQ")) {
      (*sq_number)++;
      /* more mem */	 
      if (m_realloc(sequence, *sq_number)) { mes_proc(); goto STOP; }
      sequence[*sq_number-1] = sequence_read_alloc(s);
      if (!sequence[*sq_number-1]) { mes_proc(); goto STOP; }
    }
    else {
      scanner_error(s, "unknown identifier"); goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  }

  scanner_free(&s);
  return sequence;

STOP:
  for (i = 0; i < *sq_number; i++) 
    sequence_free(&(sequence[i]));
  m_free(sequence);
  *sq_number = 0;
  return NULL;
#undef CUR_PROC
}

/*============================================================================*/
sequence_t* sequence_read_alloc(scanner_t *s){
#define CUR_PROC "sequence_read_alloc"
  int symbols = 0, lexWord = 0;
  sequence_t *sq = NULL;
  int seq_len_lex = 0;
  if( !m_calloc( sq, 1 ) ) {mes_proc();goto STOP;} 
  scanner_consume( s, '{' ); if(s->err) goto STOP; 
  while(!s->err && !s->eof && s->c - '}') {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP;
    /* array of sequences to read */
    if (!strcmp(s->id, "O")) {
      scanner_consume(s, '{'); if (s->err) goto STOP;
      sq->seq_number = 0;
      sq->total_w = 0.0;
      while( !s->eof && !s->err && s->c - '}') {
	/* another sequence --> realloc */
	if (m_realloc(sq->seq, sq->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sq->seq_len, sq->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sq->seq_label, sq->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sq->seq_id, sq->seq_number + 1))
	  { mes_proc(); goto STOP; }
        if (m_realloc(sq->seq_w, sq->seq_number + 1))
          { mes_proc(); goto STOP; }
	/* Label and ID */
	/* default */
	sq->seq_label[sq->seq_number] = -1;
	sq->seq_id[sq->seq_number] = -1.0;
	sq->seq_w[sq->seq_number] = 1;
	while (s->c == '<' || s->c == '(' || s->c == '|') {
	  if (s->c == '<') {
	    scanner_consume(s, '<'); if (s->err) goto STOP;	  
	    sq->seq_label[sq->seq_number] = scanner_get_int(s); 
	    if(s->err) goto STOP; 
	    scanner_consume(s, '>'); if (s->err) goto STOP;
	  }
	  if (s->c == '(') {
	    scanner_consume(s, '('); if (s->err) goto STOP;	  
	    sq->seq_id[sq->seq_number] = scanner_get_edouble(s); 
	    if(s->err) goto STOP; 
	    scanner_consume(s, ')'); if (s->err) goto STOP;
	  }
	  if (s->c == '|') {
	    scanner_consume(s, '|'); if (s->err) goto STOP;	  
	    sq->seq_w[sq->seq_number] = (double) scanner_get_int(s); 
	    if (sq->seq_w[sq->seq_number] <= 0) {
	      scanner_error(s, "sequence weight not positiv\n");
	      goto STOP;
	    }
	    if(s->err) goto STOP; 
	    scanner_consume(s, '|'); if (s->err) goto STOP;
	  }
	}

	sq->seq[sq->seq_number] = 
	  scanner_get_int_array(s, sq->seq_len+sq->seq_number);
	if (sq->seq_len[sq->seq_number] > MAX_SEQ_LEN) {
	  scanner_error(s, "sequence too long"); goto STOP;
	}
	scanner_consume(s, ';'); if (s->err) goto STOP;
	sq->total_w += sq->seq_w[sq->seq_number];
	sq->seq_number++;
      } /* while( !s->eof...) */
      if ((sq->seq_number == 0) || (sq->seq_number > MAX_SEQ_NUMBER)) {
	char *str = mprintf(NULL, 0,
			    "Number of sequences %ld exceeds possible range", 
			    sq->seq_number); 
	mes_prot(str);
	m_free(str);
	goto STOP;
      }
      scanner_consume(s, '}'); if (s->err) goto STOP;	    
    }  
    /* all possible seqs., sorted lexicographical */
    else if (!strcmp(s->id, "L")) {
      lexWord = 1;
      scanner_consume( s, '{' ); if(s->err) goto STOP; 
      while(!s->err && !s->eof && s->c - '}') {
	scanner_get_name(s);
	scanner_consume(s, '='); if(s->err) goto STOP; 
	if (!strcmp(s->id, "seq_len")) {   
	  seq_len_lex = scanner_get_int(s); if(s->err) goto STOP; 
	  if (seq_len_lex <= 0) {
	    mes_prot("Value for sequence length not allowed"); goto STOP;
	  }
	}
	else if (!strcmp(s->id, "symb")) {   
	  if (symbols < 0) {
	    mes_prot("Value for number of symbols not allowed"); goto STOP;
	  }
	  symbols = scanner_get_int(s); if(s->err) goto STOP; 
	}
	else {
	  scanner_error(s, "unknown identifier"); goto STOP;
	}	
	scanner_consume(s, ';'); if(s->err) goto STOP;
      }
      scanner_consume(s, '}');
      if ((seq_len_lex <= 0) || (symbols < 0)) {
	mes_prot("Values for seq. length or number of symbols not spezified");
	goto STOP;
      }
      sq = sequence_lexWords(seq_len_lex, symbols); 
      if (!sq) goto STOP;
    }  /*if (!strcmp(s->id, "L"))*/
    else {
      scanner_error(s, "unknown identifier"); goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  } /* while(!s->err && !s->eof && s->c - '}') */
  scanner_consume( s, '}'); if(s->err) goto STOP; 
  return(sq);  
STOP:
  sequence_free(&sq);
  return(NULL);
#undef CUR_PROC
} /* sequence_read_alloc */ 

/*============================================================================*/

sequence_d_t** sequence_d_read(const char *filename, int *sqd_number) {
#define CUR_PROC "sequence_d_read"
  int i;
  scanner_t *s = NULL;
  sequence_d_t **sequence = NULL;
  *sqd_number = 0;
  s = scanner_alloc(filename);  if(!s) {mes_proc(); goto STOP;}
  
  while(!s->err && !s->eof && s->c - '}') {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP; 
    /* sequence file */
	if (!strcmp(s->id, "SEQD")) {
      (*sqd_number)++;
      /* more mem */	 
      if (m_realloc(sequence, *sqd_number)) { mes_proc(); goto STOP; }
	  sequence[*sqd_number-1] = sequence_d_read_alloc(s);
      if (!sequence[*sqd_number-1]) { mes_proc(); goto STOP; }
    }
    else {
      scanner_error(s, "unknown identifier"); goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  }
  scanner_free(&s);
  return sequence;
STOP:
  scanner_free(&s);
  for (i = 0; i < *sqd_number; i++) 
    sequence_d_free(&(sequence[i]));
  m_free(sequence);
  *sqd_number = 0;
  return NULL;
#undef CUR_PROC
} /* sequence_d_read */ 

/*============================================================================*/
sequence_d_t* sequence_d_read_alloc(scanner_t *s){
#define CUR_PROC "sequence_d_read_alloc"
  sequence_d_t *sqd = NULL;
  if( !m_calloc( sqd, 1 ) ) {mes_proc();goto STOP;} 
  scanner_consume( s, '{' ); if(s->err) goto STOP; 
  while(!s->err && !s->eof && s->c - '}') {
    scanner_get_name(s);
    scanner_consume(s, '='); if(s->err) goto STOP;
    /* array of sequences to read */
    if (!strcmp(s->id, "O")) {
      scanner_consume(s, '{'); if (s->err) goto STOP;
      sqd->seq_number = 0;
      sqd->total_w = 0.0;
      while( !s->eof && !s->err && s->c - '}') {
	/* another sequence --> realloc */
	if (m_realloc(sqd->seq, sqd->seq_number+1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sqd->seq_len, sqd->seq_number+1))
	  { mes_proc(); goto STOP; }
	if (m_realloc(sqd->seq_label, sqd->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sqd->seq_id, sqd->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	if (m_realloc(sqd->seq_w, sqd->seq_number + 1)) 
	  { mes_proc(); goto STOP; }
	/* Label and ID and weight */
	/* default */
	sqd->seq_label[sqd->seq_number] = -1;
	sqd->seq_id[sqd->seq_number] = -1.0;
	sqd->seq_w[sqd->seq_number] = 1;
	while (s->c == '<' || s->c == '(' || s->c == '|') {
	  if (s->c == '<') {
	    scanner_consume(s, '<'); if (s->err) goto STOP;	  
	    sqd->seq_label[sqd->seq_number] = scanner_get_int(s); 
	    if(s->err) goto STOP; 
	    scanner_consume(s, '>'); if (s->err) goto STOP;
	  }
	  if (s->c == '(') {
	    scanner_consume(s, '('); if (s->err) goto STOP;	  
	    sqd->seq_id[sqd->seq_number] = scanner_get_edouble(s); 
	    if(s->err) goto STOP; 
	    scanner_consume(s, ')'); if (s->err) goto STOP;
	  }
	  if (s->c == '|') {
	    scanner_consume(s, '|'); if (s->err) goto STOP;	  
	    sqd->seq_w[sqd->seq_number] = (double) scanner_get_int(s);
	    if (sqd->seq_w[sqd->seq_number] < 0) {
	      scanner_error(s, "negativ sequence weight\n");
	      goto STOP;
	    }
	    if(s->err) goto STOP; 
	    scanner_consume(s, '|'); if (s->err) goto STOP;
	  }
	}
	sqd->seq[sqd->seq_number] = 
	  scanner_get_double_earray(s, sqd->seq_len + sqd->seq_number);
	if (sqd->seq_len[sqd->seq_number] > MAX_SEQ_LEN) {
	  scanner_error(s, "sequence too long"); goto STOP;
	}
	scanner_consume(s, ';'); if (s->err) goto STOP;
	sqd->total_w += sqd->seq_w[sqd->seq_number];
	sqd->seq_number++;
      } /* while( !s->eof...) */
      if ((sqd->seq_number == 0) || (sqd->seq_number > MAX_SEQ_NUMBER)) {
	char *str = mprintf(NULL, 0,
			    "Number of sequences %ld exceeds possible range", 
			    sqd->seq_number); 
	mes_prot(str);
	m_free(str);
	goto STOP;
      }
      scanner_consume(s, '}'); if (s->err) goto STOP;	    
    }  
    else {
      scanner_error(s, "unknown identifier"); goto STOP;
    }
    scanner_consume( s, ';' ); if(s->err) goto STOP;
  } /* while(!s->err && !s->eof && s->c - '}') */
  scanner_consume( s, '}'); if(s->err) goto STOP; 
  return(sqd);  
STOP:
  sequence_d_free(&sqd);
  return(NULL);
#undef CUR_PROC
} /* sequence_d_read_alloc */ 

/*============================================================================*/

/* Truncate Sequences in a given sqd_field; useful for Testing;
   returns truncated sqd_field; 
   trunc_ratio 0: no truncation
   trunc_ratio 1: truncation (mean truncation faktor = 0.5)
   trunc_ratio -1: 100 % truncation
*/

sequence_d_t **sequence_d_truncate(sequence_d_t **sqd_in, int sqd_fields, 
				  double trunc_ratio, int seed) {  
#define CUR_PROC "sequence_d_truncate"
  sequence_d_t **sq;
  int i, j, trunc_len;
  /* Hack, use -1 for complete truncation */
  if ((0 > trunc_ratio || 1 < trunc_ratio) && trunc_ratio != -1) {
    mes_prot("Error: trunc_ratio not valid\n"); goto STOP;
  }
  if( !m_calloc(sq, sqd_fields) ) {mes_proc();goto STOP;}

  gsl_rng_init();
  gsl_rng_set(RNG,seed);

  for (i = 0; i < sqd_fields; i++) {
    sq[i] = sequence_d_calloc(sqd_in[i]->seq_number);
    sq[i]->total_w = sqd_in[i]->total_w;
    for (j = 0; j < sqd_in[i]->seq_number; j++) {
      if(!m_calloc(sq[i]->seq[j], sqd_in[i]->seq_len[j])) {
	mes_proc(); goto STOP;
      }
      /* length of truncated seq. */
      if (trunc_ratio == -1) 
	trunc_len = 0;
      else
	trunc_len = (int) ceil(( sqd_in[i]->seq_len[j] * 
				 (1 - trunc_ratio * gsl_rng_uniform(RNG)) ));
      sequence_d_copy(sq[i]->seq[j], sqd_in[i]->seq[j], trunc_len);
      if(m_realloc(sq[i]->seq[j], trunc_len)) {mes_proc(); goto STOP;}
      sq[i]->seq_len[j] = trunc_len;
      sq[i]->seq_label[j] = sqd_in[i]->seq_label[j];
      sq[i]->seq_id[j] = sqd_in[i]->seq_id[j];
      sq[i]->seq_w[j] = sqd_in[i]->seq_w[j];
    }
  } /* for all sqd_fields */
  
  return sq;
STOP:
  return NULL;
#undef CUR_PROC
}

/*============================================================================*/
sequence_d_t *sequence_d_calloc(long seq_number) {
#define CUR_PROC "sequence_dcalloc"
  int i;
  //printf("*** sequence_d_t *sequence_d_calloc, nr: %d\n",seq_number);
  sequence_d_t *sqd = NULL;
  if (seq_number > MAX_SEQ_NUMBER) {
    char *str = mprintf(NULL, 0,
			"Number of sequences %ld exceeds possible range", 
			seq_number); 
    mes_prot(str);
    m_free(str);
    goto STOP;
  }
  if(!m_calloc(sqd, 1)) {mes_proc(); goto STOP;}
  if(!m_calloc(sqd->seq, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sqd->seq_len, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sqd->seq_label, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sqd->seq_id, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sqd->seq_w, seq_number)) {mes_proc(); goto STOP;}
  sqd->seq_number = seq_number;
  
  sqd-> total_w = 0.0;
  for (i = 0; i < seq_number; i++) {
    sqd->seq_label[i] = -1;
    sqd->seq_id[i] = -1.0;
    sqd->seq_w[i] = 1;    
  }
    
  return sqd;
 STOP:
  sequence_d_free(&sqd);
  return NULL;
#undef CUR_PROC
} /* sequence_d_calloc */

/*============================================================================*/
sequence_t *sequence_calloc(long seq_number) {
#define CUR_PROC "sequence_calloc"
  int i;
  sequence_t *sq = NULL;
  if (seq_number > MAX_SEQ_NUMBER) {
    char *str = mprintf(NULL, 0,
			"Number of sequences %ld exceeds possible range", 
			seq_number); 
    mes_prot(str);
    m_free(str);
    goto STOP;
  }
  if(!m_calloc(sq, 1)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->seq, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->states, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->seq_len, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->seq_label, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->seq_id, seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(sq->seq_w, seq_number)) {mes_proc(); goto STOP;}
  sq->seq_number = seq_number;
  for (i = 0; i < seq_number; i++) {
    sq->seq_label[i] = -1;
    sq->seq_id[i] = -1.0;
    sq->seq_w[i] = 1;
  }
  return sq;
 STOP:
  sequence_free(&sq);
  return NULL;
#undef CUR_PROC
} /* sequence_calloc */


/*============================================================================*/
sequence_t *sequence_lexWords(int n, int M) {
# define CUR_PROC "sequence_lexWords"
  
  sequence_t *sq = NULL;
  long seq_number, cnt = 0;
  int j = n - 1;
  int i;
  int *seq;
  if ((n < 0) || (M <= 0)) { mes_proc(); goto STOP; }
  seq_number = (long) pow((double) M, (double) n);
  sq = sequence_calloc(seq_number);
  if (!sq) { mes_proc(); goto STOP; }
  for (i = 0; i < seq_number; i++) {
    if (!m_calloc(sq->seq[i], n)) { mes_proc(); goto STOP; }
    sq->seq_len[i] = n;
    sq->seq_id[i] = i;
  }

  if (!m_calloc(seq, n)) { mes_proc(); goto STOP; }
  while (!(j < 0)) {
    sequence_copy(sq->seq[cnt], seq, n);
    j = n - 1;
    while (seq[j] == M - 1) {
      seq[j] = 0;
      j--;
    }
    seq[j]++;
    cnt++;
  }
  m_free(seq);
  return sq;
STOP:
  sequence_free(&sq);
  return NULL;
# undef CUR_PROC
} /* sequence_lewWords */

/*============================================================================*/
int sequence_max_symbol(sequence_t *sq) {
  long i, j;
  int max_symb = -1;
  for (i = 0; i < sq->seq_number; i++)
    for (j = 0; j < sq->seq_len[i]; j++) {
      if (sq->seq[i][j] > max_symb)
	max_symb = sq->seq[i][j];
    }
  return max_symb;
} /* sequence_max_symbol */

/*============================================================================*/
void sequence_copy(int *target, int *source, int len) {
  int i;
  for (i = 0; i < len; i ++)
    target[i] = source[i];
} /* sequence_copy */


/*============================================================================*/
void sequence_d_copy(double *target, double *source, int len) {
  int i;
  for (i = 0; i < len; i ++)
    target[i] = source[i];
} /* sequence_copy */


/*============================================================================*/
int sequence_add(sequence_t *target, sequence_t *source) {
#define CUR_PROC "sequence_add"

  int res = -1;
  int **old_seq       = target->seq;
  int **old_seq_st    = target->states;
  int *old_seq_len    = target->seq_len;
  long *old_seq_label = target->seq_label;
  double *old_seq_id = target->seq_id;
  double *old_seq_w = target->seq_w;
  long old_seq_number = target->seq_number;
  long i;

  target->seq_number = old_seq_number + source->seq_number;
  target->total_w += source->total_w;

  if(!m_calloc(target->seq, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->states, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_len, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_label, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_id, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_w, target->seq_number)) {mes_proc(); goto STOP;}

  for (i = 0; i < old_seq_number; i++) {
    target->seq[i] = old_seq[i];
    target->states[i] = old_seq_st[i];
    target->seq_len[i] = old_seq_len[i];
    target->seq_label[i] = old_seq_label[i];
    target->seq_id[i] = old_seq_id[i];
    target->seq_w[i] = old_seq_w[i];
  }

  for (i = 0; i < (target->seq_number - old_seq_number); i++) {
    if(!m_calloc(target->seq[i+old_seq_number], source->seq_len[i])) 
      {mes_proc(); goto STOP;}

    sequence_copy(target->seq[i+old_seq_number], source->seq[i], 
		  source->seq_len[i]);

    if(!m_calloc(target->states[i+old_seq_number], source->seq_len[i])) 
      {mes_proc(); goto STOP;}

    /*sequence_copy(target->states[i+old_seq_number], source->states[i], 
      source->seq_len[i]);*/

    target->seq_len[i+old_seq_number] = source->seq_len[i];
    target->seq_label[i+old_seq_number] = source->seq_label[i];
    target->seq_id[i+old_seq_number] = source->seq_id[i];
    target->seq_w[i+old_seq_number] = source->seq_w[i];
  } 
 
  /*m_free(old_seq);
    m_free(old_seq_st);*/
  m_free(old_seq_len);
  m_free(old_seq_label);
  m_free(old_seq_id);
  m_free(old_seq_w);
  res = 0;
 STOP:
  return res;
#undef CUR_PROC
}


/*============================================================================*/
int sequence_d_add(sequence_d_t *target, sequence_d_t *source) {
#define CUR_PROC "sequence_d_add"

  int res = -1;
  double **old_seq    = target->seq;
  int *old_seq_len    = target->seq_len;
  long *old_seq_label = target->seq_label;
  double *old_seq_id  = target->seq_id;
  double *old_seq_w  = target->seq_w;
  long old_seq_number = target->seq_number;  
  long i;

  target->seq_number = old_seq_number + source->seq_number;
  target->total_w += source->total_w;

  if(!m_calloc(target->seq, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_len, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_label, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_id, target->seq_number)) {mes_proc(); goto STOP;}
  if(!m_calloc(target->seq_w, target->seq_number)) {mes_proc(); goto STOP;}

  for (i = 0; i < old_seq_number; i++) {
    target->seq[i] = old_seq[i];
    target->seq_len[i] = old_seq_len[i];
    target->seq_label[i] = old_seq_label[i];
    target->seq_id[i] = old_seq_id[i];
    target->seq_w[i] = old_seq_w[i];
  }

  for (i = 0; i < (target->seq_number - old_seq_number); i++) {
    if(!m_calloc(target->seq[i+old_seq_number], source->seq_len[i])) 
      {mes_proc(); goto STOP;}

    sequence_d_copy(target->seq[i+old_seq_number], source->seq[i], 
		    source->seq_len[i]);
    target->seq_len[i+old_seq_number] = source->seq_len[i];
    target->seq_label[i+old_seq_number] = source->seq_label[i];
    target->seq_id[i+old_seq_number] = source->seq_id[i];
    target->seq_w[i+old_seq_number] = source->seq_w[i];
  } 
 
  m_free(old_seq);
  m_free(old_seq_len);
  m_free(old_seq_label);
  m_free(old_seq_id);
  m_free(old_seq_w);
  res = 0;

 STOP:
  return res;
#undef CUR_PROC
}

/*============================================================================*/
int sequence_check(sequence_t *sq, int max_symb) {
#define CUR_PROC "sequence_check"
  int i, j;
  for (j = 0; j < sq->seq_number; j++) {
    for (i = 0; i < sq->seq_len[j]; i++) {
      if ((sq->seq[j][i] >= max_symb) || (sq->seq[j][i] < 0)) {
	char *str =  
	  mprintf(NULL, 0, "Wrong symbol \'%d\' in sequence %d at Pos. %d;\
                            Should be within [0..%d]\n",
		  sq->seq[j][i], j + 1, i + 1, max_symb - 1);
	mes_prot(str);	
	m_free(str); return (-1);
      }
    }
  }
  return 0;
#undef CUR_PROC
} /* sequence_check */

/*============================================================================*/
int sequence_best_model(model **mo, int model_number, int *sequence, 
			int seq_len, double *log_p) {
# define CUR_PROC "seqence_best_model"
  double log_ptmp;
  int model_index, i;
  *log_p = -DBL_MAX;
  model_index = -1;
  for (i = 0; i < model_number; i++) {    
    foba_logp(mo[i], sequence, seq_len, &log_ptmp);
    if (log_ptmp != +1 && log_ptmp > *log_p) {
      *log_p = log_ptmp;
      model_index = i;
    }
  }
  if (*log_p == -DBL_MAX)
    *log_p = +1;
  return (model_index);
# undef CUR_PROC
} /* sequence_best_model */


/*============================================================================*/
void sequence_print(FILE *file, sequence_t *sq) {
  int i, j;
  fprintf(file, "SEQ = {\n\tO = {\n");
  for (i = 0; i < sq->seq_number; i++) {
    if (sq->seq_id[i] != -1.0)
      fprintf(file, "\t(%10.0f)", sq->seq_id[i]);
    if (sq->seq_label[i] != -1)
      fprintf(file, "\t<%ld>", sq->seq_label[i]);
    if (sq->seq_w[i] != 1)
      fprintf(file, "\t|%.0f|", sq->seq_w[i]);
    fprintf(file,"\t");
    if (sq->seq_len[i] > 0) {
      fprintf(file, "%d", sq->seq[i][0]);
      for (j = 1; j < sq->seq_len[i]; j++)
	fprintf(file, ", %d", sq->seq[i][j]);
    }
    fprintf(file, ";\n");
  }
  fprintf(file, "\t};\n};\n\n");
} /* sequence_print */

/*============================================================================*/
/**/
void sequence_print_xml(FILE *file, sequence_t *sq) {
  int i, j;
  /* coding missing */
  fprintf(file, "<Sequences type=\"int\" >\n");
  fprintf(file, " <DiscretePD>\n");
  for (i = 0; i < sq->seq_number; i++) {
    fprintf(file, "  %.0f <Sequence", sq->seq_w[i]);
    if (sq->seq_id[i] != -1.0)
      fprintf(file, " id=\"seq%f\" ", sq->seq_id[i]);
    fprintf(file,">");
    if (sq->seq_label[i] != -1)
      fprintf(file, "<Label>%ld</Label>", sq->seq_label[i]);
    if (sq->seq_len[i] > 0) {
      fprintf(file, "<!-- Length: %d -->", sq->seq_len[i]);
      for (j = 0; j < sq->seq_len[i]; j++)
	fprintf(file, " %d", sq->seq[i][j]);
    }
    fprintf(file, "  </Sequence>\n");
  }
  fprintf(file, " </DiscretePD>\n");
  fprintf(file, "</Sequences>\n");
} /* sequence_print_xml */

void sequence_d_print_xml(FILE *file, sequence_d_t *sq) {
  int i, j;
  /* coding missing */
  fprintf(file, "<Sequences type=\"int\" >\n");
  fprintf(file, " <DiscretePD>\n");
  for (i = 0; i < sq->seq_number; i++) {
    fprintf(file, "  %.0f <Sequence", sq->seq_w[i]);
    if (sq->seq_id[i] != -1.0)
      fprintf(file, " id=\"seq%f\" ", sq->seq_id[i]);
    fprintf(file,">");
    if (sq->seq_label[i] != -1)
      fprintf(file, "<Label>%ld</Label>", sq->seq_label[i]);
    if (sq->seq_len[i] > 0) {
      fprintf(file, "<!-- Length: %d -->", sq->seq_len[i]);
      for (j = 0; j < sq->seq_len[i]; j++)
	fprintf(file, " %f", sq->seq[i][j]);
    }
    fprintf(file, "  </Sequence>\n");
  }
  fprintf(file, " </DiscretePD>\n");
  fprintf(file, "</Sequences>\n");
} /* sequence_print_xml */

/*============================================================================*/

void sequence_mathematica_print(FILE *file, sequence_t *sq, char *name) {
  int i;
  fprintf(file, "%s = {\n", name);
  for (i = 0; i < sq->seq_number - 1; i++)
    vector_i_print(file, sq->seq[i], sq->seq_len[i], 
      "{", ",", "},");
  /* no comma after last seq. */
  vector_i_print(file, sq->seq[sq->seq_number - 1], 
		 sq->seq_len[sq->seq_number - 1], 
		 "{", ",", "}"); 
  fprintf(file, "};\n");
} /* sequence_d_mathematica_print */

/*============================================================================*/

void sequence_d_gnu_print(FILE *file, sequence_d_t *sqd) {
  int i, j;
  for (i = 0; i < sqd->seq_number; i++) {
    for (j = 0; j < sqd->seq_len[i]; j++)
      fprintf(file, "%.8f\n", sqd->seq[i][j]);
    fprintf(file, "\n\n");
  }
}

/*============================================================================*/
void sequence_d_print(FILE *file, sequence_d_t *sqd, int discrete) {
  int i, j;
  fprintf(file, "SEQD = {\n\tO = {\n");
  for (i = 0; i < sqd->seq_number; i++) {
    if (sqd->seq_id[i] != -1.0)
      fprintf(file, "\t(%10.0f)", sqd->seq_id[i]);
    if (sqd->seq_label[i] != -1)
      fprintf(file, "\t<%ld>", sqd->seq_label[i]);
    if (sqd->seq_w[i] != 1)
      fprintf(file, "\t|%.0f|", sqd->seq_w[i]);
    fprintf(file,"\t");
    if (sqd->seq_len[i] > 0) {
      if (discrete)
	fprintf(file, "%3.0f", sqd->seq[i][0]);
      else {
	if (sqd->seq[i][0] > 500) 
	  fprintf(file, "%8.0f", sqd->seq[i][0]);
	else
	  fprintf(file, "%8.2f", sqd->seq[i][0]);
      }
      for (j = 1; j < sqd->seq_len[i]; j++) {
	if (discrete)
	  fprintf(file, ", %3.0f", sqd->seq[i][j]);
	else {
	  if (sqd->seq[i][j] > 500) 
	    fprintf(file, ", %8.0f", sqd->seq[i][j]);
	  else
	    fprintf(file, ", %8.2f", sqd->seq[i][j]);
	}
      }
    }
    fprintf(file, ";\n");
  }
  fprintf(file, "\t};\n};\n\n");
} /* sequence_print */

/*============================================================================*/

void sequence_d_mathematica_print(FILE *file, sequence_d_t *sqd, char *name) {
  int i;
  fprintf(file, "%s = {\n", name);
  for (i = 0; i < sqd->seq_number - 1; i++)
      vector_d_print(file, sqd->seq[i], sqd->seq_len[i], 
      "{", ",", "},");
  /* no comma after last seq. */
  vector_d_print(file, sqd->seq[sqd->seq_number - 1], 
		 sqd->seq_len[sqd->seq_number - 1], 
		 "{", ",", "}"); 
  fprintf(file, "};\n");
} /* sequence_d_mathematica_print */

/*============================================================================*/
void sequence_clean(sequence_t *sq) {
  /* keep data, only delete references */
  m_free(sq->seq); 
  m_free(sq->seq_len);
  m_free(sq->seq_label);
  m_free(sq->seq_id);
  m_free(sq->seq_w);
  
  sq->seq_number = 0;
  sq->total_w = 0.0;
} /* sequence_clean */

/*============================================================================*/
void sequence_d_clean(sequence_d_t *sqd) {
  /* keep data, only delete references */
  m_free(sqd->seq);
  m_free(sqd->seq_len);
  m_free(sqd->seq_label);
  m_free(sqd->seq_id);
  m_free(sqd->seq_w);
  sqd->seq_number = 0;
  sqd->total_w = 0.0;
} /* sequence_d_clean */

/*============================================================================*/
int sequence_free(sequence_t **sq) {
# define CUR_PROC "sequence_free"
  mes_check_ptr(sq, return(-1));
  if( !*sq ) return(0);
  matrix_i_free(&(*sq)->seq, (*sq)->seq_number);
  /* The allocation of state must be fixed */
  /*** Added attribute to the sequence_t
  if (&(*sq)->states) { 
    matrix_i_free(&(*sq)->states, (*sq)->seq_number);
   }
  ***/
  m_free((*sq)->seq);
  m_free((*sq)->seq_len);
  m_free((*sq)->seq_label);
  m_free((*sq)->seq_id);
  m_free((*sq)->seq_w);
  m_free(*sq);
  return 0;
# undef CUR_PROC
} /* sequence_free */


/*============================================================================*/
int sequence_d_free(sequence_d_t **sqd) {
# define CUR_PROC "sequence_d_free"
  mes_check_ptr(sqd, return(-1));
  if( !*sqd ) return(0);
  
  // sequence_d_print(stdout,*sqd,0);
		  
  matrix_d_free(&(*sqd)->seq, (*sqd)->seq_number);
  m_free((*sqd)->seq_len);
  m_free((*sqd)->seq_label);
  m_free((*sqd)->seq_id);
  m_free((*sqd)->seq_w);
  m_free(*sqd);
  return 0;
# undef CUR_PROC
} /* sequence_d_free */

/*============================================================================*/
sequence_d_t *sequence_d_create_from_sq(const sequence_t *sq) {
# define CUR_PROC "sequence_d_create_from_sq"
  int j, i;
  sequence_d_t *sqd = NULL; /* target seq. array */
  if (!(sqd = sequence_d_calloc(sq->seq_number))) { mes_proc(); goto STOP; }
  for (j = 0; j < sq->seq_number; j++) {
    if (!m_calloc(sqd->seq[j], sq->seq_len[j])) { mes_proc(); goto STOP; } 
    for (i = 0; i < sq->seq_len[j]; i++)
      sqd->seq[j][i] = (double)(sq->seq[j][i]);
    sqd->seq_len[j] = sq->seq_len[j];
    sqd->seq_label[j] = sq->seq_label[j];
    sqd->seq_id[j] = sq->seq_id[j];
    sqd->seq_w[j] = sq->seq_w[j];
  }
  sqd->seq_number = sq->seq_number;
  sqd->total_w = sq->total_w;
  return(sqd);
STOP:
  sequence_d_free(&sqd);
  return NULL;
#undef CUR_PROC
} /* sequence_d_create_from_sq */ 

/*============================================================================*/
sequence_t *sequence_create_from_sqd(const sequence_d_t *sqd) {
# define CUR_PROC "sequence_create_from_sqd"
  int j, i;
  sequence_t *sq = NULL; /* target seq. array */
  if (!(sq = sequence_calloc(sqd->seq_number))) { mes_proc(); goto STOP; }
  for (j = 0; j < sqd->seq_number; j++) {
    if (!m_calloc(sq->seq[j], sqd->seq_len[j])) { mes_proc(); goto STOP; } 
    for (i = 0; i < sqd->seq_len[j]; i++) {
      sq->seq[j][i] = m_int(fabs(sqd->seq[j][i]));
    }
    sq->seq_len[j] = sqd->seq_len[j];
    sq->seq_label[j] = sqd->seq_label[j];
    sq->seq_id[j] = sqd->seq_id[j];
    sq->seq_w[j] = sqd->seq_w[j];
  }
  sq->seq_number = sqd->seq_number;
  sq->total_w = sqd->total_w;
  return(sq);
STOP:
  sequence_free(&sq);
  return NULL;
#undef CUR_PROC
} /* sequence_create_from_sqd */ 

/*============================================================================*/

int sequence_max_len(const sequence_t *sqd) {
  int i, max_len = 0;
  for (i = 0; i < sqd->seq_number; i++)
    if (max_len < sqd->seq_len[i])
      max_len = sqd->seq_len[i];
  return max_len;
} /* sequence_max_len */


/*============================================================================*/

int sequence_d_max_len(const sequence_d_t *sqd) {
  int i, max_len = 0;
  for (i = 0; i < sqd->seq_number; i++)
    if (max_len < sqd->seq_len[i])
      max_len = sqd->seq_len[i];
  return max_len;
} /* sequence_d_max_len */

/*============================================================================*/

sequence_d_t *sequence_d_mean(const sequence_d_t *sqd) {
# define CUR_PROC "sequence_d_mean"
  int i, j, max_len;
  sequence_d_t *out_sqd = NULL;

  max_len = sequence_d_max_len(sqd);
  if (!(out_sqd = sequence_d_calloc(1))) { mes_proc(); goto STOP; }
  if (!m_calloc(out_sqd->seq[0], max_len)) { mes_proc(); goto STOP; } 
  out_sqd->seq_len[0] = max_len;
  
  for (i = 0; i < sqd->seq_number; i++) 
    for (j = 0; j < sqd->seq_len[i]; j++)
      out_sqd->seq[0][j] += sqd->seq[i][j];

  for (j = 0; j < max_len; j++)
    out_sqd->seq[0][j] /= sqd->seq_number;

  return out_sqd;
STOP:
  sequence_d_free(&out_sqd);
  return NULL;
# undef CUR_PROC
} /* sequence_d_mean */

/*============================================================================*/

double **sequence_d_scatter_matrix(const sequence_d_t *sqd, int *dim) {
# define CUR_PROC "sequence_d_scatter_matrix"
  int *count, k,l,i,j;
  double **W, *mean;

  *dim = sequence_d_max_len(sqd);
  if (!(W = matrix_d_alloc(*dim, *dim))) {mes_proc(); goto STOP;}
  
  /* Mean over all sequences. Individual counts for each dimension */
  if (!m_calloc(mean, *dim)) { mes_proc(); goto STOP; } 
  if (!m_calloc(count, *dim)) { mes_proc(); goto STOP; } 
  for (i = 0; i < sqd->seq_number; i++) {
    for (l = 0; l < sqd->seq_len[i]; l++) {
      mean[l] += sqd->seq[i][l];
      count[l]++;
    }
  }
  for (l = 0; l < *dim; l++)
    mean[l] /= count[l];
  /* scatter matrix (upper triangle) */
  for (j = 0; j < sqd->seq_number; j++) {
    for (k = 0; k < *dim; k++) {
      for (l = k; l < *dim; l++) { 
	if (sqd->seq_len[j] > l) 
	  W[k][l] += (sqd->seq[j][k] - mean[k])*(sqd->seq[j][l] - mean[l]);
      }
    }
  }
  /* norm with counts, set lower triangle */
  for (k = 0; k < *dim; k++) {
    for (l = *dim-1; l >= 0; l--) { 
      if (l >= k) W[k][l] /= (double)count[l];
      else W[k][l] = W[l][k];
    }
  }
  return W;
STOP:
  matrix_d_free(&W, *dim);
  return NULL;
# undef CUR_PROC
} /* sequence_d_scatter_matrix */

/*============================================================================*/
/* dummy function at the moment */

int sequence_d_class(const double *O, int index, double *osum) {
#define CUR_PROC "sequence_d_class"
  
  return 0; 
# undef CUR_PROC
} /* sequence_d_class */

/*============================================================================*/

/* divide given field of seqs. randomly into two different fields. Also do 
   allocating. train_ratio determines approx. the fraction of seqs. that go
   into the train_set and test_set resp.
*/

int sequence_d_partition(sequence_d_t *sqd, sequence_d_t * sqd_train, 
			 sequence_d_t *sqd_test, double train_ratio) {
#define CUR_PROC "sequence_d_partition"

  double p;
  sequence_d_t *sqd_dummy = NULL;
  int i;
  long total_seqs, cur_number;

  total_seqs = sqd->seq_number;
  if (total_seqs <= 0) {
    mes_prot("Error: number of seqs. less or equal zero\n");
    goto STOP;
  }
  /* waste of memory but avoids to many reallocations */
  sqd_dummy = sqd_train;
  for (i = 0; i < 2; i++) {
    if(!m_calloc(sqd_dummy->seq, total_seqs)) {mes_proc(); goto STOP;}
    if(!m_calloc(sqd_dummy->seq_len, total_seqs)) {mes_proc(); goto STOP;}
    if(!m_calloc(sqd_dummy->seq_label, total_seqs)) {mes_proc(); goto STOP;}
    if(!m_calloc(sqd_dummy->seq_id, total_seqs)) {mes_proc(); goto STOP;}
    if(!m_calloc(sqd_dummy->seq_w, total_seqs)) {mes_proc(); goto STOP;}
    sqd_dummy->seq_number = 0;
    sqd_dummy = sqd_test;
  }

  for (i = 0; i < total_seqs; i++) {
    p = gsl_rng_uniform(RNG);
    if (p <= train_ratio)
      sqd_dummy = sqd_train;
    else
      sqd_dummy = sqd_test;
    cur_number = sqd_dummy->seq_number;
    if(!m_calloc(sqd_dummy->seq[cur_number], sqd->seq_len[i])) 
      {mes_proc(); goto STOP;}
    /* copy all entries */
    sequence_d_copy_all(sqd_dummy, cur_number, sqd, i);
    sqd_dummy->seq_number++;
  }
  
  /* reallocs*/
  sqd_dummy = sqd_train;
  for (i = 0; i < 2; i++) {
    if(m_realloc(sqd_dummy->seq, sqd_dummy->seq_number)) {mes_proc(); goto STOP;}
    if(m_realloc(sqd_dummy->seq_len, sqd_dummy->seq_number)) {mes_proc(); goto STOP;}
    if(m_realloc(sqd_dummy->seq_label, sqd_dummy->seq_number)) {mes_proc(); goto STOP;}
    if(m_realloc(sqd_dummy->seq_id, sqd_dummy->seq_number)) {mes_proc(); goto STOP;}
    if(m_realloc(sqd_dummy->seq_w, sqd_dummy->seq_number)) {mes_proc(); goto STOP;}
    sqd_dummy = sqd_test;
  }
  return 0;

 STOP:
  return -1;
#undef CUR_PROC
}

/*============================================================================*/

void sequence_d_copy_all(sequence_d_t *target, long t_num, 
			 sequence_d_t *source, long s_num) {

    sequence_d_copy(target->seq[t_num], source->seq[s_num], source->seq_len[s_num]);  
    target->seq_len[t_num] = source->seq_len[s_num];
    target->seq_label[t_num] = source->seq_label[s_num];
    target->seq_id[t_num] = source->seq_id[s_num];
    target->seq_w[t_num] = source->seq_w[s_num];
}


/*============================================================================*/
/* Likelihood function in a mixture model:
   sum_k w^k log( sum_c (alpha_c p(O^k | lambda_c)))
*/

int sequence_d_mix_like(smodel **smo, int  smo_number, sequence_d_t *sqd,
			double *like) {
#define CUR_PROC "sequence_d_mix_like"
  int i, k, error_seqs = 0;
  double seq_like = 0.0, log_p;
  
  *like = 0.0;

  for (i = 0; i < sqd->seq_number; i++) {
    seq_like = 0.0;
    for (k = 0; k < smo_number; k++) {
      if (sfoba_logp(smo[k], sqd->seq[i], sqd->seq_len[i], &log_p) != -1) {
	if (log_p > -100)
	  seq_like += exp(log_p) * smo[k]->prior; 
      }
    }
    /* no model fits */
    if (seq_like == 0.0) {
      error_seqs++;
      *like += (PENALTY_LOGP * sqd->seq_w[i]); 
    }
    else
      *like += (log(seq_like) * sqd->seq_w[i]);
  }

  return error_seqs;

#undef CUR_PROC
} /* sequence_d_mix_like */


