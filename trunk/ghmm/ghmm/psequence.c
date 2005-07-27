#include "psequence.h"
#include <stdlib.h>
#include "mes.h"
#include "matrix.h"

mysequence * init_mysequence(int length, int number_of_alphabets, int number_of_d_seqs) {
#define CUR_PROC "init_mysequence"
  int i;
  mysequence * seq;
  seq = (mysequence*)malloc(sizeof(mysequence));
  if (!seq) {mes_proc(); goto STOP;}
  seq->length = length;
  seq->number_of_alphabets = number_of_alphabets;
  seq->number_of_d_seqs = number_of_d_seqs;
  seq->seq = NULL;
  seq->d_value = NULL;
  if (number_of_alphabets > 0) 
    seq->seq = matrix_i_alloc(number_of_alphabets, length); 
  if (number_of_d_seqs > 0)
    seq->d_value = matrix_d_alloc(number_of_d_seqs, length);
  return seq;
STOP:
  free_mysequence(seq);
  return NULL;
#undef CUR_PROC
}

int free_mysequence(mysequence * seq, int number_of_alphabets, int number_of_d_seqs) {
#define CUR_PROC "free_mysequence"
  int i;
  mes_check_ptr(seq, return(-1));
  if ( seq == NULL ) return(0);
  if (seq->seq != NULL) {
    for (i=0; i<number_of_alphabets; i++) 
      m_free(seq->seq[i]);
    m_free(seq->seq);
  }
  if (seq->d_value != NULL) {
    for (i=0; i<number_of_d_seqs; i++)
      m_free(seq->d_value[i]);
    m_free(seq->d_value);
  }
  m_free(seq);
  return 0;
#undef CUR_PROC
}

void set_discrete_mysequence(mysequence * seq_pointer, int index, int * int_seq) {
  seq_pointer->seq[index] = int_seq;
}
  
void set_continuous_mysequence(mysequence * seq_pointer, int index, double * d_seq) {
  seq_pointer->d_value[index] = d_seq;
}

int * get_discrete_mysequence(mysequence * seq_pointer, int index){
  return seq_pointer->seq[index];
}

double * get_continuous_mysequence(mysequence * seq_pointer, int index){
  return seq_pointer->d_value[index];
}

mysequence * slice_mysequence(mysequence * seq_pointer, int start, int stop){
  if (stop > seq_pointer->length) {
    fprintf(stderr, "Slice: sequence index (%i) out of bounds (%i)\n", 
	    stop, seq_pointer->length);
  }
  mysequence * slice = init_mysequence(stop - start, 
				       seq_pointer->number_of_alphabets,
				       seq_pointer->number_of_d_seqs);
  int i, j;
  for (i=start; i<stop; i++){
    for (j=0; j<slice->number_of_alphabets; j++)
      slice->seq[j][i-start] = seq_pointer->seq[j][i];
    for (j=0; j<slice->number_of_d_seqs; j++)
      slice->d_value[j][i-start] = seq_pointer->d_value[j][i];
  }
  return slice;
}

int get_char_mysequence(mysequence * seq_pointer, int alphabet, int index){
  if (alphabet < seq_pointer->number_of_alphabets) {
    if (index < 0)
      return -1;
    if (index < seq_pointer->length) {
      return seq_pointer->seq[alphabet][index];
    }
    else {
      fprintf(stderr, "index (%i) larger than length (%i)!", index, seq_pointer->length);
      return -1;
    }
  }
  else {
    fprintf(stderr, "alphabet (%i) larger than number of alphabets (%i)!",
	    alphabet, seq_pointer->number_of_alphabets);
    return -1;
  }
}

double get_double_mysequence(mysequence * seq_pointer, int seq_index, int index){
  if (seq_index < seq_pointer->number_of_d_seqs) {
    if (index < 0)
      return 0;
    if (index < seq_pointer->length) {
      return seq_pointer->d_value[seq_index][index];
    }
    else {
      fprintf(stderr, "index (%i) larger than length (%i)!", index, seq_pointer->length);
      return 0;
    }
  }
  else {
    fprintf(stderr, "seq_index (%i) larger than number of seq_indexs (%i)!",
	    seq_index, seq_pointer->number_of_d_seqs);
    return 0;
  }
} 
