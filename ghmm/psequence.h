#ifndef MYSEQUENCE_H
#define MYSEQUENCE_H
struct mysequence {
  /** for each alphabet in model->number_of_alphabets there is one int seq **/
  int ** seq;
  /** number of alphabets (same as in model) **/
  int number_of_alphabets;
  /** for each sequence position there are also double values (e.g) Ka **/
  double ** d_value;
  /** number of continous sequences **/
  int number_of_d_seqs;
  /** length of the sequence **/
  int length;
};

typedef struct mysequence mysequence;

mysequence * init_mysequence(int length, int number_of_alphabets, int number_of_d_seqs);

int free_mysequence(mysequence * seq, int number_of_alphabets, int number_of_d_seqs);

void set_discrete_mysequence(mysequence * seq_pointer, int index, int * int_seq);

int * get_discrete_mysequence(mysequence * seq_pointer, int index);

void set_continuous_mysequence(mysequence * seq_pointer, int index, double * d_seq);

double * get_continuous_mysequence(mysequence * seq_pointer, int index);

mysequence * slice_mysequence(mysequence * seq_pointer, int start, int stop);

int get_char_mysequence(mysequence * seq_pointer, int alphabet, int index);

double get_double_mysequence(mysequence * seq_pointer, int seq_index, int index);
#endif
