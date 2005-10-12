#ifndef GHMM_PSEQUENCE_H
#define GHMM_PSEQUENCE_H
struct psequence {
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

typedef struct psequence psequence;

psequence * ghmm_dpseq_init(int length, int number_of_alphabets, int number_of_d_seqs);

int ghmm_dpseq_free(psequence * seq, int number_of_alphabets, int number_of_d_seqs);

void ghmm_dpseq_set_discrete(psequence * seq_pointer, int index, int * int_seq);

int * ghmm_dpseq_get_discrete(psequence * seq_pointer, int index);

void ghmm_dpseq_set_continuous(psequence * seq_pointer, int index, double * d_seq);

double * ghmm_dpseq_get_continuous(psequence * seq_pointer, int index);

psequence * ghmm_dpseq_slice(psequence * seq_pointer, int start, int stop);

int ghmm_dpseq_get_char(psequence * seq_pointer, int alphabet, int index);

double ghmm_dpseq_get_double(psequence * seq_pointer, int seq_index, int index);
#endif
