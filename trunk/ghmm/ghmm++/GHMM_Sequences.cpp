/*
 * created: 29 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ghmm/sequence.h"
#include "ghmm++/GHMM_Sequences.h"
#include "ghmm++/GHMM_Sequence.h"
#ifdef HAVE_CSTDLIB
#  include <cstdlib>
#else
#  include <stdlib.h>
#endif

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Sequences::GHMM_Sequences(XMLIO_Attributes& attrs) {
  bool found = false;

  if (attrs["type"] == "continuous") {
    found = true;
    sequence_type = GHMM_DOUBLE;
  }

  if (attrs["type"] == "discrete") {
    found = true;
    sequence_type = GHMM_INT;
  }

  if (! found) {
    fprintf(stderr,"unrecognized sequences type: type=\"%s\"\nchoose out of: continuous, discrete\n",attrs["type"].c_str());
    exit(1);
  }

  c_i_sequences = NULL;
  c_d_sequences = NULL;
  last_weight   = 1;
}


GHMM_Sequences::GHMM_Sequences(GHMM_SequenceType my_sequence_type) {
  sequence_type = my_sequence_type;

  c_i_sequences = NULL;
  c_d_sequences = NULL;
  last_weight   = 1;
}


GHMM_Sequences::GHMM_Sequences(sequence_t* seq) {
  sequence_type = GHMM_INT;

  c_i_sequences = seq;
  c_d_sequences = NULL;
  last_weight   = 1;
}


GHMM_Sequences::GHMM_Sequences(sequence_d_t* seq) {
  sequence_type = GHMM_DOUBLE;

  c_i_sequences = NULL;
  c_d_sequences = seq;
  last_weight   = 1;
}


GHMM_Sequences::~GHMM_Sequences() {
  clean();
}


const char* GHMM_Sequences::toString() const {
  return "GHMM_Sequences";
}


/** Truncate double sequences in a given sequence array. 
    Useful for Testing;
   @return truncated sqd_field; 
   @param sqd\_in sequence arrays for truncation
   @param sqd\_arrays number of sequence arrays
   @param  trunc\_ratio 0 means  no truncation, 1 max. truncation
   @param seed rng seed
*/

//sequence_d_t **GHMM_Sequences::sequence_d_truncate(sequence_d_t **sqd_in, int sqd_arrays, 
//						  double trunc_ratio, int seed) {
//}


/**
   Determine best model for a given integer sequence. 
   Choose from the set of models the 
   one with the highest likelihood for the given sequence.
   @param mo            array of models
   @param model\_number  number of models
   @param sequence      sequence
   @param seq\_len      length of sequence
   @param log\_p         log likelihood of the sequence given the best model
   @return index of best\_model (between 0 and model\_number - 1)
*/
//int GHMM_Sequences::sequence_best_model(model **mo, int model_number, int *sequence, 
//				       int seq_len, double *log_p) {
//}


//int GHMM_Sequences::check(int max_symb) {
//  if (sequence_type == GHMM_INT)
//    return sequence_check(c_t_sequence,max_symb);

//  only_for_ints;
//}


int GHMM_Sequences::add(GHMM_Sequences* source) {
  /* sequences should contain same type of data. */
  if (source->sequence_type != sequence_type)
    return -1;

  if (sequence_type == GHMM_INT)
    return sequence_add(c_i_sequences,source->c_i_sequences);

  return sequence_d_add(c_d_sequences,source->c_d_sequences);
}


int GHMM_Sequences::add(GHMM_Sequence* source) {
  /* sequences should contain same type of data. */
  if (source->sequence_type != sequence_type)
    return -1;

  if (sequence_type == GHMM_INT)
    return sequence_add(c_i_sequences,source->c_i_sequences);

  return sequence_d_add(c_d_sequences,source->c_d_sequences);
}


void GHMM_Sequences::print(FILE *file, int discrete) const {
  if (sequence_type == GHMM_INT)
    sequence_print(file,c_i_sequences);
  if (sequence_type == GHMM_DOUBLE)
    sequence_d_print(file,c_d_sequences,discrete);
}

/**
  Prints one array of integer sequences in a xml file
  @param file       output file
  @param sequence   array of sequences
  */
//void GHMM_Sequences::sequence_print_xml(FILE *file, sequence_t *sequence) {
//}

/**
   Prints one array of integer sequences in Mathematica format.
   (List of lists)
   @param file       output file
   @param sq    array of sequences
   @param name arbitrary sequence name for usage in Mathematica.
 */
//void GHMM_Sequences::sequence_mathematica_print(FILE *file, sequence_t *sq, char *name) {
//}


/**
   Prints one array of double sequences in Mathematica format.
   (List of lists)
   @param file       output file
   @param sqd    array of sequences
   @param name arbitrary sequence name for usage in Mathematica.
 */
//void GHMM_Sequences::sequence_d_mathematica_print(FILE *file, sequence_d_t *sqd, char *name) {
//}

/** Output of double sequences suitable for gnuplot. One symbol per line,
    sequences seperated by double newline.
    @param file output file
    @param sqd array of double sequences
*/
//void GHMM_Sequences::sequence_d_gnu_print(FILE *file, sequence_d_t *sqd) {
//}


/**
   Return biggest symbol in an interger sequence.
   @param sq sequence structure
   @return max value
 */
//int GHMM_Sequences::sequence_max_symbol(sequence_t *sq) {
//}

/**
   Memory allocation for an integer sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
//sequence_t *GHMM_Sequences::sequence_calloc(long seq_number) {
//}

/**
   Memory allocation for a double  sequence struct. Allocates arrays of lenght
   seq\_number. NO allocation for the actual sequence, since its length is 
   unknown.
   @param seq\_number:  number of sequences
   @return:     pointer of sequence struct
*/
//sequence_d_t *GHMM_Sequences::sequence_d_calloc(long seq_number) {
//}

/**
   Copies array of integer sequences to double sequences.
   @return       double sequence struct (target)
   @param sq    integer sequence struct (source)
   */
//sequence_d_t *GHMM_Sequences::sequence_d_create_from_sq(const sequence_t *sq) {
//}

/**
   Copies array of double sequences into an array of integer
   sequences. Truncates positions after decimal point.
   @return       integer sequence struct (target)
   @param sq    double sequence struct (source)
   */
//sequence_t *GHMM_Sequences::sequence_create_from_sqd(const sequence_d_t *sqd) {
//}


int GHMM_Sequences::max_len() const {
  if (sequence_type == GHMM_INT)
    return sequence_max_len(c_i_sequences);

  return sequence_d_max_len(c_d_sequences);
}

/**
  Calculates a mean sequence of a given array of double sequences.
  Missing values of shorter sequences a assumed to be zero.
  @param sqd sequence struct
  @return pointer of sequence struct containing the mean sequence
  */
//sequence_d_t *GHMM_Sequences::sequence_d_mean(const sequence_d_t *sqd) {
//}

/**
   Calculates the scatter matrix of an array of double sequences. 
   Missing parts of short sequences are NOT taken into account.
   @return        scatter matrix
   @param sqd     sequence struct
   @param sqd     (calculated) dimension of scatter matrix
  */
//double **GHMM_Sequences::sequence_d_scatter_matrix(const sequence_d_t *sqd, int *dim) {
//}

/**
   Calculates transition class for a given double sequence
   at a specified position. Very application specific!!! Currently 
   implemented only dummy function: allways returns 0 which
   means no usage of multiple transition classes.
   @param O double sequence
   @param index position for class calculation
   @param osum sum of symbols upto index
   @return currently always 0
 */
//int GHMM_Sequences::sequence_d_class(const double *O, int index, double *osum) {
//}

/** Divides randomly a given array of double sequences into two sets. 
    Useful if a training and test set is needed. Memory allocation is done 
    here.
    @param sqd input sequence array
    @param sqd\_train training sequences
    @param sqd\_test test sequences
    @param train\_ratio ratio of number of train vs number of test sequences
    @return 0 for success, -1 for error
*/
//int GHMM_Sequences::sequence_d_partition(sequence_d_t *sqd, sequence_d_t * sqd_train, 
//					 sequence_d_t *sqd_test, double train_ratio) {
//}


/** 
    Copies all entries from one sequence in a source array to a target array.
    No memory allocation here.
    @param target double sequence target
    @param source double sequence source
    @param t_num position in target array
    @param s_num position in source array
*/
//void GHMM_Sequences::sequence_d_copy_all(sequence_d_t *target, long t_num, 
//					sequence_d_t *source, long s_num) {
//}

/** Log-Likelihood function in a mixture model:
    (mathe mode?)
    $\sum_k w^k \log( \sum_c (\alpha_c p(O^k | \lambda_c)))$
    @param smo pointer to array of smodels
    @param smo\_number number of models
    @param sqd sequence struct
    @param like log likelihood
*/
//int GHMM_Sequences::sequence_d_mix_like(smodel **smo, int  smo_number, sequence_d_t *sqd, double *like) {
//}


void GHMM_Sequences::clean() {
  if (c_i_sequences)
    sequence_free(&c_i_sequences);
  if (c_d_sequences)
    sequence_d_free(&c_d_sequences);

  c_i_sequences = NULL;
  c_d_sequences = NULL;

  clean_cpp();
}


int* GHMM_Sequences::getIntSequence(int index) const {
  if (!c_i_sequences) {
    fprintf(stderr,"GHMM_Sequences::getIntSequence(int) object does not contain int sequences.\n");
    exit(1);
  }

  return c_i_sequences->seq[index];
}


double* GHMM_Sequences::getDoubleSequence(int index) const {
  if (!c_d_sequences) {
    fprintf(stderr,"GHMM_Sequences::getDoubleSequence(int) object does not contain double sequences.\n");
    exit(1);
  }

  return c_d_sequences->seq[index];
}


unsigned int GHMM_Sequences::getLength(int index) const {
  if (c_d_sequences)
    return c_d_sequences->seq_len[index];

  if (c_i_sequences)
    return c_i_sequences->seq_len[index];

  fprintf(stderr,"GHMM_Sequences::getLength(int) object does not contain int sequences.\n");
  exit(1);

  return 0;
}


XMLIO_Element* GHMM_Sequences::XMLIO_startTag(const string& my_tag, XMLIO_Attributes& my_attributes) {
  if (my_tag == "sequence") {
    GHMM_Sequence* sequence = new GHMM_Sequence(sequence_type,last_weight);
    sequences.push_back(sequence);

    last_weight = 1;

    return sequence;
  }

  fprintf(stderr,"unexpected tag: %s in <sequences>\n",my_tag.c_str());
  exit(1);

  return NULL;
}


void GHMM_Sequences::XMLIO_finishedReading() {
  if (sequence_type == GHMM_INT)
    c_i_sequences = sequence_calloc(0);

  if (sequence_type == GHMM_DOUBLE)
    c_d_sequences = sequence_d_calloc(0);

  unsigned int i;
  for (i = 0; i < sequences.size(); ++i)
    add(sequences[i]);

  clean_cpp();
}


void GHMM_Sequences::clean_cpp() {
  unsigned int i;
  for (i = 0; i < sequences.size(); ++i)
    SAFE_DELETE(sequences[i]);
  sequences.clear();
}


void GHMM_Sequences::read(const string& filename) {
  int seq_number;
  int i;

  /* first clean up old data. */
  clean();

  if (sequence_type == GHMM_DOUBLE) {
    sequence_d_t** seq = sequence_d_read(filename.c_str(),&seq_number);

    if (seq_number > 0)
      copyFromSequences(seq[0]);

    for (i = 0; i < seq_number; ++i)
      sequence_d_free(&seq[i]);

    free(seq);
  }
  else {
    sequence_t** seq = sequence_read(filename.c_str(),&seq_number);

    if (seq_number > 0)
      copyFromSequences(seq[0]);

    for (i = 0; i < seq_number; ++i)
      sequence_free(&seq[i]);

    free(seq);
  }
}


void GHMM_Sequences::copyFromSequences(sequence_d_t* seq) {
  int i;
  int j;

  clean();

  c_d_sequences = sequence_d_calloc(seq->seq_number);

  for (i = 0; i < seq->seq_number; ++i) {
    c_d_sequences->seq[i] = (double*) malloc(sizeof(double) * seq->seq_len[i]);
    for (j = 0; j < seq->seq_len[i]; ++j)
      c_d_sequences->seq[i][j] = seq->seq[i][j];

    c_d_sequences->seq_len[i]   = seq->seq_len[i];
    c_d_sequences->seq_label[i] = seq->seq_label[i];
    c_d_sequences->seq_id[i]    = seq->seq_id[i];
    c_d_sequences->seq_w[i]     = seq->seq_w[i];
  }
  c_d_sequences->seq_number = seq->seq_number;
  c_d_sequences->total_w    = seq->total_w;
}


void GHMM_Sequences::copyFromSequences(sequence_t* seq) {
  int i;
  int j;

  clean();

  c_i_sequences = sequence_calloc(seq->seq_number);

  for (i = 0; i < seq->seq_number; ++i) {
    c_i_sequences->seq[i] = (int*) malloc(sizeof(double) * seq->seq_len[i]);
    for (j = 0; j < seq->seq_len[i]; ++j)
      c_d_sequences->seq[i][j] = seq->seq[i][j];

    c_i_sequences->seq_len[i]   = seq->seq_len[i];
    c_i_sequences->seq_label[i] = seq->seq_label[i];
    c_i_sequences->seq_id[i]    = seq->seq_id[i];
    c_i_sequences->seq_w[i]     = seq->seq_w[i];
  }
  c_i_sequences->seq_number = seq->seq_number;
  c_i_sequences->total_w    = seq->total_w;
}
