/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#include "ppghmm++/GHMM_Sequence.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Sequence::GHMM_Sequence(GHMM_SequenceType my_sequence_type, int my_weight) {
  sequence_type = my_sequence_type;

  c_i_sequences = NULL;
  c_d_sequences = NULL;
  weight        = my_weight;
}


GHMM_Sequence::~GHMM_Sequence() {
  if (c_i_sequences)
    sequence_free(&c_i_sequences);
  if (c_d_sequences)
    sequence_d_free(&c_d_sequences);
  
  c_i_sequences = NULL;
  c_d_sequences = NULL;
}


const char* GHMM_Sequence::toString() const {
  return "GHMM_Sequence";
}


void GHMM_Sequence::XMLIO_finishedReading() {
  XMLIO_ArrayElement<string>::XMLIO_finishedReading();

  if (sequence_type == GHMM_INT) {
    c_i_sequences = sequence_calloc(1);

    c_i_sequences->seq[0]     = (int*) malloc(sizeof(int) * size());
    c_i_sequences->seq_len[0] = size();
    c_i_sequences->seq_w[0]   = weight;
  }

  if (sequence_type == GHMM_DOUBLE) {
    c_d_sequences = sequence_d_calloc(1);

    c_d_sequences->seq[0]     = (double*) malloc(sizeof(double) * size());
    c_d_sequences->seq_len[0] = size();
    c_d_sequences->seq_w[0]   = weight;
  }

  unsigned int i;
  for (i = 0; i < size(); ++i) {
    if (sequence_type == GHMM_DOUBLE)
      c_d_sequences->seq[0][i] = atof((*this)[i].c_str());
    if (sequence_type == GHMM_INT)
      c_d_sequences->seq[0][i] = atoi((*this)[i].c_str());
  }

  /* clear array. */
  clear();
}
