/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 * 
 * __copyright__
 */

#include "ghmm++/GHMM_Sequence.h"
#include "ghmm++/GHMM_Alphabet.h"


#ifdef HAVE_NAMESPACES
using namespace std;
#endif


GHMM_Sequence::GHMM_Sequence(GHMM_SequenceType my_sequence_type, int len, int my_weight) {
  switch (my_sequence_type) {

  case GHMM_INT:
    init_INT(NULL,len,my_weight);
    break;

  case GHMM_DOUBLE:
    init_DOUBLE(len,my_weight);
    break;
  }
}


GHMM_Sequence::GHMM_Sequence(GHMM_Alphabet* my_alphabet, int len, int my_weight) {
  init_INT(my_alphabet,len,my_weight);
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
  unsigned int i;

  XMLIO_ArrayElement<string>::XMLIO_finishedReading();

  if (sequence_type == GHMM_INT) {
    string seq;
    
    c_i_sequences->seq_w[0] = weight;
    
    for (i = 0; i < size(); ++i)
      seq += (*this)[i];
    
    resize(seq.length());
    
    for (i = 0; i < seq.length(); ++i)
      c_i_sequences->seq[0][i] = alphabet->getIndex(seq.substr(i,1));
  }

  if (sequence_type == GHMM_DOUBLE) {
    resize(size());

    c_d_sequences->seq_w[0] = weight;

    for (i = 0; i < size(); ++i)
      c_d_sequences->seq[0][i] = atof((*this)[i].c_str());
  }

  /* clear array. */
  clear();
}


void GHMM_Sequence::resize(int new_len) {
  int new_alloc_len = MAX(1,new_len);

  switch (sequence_type) {

  case GHMM_INT:
    c_i_sequences->seq[0]     = (int*) realloc(c_i_sequences->seq[0],sizeof(int) * new_alloc_len);
    c_i_sequences->seq_len[0] = new_len;
    break;

  case GHMM_DOUBLE:
    c_d_sequences->seq[0]     = (double*) realloc(c_d_sequences->seq[0],sizeof(double) * new_alloc_len);
    c_d_sequences->seq_len[0] = new_len;
    break;
  }
}


void GHMM_Sequence::init() {
  c_i_sequences = NULL;
  c_d_sequences = NULL;
  alphabet      = NULL;
}


void GHMM_Sequence::init_INT(GHMM_Alphabet* my_alphabet, int len, int my_weight) {
  init();

  sequence_type = GHMM_INT;
  c_i_sequences = sequence_calloc(1);
  alphabet      = my_alphabet;

  resize(len);
  weight = my_weight;
}


void GHMM_Sequence::init_DOUBLE(int len, int my_weight) {
  init();

  sequence_type = GHMM_DOUBLE;
  c_d_sequences = sequence_d_calloc(1);

  resize(len);
  weight = my_weight;
}


void GHMM_Sequence::setDouble(int index, double value) {
  if (c_d_sequences)
    if (c_d_sequences->seq_len[0] > index)
      c_d_sequences->seq[0][index] = value;
}


void GHMM_Sequence::setInt(int index, int value) {
  if (c_i_sequences)
    if (c_i_sequences->seq_len[0] > index)
      c_i_sequences->seq[0][index] = value;
}


void GHMM_Sequence::print(FILE *file, int discrete) const {
  if (sequence_type == GHMM_INT)
    sequence_print(file,c_i_sequences);
  if (sequence_type == GHMM_DOUBLE)
    sequence_d_print(file,c_d_sequences,discrete);
}
