#ifndef KBEST_H
#define KBEST_H

#ifdef __cplusplus
extern "C" {
#endif

/**
   Calculates the most probable labeling for the given sequence in the given
   model using k-best decoding.
   Labels must be from interval [0:max_label] without gaps!!! (not checked)
   Model must not have silent states. (checked in Python wrapper)
   @return array of labels (internal representation)
   @param mo:         pointer to a model
   @param o_seq:      output sequence (array of internal representation chars)
   @param seq_len:    length of output sequence
   @param k:          number of hypotheses to keep for each state
   @param log_p:      variable reference to store the log prob. of the labeling
 */
int* kbest(model* mo, int* o_seq, int seq_len, int k, double* log_p);

#ifdef __cplusplus
}
#endif

#endif
