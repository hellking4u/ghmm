#include "pclasschange.h"

/* classification functions for the pairhmm genefinder 
   transition class 0 is the default class */

/* if the sum of the ka values is less than the threshold return 1 */
int lt_sum(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) {
  /* cast the user data */
  threshold_user_data * td = (threshold_user_data *)user_data;
  if (get_double_psequence(X, td->seq_index, index_x + td->offset_x) + 
      get_double_psequence(Y, td->seq_index, index_y + td->offset_y) < 
      td->threshold)
    return 1;
  else
    return 0;
}

/* if the sum of the ka values is greater than the threshold return 1 */
int gt_sum(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) {
  /* cast the user data */
  threshold_user_data * td = (threshold_user_data *)user_data;
  if (get_double_psequence(X, td->seq_index, index_x + td->offset_x) + 
      get_double_psequence(Y, td->seq_index, index_y + td->offset_y) >
      td->threshold)
    return 1;
  else
    return 0;
}

/* reads int sequences of boolean precomputed classification. If both positions
   are true return true else false */
int boolean_and(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) {
  boolean_user_data * td = (boolean_user_data *)user_data;
  if (get_char_psequence(X, td->seq_index, index_x + td->offset_x) && 
      get_char_psequence(Y, td->seq_index, index_y + td->offset_y))
    return 1;
  else
    return 0;
}

/* reads int sequences of boolean precomputed classification. If one position
   is true return true else false */
int boolean_or(pmodel * mo, psequence * X, psequence * Y, int index_x, int index_y, void * user_data) {
  boolean_user_data * td = (boolean_user_data *)user_data;
  if (get_char_psequence(X, td->seq_index, index_x + td->offset_x) ||
      get_char_psequence(Y, td->seq_index, index_y + td->offset_y))
    return 1;
  else
    return 0;
}

void set_to_lt_sum(pclass_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y) {
  if (pccc) {
    threshold_user_data * td;
    m_calloc(td, 1);
    td->seq_index = seq_index;
    td->threshold = threshold;
    td->offset_x = offset_x;
    td->offset_y = offset_y;
    pccc->user_data = td;
    pccc->get_class = &lt_sum;
  }
  else
    fprintf(stderr, "set_to_lt_sum_ka: No class change context\n");
}

void set_to_gt_sum(pclass_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y) {
  if (pccc){
    threshold_user_data * td;
    m_calloc(td, 1);
    td->seq_index = seq_index;
    td->threshold = threshold;
    td->offset_x = offset_x;
    td->offset_y = offset_y;
    pccc->user_data = td;
    pccc->get_class = &gt_sum;
  }
  else
    fprintf(stderr, "set_to_gt_sum_deltaka: No class change context\n");
}

void set_to_boolean_and(pclass_change_context * pccc, int seq_index, int offset_x, int offset_y) {
  if (pccc){
    boolean_user_data * td;
    m_calloc(td, 1);
    td->seq_index = seq_index;
    td->offset_x = offset_x;
    td->offset_y = offset_y;
    pccc->user_data = td;
    pccc->get_class = &boolean_and;
  }
  else
    fprintf(stderr, "set_to_boolean_and: No class change context\n");
}

void set_to_boolean_or(pclass_change_context * pccc, int seq_index, int offset_x, int offset_y) {
  if (pccc){
    boolean_user_data * td;
    m_calloc(td, 1);
    td->seq_index = seq_index;
    td->offset_x = offset_x;
    td->offset_y = offset_y;
    pccc->user_data = td;
    pccc->get_class = &boolean_or;
  }
  else
    fprintf(stderr, "set_to_boolean_and: No class change context\n");
}

