#include <ghmm/pmodel.h>
#include <ghmm/mysequence.h>
#include <ghmm/mes.h>

#ifndef PCLASSCHANGE_H
#define PCLASSCHANGE_H

struct threshold_user_data {
  /** which double value in myseq **/
  int seq_index;
  /** cut off value **/
  double threshold;
  /** add this to the index in sequence X **/
  int offset_x;
  /** add this to the index in sequence Y **/
  int offset_y;
};
typedef struct threshold_user_data threshold_user_data;

struct boolean_user_data {
  /** which double value in myseq **/
  int seq_index;
  /** add this to the index in sequence X **/
  int offset_x;
  /** add this to the index in sequence Y **/
  int offset_y;
};
typedef struct boolean_user_data boolean_user_data;

int gt_sum(pmodel * mo, mysequence * X, mysequence * Y, int index_x, int index_y, void * user_data);

int lt_sum(pmodel * mo, mysequence * X, mysequence * Y, int index_x, int index_y, void * user_data);

int boolean_and(pmodel * mo, mysequence * X, mysequence * Y, int index_x, int index_y, void * user_data);

int boolean_or(pmodel * mo, mysequence * X, mysequence * Y, int index_x, int index_y, void * user_data);

void set_to_lt_sum(pclass_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y);

void set_to_gt_sum(pclass_change_context * pccc, int seq_index, double threshold, int offset_x, int offset_y);

void set_to_boolean_and(pclass_change_context * pccc, int seq_index, int offset_x, int offset_y);

void set_to_boolean_or(pclass_change_context * pccc, int seq_index, int offset_x, int offset_y);
#endif
