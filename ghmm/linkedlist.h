#ifndef LINKEDLIST_H
#define LINKEDLIST_H

/* integer element */
struct i_el {
  int val;
  struct i_el * next;
};
typedef struct i_el i_el;

/* list header */
struct i_list {
  i_el * first;
  i_el * last;
  int length;
};
typedef struct i_list i_list;

i_list * ighmm_list_init_list();

int ighmm_list_free(i_list * list);

void ighmm_list_append(i_list * list, int val);

void ighmm_list_insert(i_list * list, int val);

int * ighmm_list_to_array(i_list * list);

i_el * ighmm_list_init_el(int val);
#endif
