#include "linkedlist.h"
#include <stdlib.h>
#include "mes.h"

void i_list_append(i_list * list, int val){
  i_el * last;
  i_el * el;
  el = init_i_el(val);
  if (list->first == NULL) {
    list->first = el;
    list->last = el;
  }
  else {
    last = list->last;
    last->next = el;
    list->last = el;
  }
  list->length++;
}

void i_list_insert(i_list * list, int val) {
  i_el * first;
  i_el * el;
  el = init_i_el(val);
  if (list->first == NULL) {
    list->first = el;
    list->last = el;
  }
  else {
    first = list->first;
    el->next = first;
    list->first = el;
  }
  list->length++;
}

void i_list_print(i_list * list) {
  i_el * el = list->first;
  printf("LIST : ");
  while(el != NULL) {
    printf("%i, ", el->val);
    el = el->next;
  }
  printf("\n");
}

int * i_list_to_array(i_list * list) {
#define CUR_PROC "i_list_to_array"
  int * array;
  if (!m_calloc(array, list->length)) {mes_proc(); goto STOP;}
  i_el * el;
  int counter = 0;
  el = list->first;
  while(el != NULL) {
    array[counter] = el->val;
    el = el->next;
    counter++;
  }
  return array;
STOP:
  free(array);
  return NULL;
#undef CUR_PROC
}

i_list * init_i_list() {
#define CUR_PROC "init_i_list"
  i_list * list;
  list = m_calloc(list, 1);
  if (!list) {mes_proc(); goto STOP;}
  list->first = NULL;
  list->last = NULL;
  list->length = 0;
  return list;
STOP:
  free_i_list(list);
  return NULL;
#undef CUR_PROC
}

i_el * init_i_el(int val) {
#define CUR_PROC "init_i_el"
  i_el * el;
  el = m_calloc(el, 1);
  if (!el) {mes_proc(); goto STOP;}
  el->next = NULL;
  el->val = val;
  return el;
STOP:
  free(el);
  return NULL;
#undef CUR_PROC
}

int free_i_list(i_list * list) {
  i_el * el;
  i_el * next;
  el = list->first;
  while(el != NULL) {
    next = el->next;
    free(el);
    el = next;
  }
}

