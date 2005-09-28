/*******************************************************************************
  author       : Achim Gädke
  filename     : ghmm/tests/sequences_test.c
  created      : DATE: Thu 26. June 2001
  $Id$

__copyright__

*****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ghmm/sequence.h>

void sequence_alloc_print(void)
{
  sequence_t* seq_array;
  int i;

  seq_array= sequence_calloc(1);
  seq_array->seq_len[0]=10;
#ifdef GHMM_OBSOLETE
  seq_array->seq_label[0]=100;
#endif /* GHMM_OBSOLETE */
  seq_array->seq_id[0]=101.0;
  seq_array->seq[0]=(int*)malloc(seq_array->seq_len[0]*sizeof(int));

  for (i=0; i<seq_array->seq_len[0]; i++)
    seq_array->seq[0][i]=1;

  sequence_print_xml(stdout,seq_array);

  sequence_free(&seq_array);
}

int main()
{
  sequence_alloc_print();
  return 0;
}
