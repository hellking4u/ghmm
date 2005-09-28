/*
  author       : David Posada <dposada@variagenics.com>
  filename     : ghmm/tests/sequences_old_format.c
  created      : DATE: Februar 2002
  $Id$
  __copyright__
*/

/* should be corrected with ghmm/sequence.c version 1.9 */
#include <stdio.h>
#include <ghmm/sequence.h>

#include <ghmm/obsolete.h>
   
int main()
{
#ifdef GHMM_OBSOLETE
    int test_result=0;
    const char* double_sequences_file="data/test100.sqd";
    sequence_d_t **sqd = NULL;
    int sqd_number;
    const char* int_sequences_file="data/sequences_old_format.sq";    
    sequence_t **data = NULL;
    int data_number;

    /* read double sequences (this works fine)*/
    fprintf(stderr,"reading double sequences from %s ...",double_sequences_file);
    sqd=sequence_d_read((char*)double_sequences_file, &sqd_number);
    if (sqd==NULL) {
      test_result=1;
      fprintf(stdout, " Failed\n");
    }
    else {
      fprintf(stdout," Done\n");
      sequence_d_free(sqd);
    }


     /* read int sequences (this gives a segmentation fault)*/
    fprintf(stderr,"reading int sequences from %s ...",int_sequences_file);
    data=sequence_read((char*)int_sequences_file,&data_number);
    if (data==NULL) {
      test_result=1;
      fprintf(stdout, " Failed\n");
    }
    else {
      fprintf(stdout," Done\n");
      sequence_free(data);
    }
    return test_result;
#else /* GHMM_OBSOLETE */
    return 0;
#endif /* GHMM_OBSOLETE */
}
