#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/sdmodel.h>

int cp_class_change(int *seq, int len) {
  int sum = 0;
  int i;
  for(i=0;i<=len;i++){
	  sum += seq[i];
  }
  //printf("sum = %d\n",sum);
  if (sum >= 6) {
    //printf("\n++++++++++++++++++++++++++++++++Switching class .... ");    
    return 1;
  }
  else {
    return 0;
  } 
} 		
		

void setSwitchingFunction( sdmodel *smd ) {
  smd->get_class = cp_class_change;
}

