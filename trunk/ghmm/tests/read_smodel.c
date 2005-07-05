#include <stdio.h>
#include <ghmm/smodel.h>

int main()
{
  smodel **model_array;
  int model_counter;
  char *inFileName = "sample.smo";

  ghmm_rng_init();  /* Important! initialise rng  */
  
  model_array = smodel_read(inFileName,&model_counter);
  printf("Read %d model(s) from '%s'\n", model_counter, inFileName);
  smodel_free(model_array);
}
  
