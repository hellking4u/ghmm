#include <stdio.h>
#include <stdlib.h>
#include <ghmm/rng.h>
#include <ghmm/sequence.h>
#include <ghmm/sdmodel.h>

int cp_class_change( smodel *smo, int *seq, int k, int t);
void setSwitchingFunction( smodel *smd );

int python_class_change( smodel* smo, int* seq, int k, int t );
void setPythonSwitching( smodel *smd, char* python_module, char* python_function);
