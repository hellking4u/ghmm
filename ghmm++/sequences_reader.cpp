#include "sequences.h"

int main() {
  sequence_t* seq;
  char* file="ghmm.xml";
  seq=return_sequences(file);
  sequence_print_xml(stdout,seq);
  sequence_free(&seq);
  return 0;
}
