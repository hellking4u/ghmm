#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "sequences.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main() {
  sequence_t* seq;
  char* file="ghmm.xml";
  seq=return_sequences(file);
  sequence_print_xml(stdout,seq);
  sequence_free(&seq);
  return 0;
}
