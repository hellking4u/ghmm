#include "ppghmm++/GHMM_Sequences.h"
#include "ppghmm++/GHMM_ContinuousModel.h"
#include "ppghmm++/GHMM_Document.h"
#ifdef HAVE_CMATH
#  include <cmath>
#else
#  include <math.h>
#endif

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main(int argc, char* argv[]) {
  unsigned int i;
  int j;
  double logp;
  
  GHMM_Sequences* trainingsseq;
  GHMM_Sequences* genseq;
  GHMM_ContinuousModel* smo;
  GHMM_Document doc;

  doc.open(argv[1],"r");
  doc.readDocument();
  doc.close();

  smo          = doc.getContinuousModel();
  trainingsseq = doc.getSequences();

  smo->reestimate_baum_welch(trainingsseq,&logp,0.00001,1000);

  smo->print(stdout);
  for (j = 0; j < 200; j++) {
    genseq = smo->generate_sequences(0,20,1,1,20);
    for (i = 0; i < genseq->getLength(0); i++) 
      printf("%8f\n",(genseq->getDoubleSequence(0)[i]));
  }

  delete genseq;

  //printf("Generierte Sequenzen :\n");
  // sequence_d_print(stdout,genseq,0);
  //printf("Eingabesequenzen: \n");
  //sequence_d_print(stdout,seq[0],0);

  //smodel_print(stdout,smo[0]);

 return 1;  
}
