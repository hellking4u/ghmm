#include "ghmm++/GHMM.h"

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

  if (argc<2) {
    fprintf(stderr,"generiere_nvs needs a file as argument\n");
    exit(1);
  }

  if (0!=doc.open(argv[1],"r")) {
    fprintf(stderr,"generiere_nvs could not open file %s\n",argv[1]);
    exit(1);
  }

  doc.readDocument();
  doc.close();

  smo          = doc.getContinuousModel();
  trainingsseq = doc.getSequences();

  smo->reestimate_baum_welch(trainingsseq,&logp,0.00001,35);

  smo->print(stdout);
  for (j = 0; j < 200; j++) {
    genseq = smo->generate_sequences(0,20,1,1,20);
    for (i = 0; i < genseq->getLength(0); i++) 
      printf("%8f\n",(genseq->getDoubleSequence(0)[i]));
  }

  if (genseq)
    delete genseq;

  delete smo;
  delete trainingsseq;

  //printf("Generierte Sequenzen :\n");
  // sequence_d_print(stdout,genseq,0);
  //printf("Eingabesequenzen: \n");
  //sequence_d_print(stdout,seq[0],0);

  //smodel_print(stdout,smo[0]);

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

 return 0;  
}
