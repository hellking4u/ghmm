#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main(int argc, char* argv[]) {
  unsigned int i;
  int j;
  double logp;
  
  GHMM_Sequences trainingsseq(GHMM_DOUBLE);
  GHMM_Sequences* genseq;
  GHMM_ContinuousModel smo;

  if (argc < 3) {
    fprintf(stderr,"usage: generiere_nvs_old_files nvs.mod nvs.sqd\n");
    exit(1);
  }

  smo.read(argv[1]);
  trainingsseq.read(argv[2]);

  smo.reestimate_baum_welch(&trainingsseq,&logp,0.00001,35);

  smo.print(stdout);
  for (j = 0; j < 200; j++) {
    genseq = smo.generate_sequences(0,20,1,1,20);
    for (i = 0; i < genseq->getLength(0); i++) 
      printf("%8f\n",(genseq->getDoubleSequence(0)[i]));
  }

  delete genseq;

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
