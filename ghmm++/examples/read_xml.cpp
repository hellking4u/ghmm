/*
 *
 * created: 26 Feb 2002 by Wasinee Rungsarityotin
 * authors: Wasinee Rungsarityotin (rungsari@molgen.mpg.de)
 * file   : $Source$
 * $Id$
 * revision date   : $Date$
  ___copyright__
 
 */

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <ghmm/rng.h>
#include "ghmm++/GHMM.h"
#include "ghmm++/GHMM_GMLDoc.h"
#include <iostream>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int main(int argc, char* argv[]) {
  GHMM_Document      doc;
  // GHMM_GraphMLDoc    doc;
  GHMM_DiscreteModel *dmo;
  GHMM_Sequences* training_seq = NULL;
  GHMM_Sequences* my_output = NULL;
  GHMM_Alphabet*  alphas = NULL;

  /* Important! initialise rng  */
  GHMM_Toolkit::gsl_rng_init();


  doc.open("dna.xml", "r");
  doc.readDocument();
  doc.close();

  dmo  = doc.getDiscreteModel();

  //dmo->setTransition(0,0,0.0);
  //dmo->setTransition(0,1,1.0);

  if (dmo->check() == -1)
    return 1;

  fprintf(stdout, "initial distribution:\n");
  dmo->Pi_print(stdout,"  ",",","\n");
  fprintf(stdout,"transition matrix:\n");
  dmo->A_print(stdout,""," ","\n");
  fprintf(stdout,"observation symbol matrix:\n");
  dmo->B_print(stdout,""," ","\n");

  alphas = dmo->getAlphabet();
  assert( alphas != NULL );
  for(int i=0; i < alphas->size(); i++)
    {
      cout << "Symbol " << i << ":";
      cout << alphas->getSymbol(i) << endl;
    }

  if ( ( training_seq = doc.getSequences() ) != NULL )
    {
      training_seq->print(stdout);
      
      printf("reestimating with Baum-Welch-algorithm...");
      dmo->reestimate_baum_welch( training_seq );

      /* print the result */
      printf("Done\nthe result is:\n");  
      dmo->A_print(stdout,""," ","\n");
      fprintf(stdout,"observation symbol matrix:\n");
      dmo->B_print(stdout,""," ","\n");

      my_output = dmo->generate_sequences(1,10,10);
      my_output->print(stdout);

      doc.open("dna2.xml","w"); // Try writing to a new format!!!!
      doc.writeElement(dmo);
      doc.writeEndl();
      doc.writeEndl();
      //doc.writeElement(training_seq);
      doc.writeEndl();
      doc.close();
      
      cout << "Wrote to dna2.xml" << endl;
    }

  delete my_output;
  delete alphas;

  SAFE_DELETE(dmo);

  return 0;
}
