/*******************************************************************************
  authors      : Peter Pipenbacher
  filename     : ghmm++/examples/dna.cpp
  created      : DATE: 2002-03-14
  $Id$

  __copyright__
*******************************************************************************/

#include "ghmm++/GHMM.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif


int dna_test() {
  GHMM_Alphabet alphabet;
  alphabet.addSymbol("A");
  alphabet.addSymbol("C");
  alphabet.addSymbol("T");
  alphabet.addSymbol("G");
  alphabet.addSymbol("-");

  GHMM_DiscreteModel model(&alphabet);

  model.addState("StartRegion");
  model.addState("EndRegion");

  model.getState("StartRegion")->setInitialProbability(0.6);
  model.getState("StartRegion")->setOutputProbability("A",0.2);
  model.getState("StartRegion")->setOutputProbability("C",0.2);
  model.getState("StartRegion")->setOutputProbability("T",0.2);
  model.getState("StartRegion")->setOutputProbability("G",0.2);
  model.getState("StartRegion")->setOutputProbability("-",0.2);

  model.getState("EndRegion")->setInitialProbability(0.4);
  model.getState("EndRegion")->setOutputProbability("A",0.0);
  model.getState("EndRegion")->setOutputProbability("C",0.5);
  model.getState("EndRegion")->setOutputProbability("T",0.0);
  model.getState("EndRegion")->setOutputProbability("G",0.5);
  
  model.setTransition("EndRegion","EndRegion",0.9);
  model.setTransition("EndRegion","StartRegion",0.1);
  model.setTransition("StartRegion","EndRegion",0.1);
  model.setTransition("StartRegion","StartRegion",0.9);

  if (model.check() != 0)
    return 1;

  GHMM_Sequences seq(&alphabet);
  seq.add("AAACGCGCGCGCGCG");
  seq.add("CGCTTGGCGCGCGCGCG");
  seq.add("TAGCA-CGCGCGCGCG");
  seq.add("-CC-AGCGCGCGCGCG");
  seq.add("--CTGATTCGCGCGCGCG");
  seq.add("C-G-GGCGCGCGCGCG");
  seq.add("A--CGCGCGCGCGCG");
  seq.add("------CGCGCGCGCG");
  seq.add("CGATAACTTGCGCGCGCG");
  seq.add("CGGG--TCGCGCGCGCG");

  model.reestimate_baum_welch(&seq);
  model.A_print(stdout,""," ","\n"); 
  model.B_print(stdout,""," ","\n"); 

  GHMM_Document doc;
  doc.open("dna.xml","w");
  doc.writeElement(&model);
  doc.writeEndl();
  doc.writeEndl();
  doc.writeElement(&seq);
  doc.writeEndl();
  doc.close();

  doc.open("dna.xml","r");
  doc.readDocument();
  doc.close();

  GHMM_DiscreteModel* new_model = doc.getDiscreteModel();
  GHMM_Sequences* new_seq = doc.getSequences();

  doc.open("dna2.xml","w");
  doc.writeElement(new_model);
  doc.writeEndl();
  doc.writeEndl();
  doc.writeElement(new_seq);
  doc.writeEndl();
  doc.close();
  
  cout << "Wrote to dna2.xml" << endl;

  delete new_model;
  delete new_seq;

  //SAFE_DELETE(new_model);
  //SAFE_DELETE(new_seq);

  return 0;
}


int main() {
  /* Important! initialise rng  */
  GHMM_Toolkit::ghmm_rng_init();

  int result = dna_test();

#ifdef WIN32
  printf("\nPress ENTER\n");
  fgetc(stdin);
#endif

  return result;
}
