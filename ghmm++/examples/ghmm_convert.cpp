/*******************************************************************************
  author       : Bernd Wichern
  filename     : /zpr/bspk/src/hmm/ghmm/ghmm++/ghmm_convert.cpp
  created      : TIME: 10:13:13     DATE: Thu 13. September 2001
  last-modified: TIME: 11:33:44     DATE: Fri 14. September 2001
  $Id$
*******************************************************************************/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <xmlio/XMLIO_Document.h>
#include <ghmm++/ghmm.h>
#include <ghmm++/sequences_document.h>
#include <ghmm/sequence.h>

#ifdef HAVE_NAMESPACES
using namespace std;
#endif



int convert() {
  int i, sqd_number;
  sequence_d_t ** sqd = NULL;
  string seq_file_name;
  cerr<<"convert sequences from file: "<<endl;
  cin>>seq_file_name;
  cerr<<"Seq.File "<<seq_file_name<<endl;

  sqd = sequence_d_read((char *)seq_file_name.c_str(), &sqd_number);
  if (!sqd) {
    cerr<<"unable to read seqs. from file "<<seq_file_name<<endl;
    return -1;
  }
  
  cout<<"sqdnumber "<<sqd_number<<endl;

  // sequence_d_print(stdout, sqd[i], 0);
  /* make ghmm++ objects and print in xml format */
  sequences_document output(sqd, sqd_number);
  output.write_sequences("/dev/stdout");

  return 0;
}


int main()
{

  return convert();
}




