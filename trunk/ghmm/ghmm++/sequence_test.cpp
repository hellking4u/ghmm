#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>
#include <iostream>
#include "sequences.h"

#ifdef HAVE_NAMESPACES
using namespace std;
#endif

int sequenceReader_test()
{
  string filename="ghmm.xml";
  sequenceReader seq;
  (void)seq.read_sequences(filename);
  cout<<"read "<<seq.size()<<" sequences"<<endl;

  sequenceReader::const_iterator pos=seq.begin();
  while(pos!=seq.end())
    {
      const sequences* these_sequences=*pos;
      if (these_sequences!=NULL)
	{
	  if (these_sequences->get_type()=="int")
	    {
	      sequence_t* my_seq=these_sequences->create_sequence_t();
	      sequence_print(stdout,my_seq);
	      sequence_free(&my_seq);
	    }
	  else if (these_sequences->get_type()=="double")
	    {
	      sequence_d_t* my_seq=these_sequences->create_sequence_d_t();
	      sequence_d_print(stdout,my_seq,0);
	      sequence_d_free(&my_seq); 
	    }
	}
      ++pos;
    }
  return 0;
}

int main()
{
  return sequenceReader_test();

}





