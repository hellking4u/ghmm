/*
 * created: 29 Jan 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_SEQUENCES_H
#define _GHMM_SEQUENCES_H 1

#include <vector>
#include <xmlio/XMLIO_Element.h>
#include "ghmm/sequence.h"
#include "ppghmm++/GHMM_Types.h"

#include <ppghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Sequences;
class GHMM_Sequence;

/** This can either be a sequence of double sequences or int sequences.*/
class GHMM_Sequences: public XMLIO_Element {
  
 public:
  
  /** For reading from xml file. */
  GHMM_Sequences(XMLIO_Attributes& attrs);
  /** Constructor. */
  GHMM_Sequences(GHMM_SequenceType sequence_type);
  /** Constructor. Constructs object from C structure.
      Object now owns the sequence_t. */
  GHMM_Sequences(sequence_t* seq);
  /** Constructor. Constructs object from C structure.
      Object now owns the sequence_d_t. */
  GHMM_Sequences(sequence_d_t* seq);
  /** Destructor. */
  virtual ~GHMM_Sequences();
  
  /** Returns name of class. */
  virtual const char* toString() const;
  
  /**
     Adds all sequences, sequence lengths etc 
     from source to this object. Memory allocation is done here.
     @param source  source sequence structure
     @return -1 for error, 0 for success
  */
  int add(GHMM_Sequences* source);
  /**
     Adds all sequences, sequence lengths etc 
     from source to this object. Memory allocation is done here.
     @param source  source sequence structure
     @return -1 for error, 0 for success
  */
  int add(GHMM_Sequence* source);
  /** */
  void clean_cpp();
  /** Copy content from C sequence to this object. */
  void copyFromSequences(sequence_t* seq);
  /** Copy content from C sequence to this object. */
  void copyFromSequences(sequence_d_t* seq);
  /** Returns double sequence with given index. 
      Aborts program if object doesn't contain double sequences. */
  double* getDoubleSequence(int index) const;
  /** Returns integer sequence with given index.
      Aborts program if object doesn't contain integer sequences. */
  int* getIntSequence(int index) const;
  /** Returns length of sequence with given index. */
  unsigned int getLength(int index) const;
  /**
     Make sure that the sequences only contain allowed symbols. 
     (between 0 and max\_symbol - 1)
     @param max_symb    number of different symbols
     @return            -1 for error, 0 for no errors
  */
  //  int check(int max_symb);
  /**
     Frees allocated memory of C structs and sets pointer to NULL.
  */
  void clean();
  /** 
      Determines max length of all sequences.
      @return max sequence length
  */
  int max_len() const;
  /**
     Prints one array of integer sequences in a file.
     @param file       output file
     @param discrete   switch: 0 means double output for symbols,  
                       1 means truncate symbols to integer
  */
  void print(FILE *file, int discrete = 0) const;
  /**
     Reads one or several arrays of sequences. 
     @param filename    input filename
  */
  void read(const string& filename);

  /** */
  virtual void XMLIO_finishedReading();
  /** */
  virtual XMLIO_Element* XMLIO_startTag(const string& my_tag, XMLIO_Attributes& my_attributes);

  /** Type of current sequence. */
  GHMM_SequenceType sequence_type;

  /** C type int sequence. */
  sequence_t*   c_i_sequences;
  /** C type double sequence. */
  sequence_d_t* c_d_sequences;
  /** */
  int last_weight;

 private:
  
  vector<GHMM_Sequence*> sequences;
};
 
#ifdef HAVE_NAMESPACES
}
#endif

#include <ppghmm++/close_code.h>

#endif /* _GHMM_SEQUENCES_H */

