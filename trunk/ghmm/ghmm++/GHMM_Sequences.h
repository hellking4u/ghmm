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
#include "ghmm++/GHMM_Types.h"

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Sequences;
class GHMM_Sequence;
class GHMM_Alphabet;

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
  /** Constructor. */
  GHMM_Sequences(GHMM_Alphabet* my_alphabet);
  /** Constructor. Creates data structure from GHMM_Sequence object. */
  GHMM_Sequences(GHMM_Sequence* sequence);
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
  int add(const string& sequence);
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
  /** Returns number of sequences. */
  unsigned int getNumberOfSequences() const;
  /** Returns sequence as string. */
  string getSequence(int index) const;
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
  /** Generates all possible integer sequence of lenght n from an alphabet with
      M letters. Use lexicographical ordering. Memory allocation here.
      @param n     length of sequences
      @param M     size of alphabet
  */
  void lexWords(int n, int M);
  /** 
      Determines max length of all sequences.
      @return max sequence length
  */
  int max_len() const;
  /**
     Prints one array of sequences in a file.
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

  /** Type of current sequence. */
  GHMM_SequenceType sequence_type;
  /** C type int sequence. */
  sequence_t*   c_i_sequences;
  /** C type double sequence. */
  sequence_d_t* c_d_sequences;
  /** */
  int last_weight;


 protected:

  /** */
  virtual void XMLIO_finishedReading();
  /** Called by GHMM_Document when a start tag is received. Tag and 
      attributes are passed to this function. */
  virtual XMLIO_Element* XMLIO_startTag(const string& my_tag, XMLIO_Attributes& my_attributes);
  /** Writes the content (XML Spec[43]) of this element.
      You should use the public XMLIO_Document::write* functions.
      @return Returns the number of bytes written,
      but is negative when an error occured and 0 by default. */
  virtual const int XMLIO_writeContent(XMLIO_Document& doc);


 private:

  /** Used by constructor. */
  void init();
  /** Used by constructor. */
  void init_INT(GHMM_Alphabet* my_alphabet, sequence_t* my_c_i_sequences);
  /** Used by constructor. */
  void init_DOUBLE(sequence_d_t* my_c_d_sequences);

  /** */
  vector<GHMM_Sequence*> sequences;
  /** alphabet used */
  GHMM_Alphabet* alphabet;
  /** Does object own this alphabet? */
  bool own_alphabet;
};
 
#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_SEQUENCES_H */
