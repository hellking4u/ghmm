/*
 * created: 21 Feb 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 */

#ifndef _GHMM_SEQUENCE_H
#define _GHMM_SEQUENCE_H 1

#include <xmlio/XMLIO_ArrayElement.h>
#include "ghmm++/GHMM_Sequences.h"

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Sequence;
class GHMM_Alphabet;

/** */
class GHMM_Sequence: public XMLIO_ArrayElement<string> {

 public:

  /** Constructor. */
  GHMM_Sequence(GHMM_SequenceType my_sequence_type, int len, int weight = 1);
  /** Constructor. */
  GHMM_Sequence(GHMM_Alphabet* my_alphabet, int len, int weight = 1);
  /** Destructor. */
  virtual ~GHMM_Sequence();

  /** Returns name of class. */
  virtual const char* toString() const;

  /** */
  void resize(int new_len);
  /** */
  virtual void XMLIO_finishedReading();

  /**
     Prints one array of sequences in a file.
     @param file       output file
     @param discrete   switch: 0 means double output for symbols,  
                       1 means truncate symbols to integer
  */
  void print(FILE *file, int discrete = 0) const;
  /** Sets double value at given position in sequence. */
  void setDouble(int index, double value);
  /** Sets integer value at given position in sequence. */
  void setInt(int index, int value);

  /** Type of current sequence. */
  GHMM_SequenceType sequence_type;

  /** C type int sequence. */
  sequence_t*   c_i_sequences;
  /** C type double sequence. */
  sequence_d_t* c_d_sequences;
  /** */
  int weight;

  
 private:

  /** */
  void init();
  /** */
  void init_INT(GHMM_Alphabet* my_alphabet, int len, int my_weight);
  /** */
  void init_DOUBLE(int len, int my_weight);

  /** Alphabet for discrete sequences .*/
  GHMM_Alphabet* alphabet;
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_SEQUENCE_H */
