/*
  author: Achim Gaedke
  created: 18. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_alphabet.h
  $Id$
 */

#ifndef GHMMPP_GHMM_ALPHABET_H
#define GHMMPP_GHMM_ALPHABET_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ghmm++/ghmm_symbol.h>
#include <xmlio/XMLIO_ArrayElement.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

/** */
class ghmm_alphabet: public XMLIO_ElementArrayElement<ghmm_symbol>
{
public:
  /** */
  ghmm_alphabet();
  /** */
  ghmm_alphabet(const string& tag, XMLIO_Attributes &attributes);

  /** */
  virtual const char* toString() const;
  /** */
  virtual void print() const;
  /**
   */
  virtual const XMLIO_Attributes& XMLIO_getAttributes() const;
  /**
     the id of this alphabet for reference
   */
  string id;
 private:
  void init_members();
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_ALPHABET_H */
