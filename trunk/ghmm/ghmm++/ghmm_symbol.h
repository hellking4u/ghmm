/*
  author: Achim Gaedke
  created: 18. Juli 2001
  file: xmlio/examples/ghmm++/ghmm_symbol.h
  $Id$
 */

#ifndef GHMMPP_GHMM_SYMBOL_H
#define GHMMPP_GHMM_SYMBOL_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ghmm++/ghmm_symbol.h>
#include <xmlio/XMLIO_StringElement.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

/** */
class ghmm_symbol: public XMLIO_StringElement
{
public:
  /** */
  ghmm_symbol();
  /** */
  ghmm_symbol(int code, const string& symbol);
  /** */
  ghmm_symbol(const string& tag, XMLIO_Attributes &attributes);

  /** */
  virtual const char* toString() const;
  /** */
  virtual void print() const;
  /** */
  int code;
 private:
  /** */
  void init_members();
  /** */
  int valid;
};

#ifdef HAVE_NAMESPACES
}
#endif

#endif /* GHMMPP_GHMM_SYMBOL_H */
