/*
 * created: 06 Mar 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 *
 * __copyright__
 *
 */

#ifndef _GHMM_TOOLKIT_H
#define _GHMM_TOOLKIT_H 1

#include <string>

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Toolkit;

/** TODO: Documentation */
class GHMM_Toolkit {

 public:

  /** */
  static void ghmm_rng_init();
  /** Converts integer to string. */
  static string toString(int var);

  /** Converts integer to string. */
  static string toString(double var);


};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_TOOLKIT_H */
