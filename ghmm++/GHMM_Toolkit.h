/*
 * created: 06 Mar 2002 by Peter Pipenbacher
 * authors: Peter Pipenbacher (pipenb@zpr.uni-koeln.de)
 * file   : $Source$
 * $Id$
 */

#ifndef _GHMM_TOOLKIT_H
#define _GHMM_TOOLKIT_H 1

#include <ghmm++/begin_code.h>

#ifdef HAVE_NAMESPACES
namespace std {
#endif

class GHMM_Toolkit;

/** TODO: Documentation */
class GHMM_Toolkit {

 public:

  static void gsl_rng_init();
};

#ifdef HAVE_NAMESPACES
}
#endif

#include <ghmm++/close_code.h>

#endif /* _GHMM_TOOLKIT_H */
