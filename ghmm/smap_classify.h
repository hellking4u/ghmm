/*******************************************************************************
  author       : Bernd Wichern
  filename     : ghmm/ghmm/smap_classify.h
  created      : TIME: 13:40:31     DATE: Wed 12. January 2000
  $Id$

__copyright__

*******************************************************************************/

#ifndef SMAP_CLASSIFY_H
#define SMAP_CLASSIFY_H

#include <ghmm/smodel.h>

/**@name smap functions */
/*@{ */

///
int smap_classify(smodel **smo, double *result, int smo_number, 
		   double *O, int T);

///
int smap_bayes(smodel **smo, double *result, int smo_number, double *O, int T);

/*@} */

#endif /* SMAP_CLASSIFY_H */
